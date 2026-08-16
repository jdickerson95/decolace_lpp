#!Python
# ===================================================================
# ScriptName     PACEtomo_ronchi
# Purpose:       Ronchigram Trial acquire, FFT analysis, C3 and laser
#                X-tilt correction. Shared by PACEtomo and decolace_LPP.
# ===================================================================

import numpy as np
import serialem as sem

import PACEtomo_lafis as lafis

_sem = sem
_logger = None
_debug = False
_delay_is = 0.1
_c3_position = None  # callable -> PACEtomo position list, or the list itself
_beam_tilt_comp = True

hasXLens = True
doRonchigram = True
ronchiBaseSuffix = "_ronchi"
ronchiC3Offset = -20
ronchiDelay = 1.0
ronchiBinning = 32
ronchiPixelSize = 0.98e-4 * 2
ronchiTargetPhaseA = -1.93941993
ronchiTargetPhaseB = 1.67658165
ronchiCorrectKs = [[9.303, -0.662], [0.856, 8.680]]
ronchiPeakRadius = 100
ronchiMontage = True
ronchiCorrMatrix = [[0.212, 1.28], [1.22, -0.243]]
ronchiCorrectC3 = True
ronchiC3CorrectionFactor = 20 / 9.1
ronchiMinErrForC3Correction = 0.3
ronchiMinErrForC3CorrectionRedo = 0.5
redo_ronchi_after_C3 = True
ronchiPerPositionC3 = True
ronchiXLensTolerance = 0.000125
ronchiStartXLensX = None
ronchiStartXLensY = None
ronchiStartC3Offset = None


def _echo(text):
    if _logger is not None:
        _logger(text)
    else:
        _sem.Echo(text)


def _position():
    if _c3_position is None:
        return None
    if callable(_c3_position):
        return _c3_position()
    return _c3_position


def configure(
    sem_module=None,
    logger=None,
    has_x_lens=None,
    do_ronchigram=None,
    debug=None,
    delay_is=None,
    c3_position_store=None,
    beam_tilt_comp=None,
    **kwargs,
):
    """kwargs: ronchi* setting names matching PACEtomo.py."""
    global _sem, _logger, _debug, _delay_is, _c3_position, _beam_tilt_comp
    global hasXLens, doRonchigram
    if sem_module is not None:
        _sem = sem_module
    if logger is not None:
        _logger = logger
    if has_x_lens is not None:
        hasXLens = bool(has_x_lens)
    if do_ronchigram is not None:
        doRonchigram = bool(do_ronchigram)
    if debug is not None:
        _debug = bool(debug)
    if delay_is is not None:
        _delay_is = float(delay_is)
    if c3_position_store is not None:
        _c3_position = c3_position_store
    if beam_tilt_comp is not None:
        _beam_tilt_comp = bool(beam_tilt_comp)
    g = globals()
    for key, value in kwargs.items():
        if key in g:
            g[key] = value


def _ronchi_bin_image(image, binning=32):
    bins = [
        image[i:(image.shape[0] // binning * binning):binning, j:(image.shape[1] // binning * binning):binning]
        for i in range(binning)
        for j in range(binning)
    ]
    return np.sum(bins, axis=0)


def _ronchi_find_fourier_centered(image, padded_size=4096):
    fourier = np.fft.fftshift(np.fft.fft2(image - np.mean(image), s=(padded_size, padded_size)))
    center = np.array(image.shape) / 2
    grid_y, grid_x = np.indices((padded_size, padded_size), dtype=float)
    shifted_y = grid_y - padded_size / 2
    shifted_x = grid_x - padded_size / 2
    correction_phase_x = np.exp(2j * np.pi * shifted_x / padded_size * center[1])
    correction_phase_y = np.exp(2j * np.pi * shifted_y / padded_size * center[0])
    return fourier * correction_phase_x * correction_phase_y


def _ronchi_report_angles(ks, start_angle=-135):
    return np.mod(np.arctan2(-ks[:, 0], ks[:, 1]) * 180 / np.pi - start_angle, 360) + start_angle


def _ronchi_find_peaks(fourier, radius=100, npeaks=4):
    fourier_abs = np.abs(fourier).copy()
    peak_locations = []
    phases = []
    size = np.shape(fourier_abs)[0]
    peak_coords = [size // 2, size // 2]
    fourier_abs[
        (peak_coords[0] - radius):(peak_coords[0] + radius),
        (peak_coords[1] - radius):(peak_coords[1] + radius),
    ] = 0
    for _ in range(npeaks):
        max_idx = np.argmax(fourier_abs)
        peak_coords = np.unravel_index(max_idx, fourier_abs.shape)
        peak_locations.append(peak_coords)
        fourier_abs[
            (peak_coords[0] - radius):(peak_coords[0] + radius),
            (peak_coords[1] - radius):(peak_coords[1] + radius),
        ] = 0
        phases.append(np.angle(fourier[peak_coords]))
    return np.array(peak_locations) - size / 2, np.array(phases)


def _ronchi_find_ks_phases(corrected_fourier, pixel_size_um, npeaks=2, radius=100, binning=1, fourier_size=None):
    if fourier_size is None:
        fourier_size = corrected_fourier.shape[0]
    peaks, phases = _ronchi_find_peaks(corrected_fourier, radius=radius, npeaks=npeaks * 2)
    ordering = np.argsort(_ronchi_report_angles(peaks, start_angle=-135))
    peaks = peaks[ordering][:npeaks]
    phases = phases[ordering][:npeaks]
    ks = peaks / fourier_size * 1 / (pixel_size_um * binning)
    return ks, phases


def analyze_ronchigram(
    image,
    pixel_size_um,
    binning,
    target_phase_a,
    target_phase_b,
    correct_ks,
    peak_radius=100,
    corr_matrix=None,
    corr_scale=1e-5,
    c3_correction_factor=20 / 6.85,
):
    """FFT peak phases -> laser deflector and C3 corrections."""
    if corr_matrix is None:
        corr_matrix = [[0.212, 1.28], [1.22, -0.243]]
    correct_ks = np.asarray(correct_ks, dtype=float)
    binned = _ronchi_bin_image(np.asarray(image), binning=binning)
    image_fft = _ronchi_find_fourier_centered(binned)
    ks, phases = _ronchi_find_ks_phases(
        image_fft,
        pixel_size_um * binning,
        npeaks=2,
        radius=peak_radius,
        binning=1,
        fourier_size=image_fft.shape[0],
    )
    ks_error = ks - correct_ks
    ks_total_err = float(np.linalg.norm(ks_error))
    ks_avg_err = (ks_error[0, 0] + ks_error[1, 1]) / 2.0
    c3_correction = ks_avg_err * c3_correction_factor
    phase_err_a = np.mod(phases[0] - target_phase_a + np.pi, 2 * np.pi) - np.pi
    phase_err_b = np.mod(phases[1] - target_phase_b + np.pi, 2 * np.pi) - np.pi
    corr = np.asarray(corr_matrix, dtype=float) * corr_scale
    correction_x = phase_err_a * corr[0, 0] + phase_err_b * corr[0, 1]
    correction_y = phase_err_a * corr[1, 0] + phase_err_b * corr[1, 1]
    return {
        "ks": ks,
        "phases": phases,
        "ks_error": ks_error,
        "ks_total_err": ks_total_err,
        "ks_avg_err": ks_avg_err,
        "c3_correction": c3_correction,
        "phase_err_a": phase_err_a,
        "phase_err_b": phase_err_b,
        "correction_x": correction_x,
        "correction_y": correction_y,
    }


def applyRonchigramXtiltCorrection(correction_x, correction_y, lens_index=2):
    if not hasXLens:
        return
    xtX, xtY = _sem.ReportXLensDeflector(lens_index)
    _sem.SetXLensDeflector(lens_index, xtX + correction_x, xtY + correction_y)


def _reset_ronchi_xlens_if_out_of_tolerance(lens_index=2):
    global ronchiStartXLensX, ronchiStartXLensY
    if not hasXLens or ronchiStartXLensX is None or ronchiStartXLensY is None:
        return
    xtX, xtY = _sem.ReportXLensDeflector(lens_index)
    xtX, xtY = float(xtX), float(xtY)
    xtX_delta = xtX - lafis.lafisXtCorrectionX - ronchiStartXLensX
    xtY_delta = xtY - lafis.lafisXtCorrectionY - ronchiStartXLensY
    if abs(xtX_delta) > ronchiXLensTolerance or abs(xtY_delta) > ronchiXLensTolerance:
        _echo(
            f"WARNING: Ronchigram X lens deflector ({xtX:.6f}, {xtY:.6f}) beyond tolerance "
            f"{ronchiXLensTolerance} from start ({ronchiStartXLensX:.6f}, {ronchiStartXLensY:.6f}); resetting."
        )
        _sem.SetXLensDeflector(
            lens_index,
            ronchiStartXLensX + lafis.lafisXtCorrectionX,
            ronchiStartXLensY + lafis.lafisXtCorrectionY,
        )


def applyRonchigramC3Correction(c3_correction, baseline_offset):
    new_offset = baseline_offset + c3_correction
    _sem.SetImageDistanceOffset(new_offset)
    return new_offset


def _report_frame_basename():
    r = _sem.ReportFrameBaseName()
    use_in_frame = int(r[0])
    name = r[1] if len(r) > 1 else ""
    use_in_folder = int(r[2]) if len(r) > 2 else 0
    if isinstance(name, str) and name.lower() == "none":
        name = ""
    return use_in_frame, name, use_in_folder


def _ronchi_trial_basename(use_in_frame, current_name, use_in_folder):
    root = (current_name or "").strip()
    return root + ronchiBaseSuffix, root


def _set_ronchi_trial_frame_basename():
    use_in_frame, name, use_in_folder = _report_frame_basename()
    ronchi_name, root = _ronchi_trial_basename(use_in_frame, name, use_in_folder)
    if not root:
        _echo("WARNING: Ronchigram Trial: no frame base name; SetFrameBaseName before acquire.")
    else:
        _sem.SetFrameBaseName(0, use_in_frame, use_in_folder, ronchi_name)
        _echo(f"Ronchigram Trial: SetFrameBaseName -> {ronchi_name}")
    return use_in_frame, name, use_in_folder


def _restore_frame_basename(saved):
    use_in_frame, name, use_in_folder = saved
    _sem.SetFrameBaseName(0, use_in_frame, use_in_folder, name)
    _echo(f"Ronchigram Trial: restored frame base name -> {name or '(none)'}")


def checkRonchigramSetup():
    global ronchiStartXLensX, ronchiStartXLensY, ronchiStartC3Offset, doRonchigram
    if not hasXLens and doRonchigram:
        _sem.OKBox(
            "ERROR: doRonchigram requires XLensDeflector. "
            "Set hasXLens = True, or set doRonchigram = False."
        )
        _sem.Exit()
    if not doRonchigram:
        if not hasXLens:
            _echo("NOTE: hasXLens=False; all XLens Report/Set/Restore calls are skipped.")
        return
    if ronchiStartXLensX is None:
        ronchiStartXLensX, ronchiStartXLensY = [float(v) for v in _sem.ReportXLensDeflector(2)[:2]]
    if ronchiStartC3Offset is None:
        ronchiStartC3Offset = float(_sem.ReportImageDistanceOffset())
    _echo(
        f"NOTE: Ronchigram X lens deflector start ({ronchiStartXLensX:.6f}, {ronchiStartXLensY:.6f}), "
        f"reset tolerance {ronchiXLensTolerance}"
    )
    _echo(f"NOTE: Ronchigram C3 (ImageDistanceOffset) session start {ronchiStartC3Offset:.2f} um")
    trial_exp, *_ = _sem.ReportExposure("T")
    record_exp, *_ = _sem.ReportExposure("R")
    if trial_exp <= 0 and _sem.IsVariableDefined("warningRonchiTrial") == 0:
        _sem.Pause(
            "WARNING: Trial exposure is zero or not set. Configure Trial with a very short exposure at the same position as Record."
        )
        _sem.SetPersistentVar("warningRonchiTrial", "")
    if trial_exp >= record_exp * 0.5:
        _echo("WARNING: Trial exposure should be much shorter than Record for negligible ronchigram dose.")
    try:
        t_shift = _sem.ReportLDAreaShift("T")
        if len(t_shift) >= 2 and np.linalg.norm(np.array(t_shift[:2], dtype=float)) > 0.01:
            _echo("WARNING: Trial area position differs from Record. Set Trial LD offsets to match Record.")
    except (AttributeError, TypeError, ValueError):
        pass
    _echo(
        "NOTE: Ronchigram uses Trial at Record beam position. Set Trial LD offsets identical to Record; only exposure should differ."
    )
    _echo(f"NOTE: Ronchigram Trial temporarily appends '{ronchiBaseSuffix}' to the active frame base name.")
    if ronchiCorrectC3:
        _echo(
            f"NOTE: Ronchigram C3 correction enabled (factor {ronchiC3CorrectionFactor:.4f}, "
            f"only if |C3 correction| > {ronchiMinErrForC3Correction} um)."
        )
    if redo_ronchi_after_C3:
        _echo(
            "NOTE: Ronchigram may use up to 3 Trials before Record: "
            "1st C3 (if above threshold), 2nd phase or C3, 3rd phase-only if 2nd C3 applied."
        )
        _echo(
            f"NOTE: 2nd-Trial C3 threshold |correction| > {ronchiMinErrForC3CorrectionRedo} um "
            f"(1st-Trial threshold {ronchiMinErrForC3Correction} um)."
        )
    if ronchiPerPositionC3:
        _echo("NOTE: Ronchigram C3 offset is stored and restored per target across tilts.")


def reset_to_session_start(when=""):
    """Restore XLens(2) and C3 to values captured in checkRonchigramSetup."""
    if not doRonchigram:
        return
    suffix = f" {when}" if when else ""
    if hasXLens and ronchiStartXLensX is not None and ronchiStartXLensY is not None:
        _sem.SetXLensDeflector(2, ronchiStartXLensX, ronchiStartXLensY)
        _echo(
            f"NOTE: Reset XLensDeflector(2) to session start "
            f"({ronchiStartXLensX:.6f}, {ronchiStartXLensY:.6f}){suffix}."
        )
    if ronchiStartC3Offset is not None:
        _sem.SetImageDistanceOffset(ronchiStartC3Offset)
        _echo(f"NOTE: Reset ImageDistanceOffset to session start {ronchiStartC3Offset:.2f} um{suffix}.")


def _apply_stored_c3_offset(pos, pn):
    position = _position()
    if position is None:
        return float(_sem.ReportImageDistanceOffset())
    stored = position[pos][pn].get("c3_offset")
    if stored is not None:
        _sem.SetImageDistanceOffset(stored)
        _echo(f"Ronchigram: target {pos + 1} using stored C3 offset {stored:.2f} um")
        return float(stored)
    current = float(_sem.ReportImageDistanceOffset())
    for b in range(len(position[pos])):
        position[pos][b]["c3_offset"] = current
    _echo(f"Ronchigram: target {pos + 1} initializing C3 offset {current:.2f} um")
    return current


def _save_c3_offset_for_target(pos, pn):
    position = _position()
    if position is None:
        return
    current = float(_sem.ReportImageDistanceOffset())
    for b in range(len(position[pos])):
        position[pos][b]["c3_offset"] = current
    _echo(f"Ronchigram: target {pos + 1} saved C3 offset {current:.2f} um")


def _log_ronchi_ks(result, pass_label=""):
    prefix = f"Ronchigram{pass_label}"
    _echo(
        f"{prefix} ks (1/um): {np.array2string(result['ks'], precision=4)} | "
        f"ks error: {np.array2string(result['ks_error'], precision=4)}"
    )
    _echo(
        f"{prefix} ||ks error||: {result['ks_total_err']:.4f} (1/um) | "
        f"mean diagonal ks error: {result['ks_avg_err']:.4f} (1/um) | "
        f"recommended C3 correction: {result['c3_correction']:.2f} um"
    )


def _log_ronchi_phases(result, pass_label=""):
    phases = result["phases"]
    prefix = f"Ronchigram{pass_label}"
    _echo(
        f"{prefix} phases (rad): measured vertical={phases[0]:.3f} horizontal={phases[1]:.3f} | "
        f"targets vertical={ronchiTargetPhaseA:.3f} horizontal={ronchiTargetPhaseB:.3f}"
    )
    _echo(
        f"{prefix} phase error (rad): vertical={result['phase_err_a']:.3f} "
        f"horizontal={result['phase_err_b']:.3f} | "
        f"deflector dX={result['correction_x']:.3e} dY={result['correction_y']:.3e}"
    )


def _analyze_ronchi_image(image):
    return analyze_ronchigram(
        image,
        ronchiPixelSize,
        ronchiBinning,
        ronchiTargetPhaseA,
        ronchiTargetPhaseB,
        ronchiCorrectKs,
        peak_radius=ronchiPeakRadius,
        corr_matrix=ronchiCorrMatrix,
        c3_correction_factor=ronchiC3CorrectionFactor,
    )


def _acquire_ronchi_trial(trial_offset_baseline, pass_label=""):
    _reset_ronchi_xlens_if_out_of_tolerance()
    is_x, is_y, *_ = _sem.ReportImageShift()
    _sem.GoToLowDoseArea("T")
    _sem.SetImageShift(0, 0)
    _sem.SetImageShift(is_x, is_y)
    _sem.SetImageDistanceOffset(trial_offset_baseline + ronchiC3Offset)
    saved_basename = _set_ronchi_trial_frame_basename()
    try:
        _sem.Delay(ronchiDelay, "s")
        lafis.add_lpp_meta_to_next_mdoc()
        _sem.T()
    finally:
        _sem.SetImageDistanceOffset(trial_offset_baseline)
        _restore_frame_basename(saved_basename)
    if pass_label:
        _echo(f"Ronchigram{pass_label}: Trial image acquired.")
    return np.asarray(_sem.bufferImage("A"))


def _try_apply_ronchi_c3(result, c3_baseline_offset, pass_label="", min_err=None):
    if min_err is None:
        min_err = ronchiMinErrForC3Correction
    prefix = f"Ronchigram{pass_label}"
    if not ronchiCorrectC3:
        _echo(f"{prefix} C3: correction {result['c3_correction']:.2f} um not applied (ronchiCorrectC3=False)")
        return False
    if abs(result["c3_correction"]) <= min_err:
        _echo(
            f"{prefix} C3: skipped (|correction| {abs(result['c3_correction']):.2f} um <= "
            f"minimum {min_err} um)"
        )
        return False
    new_offset = applyRonchigramC3Correction(result["c3_correction"], c3_baseline_offset)
    _echo(
        f"{prefix} C3: ImageDistanceOffset adjusted by {result['c3_correction']:.2f} um "
        f"(now {new_offset:.2f} um)"
    )
    return True


def _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=""):
    image = _acquire_ronchi_trial(c3_baseline_offset, pass_label=pass_label)
    result = _analyze_ronchi_image(image)
    _log_ronchi_ks(result, pass_label=pass_label)
    return result


def _apply_ronchi_phase(result, pass_label=""):
    _log_ronchi_phases(result, pass_label=pass_label)
    applyRonchigramXtiltCorrection(result["correction_x"], result["correction_y"])


def doRonchigramCorrection(set_track_fn=None, pos=None, pn=None):
    if not doRonchigram:
        return
    try:
        _sem.UpdateLowDoseParams("T")
    except AttributeError:
        pass
    tilt = float(_sem.ReportTiltAngle())
    trial_exp, *_ = _sem.ReportExposure("T")
    _echo(
        f"Ronchigram: Trial acquire at tilt {tilt:.1f} deg | "
        f"C3 offset {ronchiC3Offset} um | binning {ronchiBinning} | Trial exposure {trial_exp:.4g} s"
    )
    if ronchiPerPositionC3 and pos is not None and pn is not None:
        c3_baseline_offset = _apply_stored_c3_offset(pos, pn)
    else:
        c3_baseline_offset = float(_sem.ReportImageDistanceOffset())
    try:
        result = _ronchi_trial_and_analyze(c3_baseline_offset)
        c3_changed = _try_apply_ronchi_c3(result, c3_baseline_offset)
        if c3_changed and redo_ronchi_after_C3:
            _echo("Ronchigram: phase correction deferred until after 2nd Trial.")
            c3_baseline_offset = float(_sem.ReportImageDistanceOffset())
            _echo("Ronchigram: 2nd Trial after 1st C3 change.")
            result = _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=" (2nd)")
            c3_changed_redo = _try_apply_ronchi_c3(
                result,
                c3_baseline_offset,
                pass_label=" (2nd)",
                min_err=ronchiMinErrForC3CorrectionRedo,
            )
            if c3_changed_redo:
                _echo("Ronchigram: phase correction deferred until after 3rd Trial.")
                c3_baseline_offset = float(_sem.ReportImageDistanceOffset())
                _echo("Ronchigram: 3rd Trial (phase correction only).")
                result = _ronchi_trial_and_analyze(c3_baseline_offset, pass_label=" (3rd)")
                _apply_ronchi_phase(result, pass_label=" (3rd)")
            else:
                _apply_ronchi_phase(result, pass_label=" (2nd)")
        else:
            _apply_ronchi_phase(result)
        if _debug:
            _echo(f"DEBUG: Ronchigram pixel size {ronchiPixelSize} um, peak radius {ronchiPeakRadius} px")
        if ronchiPerPositionC3 and pos is not None and pn is not None:
            _save_c3_offset_for_target(pos, pn)
    except Exception as e:
        _echo(f"WARNING: Ronchigram analysis failed: {e}. Continuing without laser correction.")
        _sem.SetImageDistanceOffset(ronchiStartC3Offset)
        if ronchiPerPositionC3 and pos is not None and pn is not None:
            _save_c3_offset_for_target(pos, pn)
    _sem.GoToLowDoseArea("R")
    if set_track_fn is not None:
        set_track_fn()


def ronchi_before_preview_align(acquire_label="preview alignment", pos=None, pn=None):
    if not doRonchigram:
        return
    _echo(f"NOTE: {acquire_label} - ronchigram before preview alignment image.")
    is_x, is_y, *_ = _sem.ReportImageShift()
    _sem.GoToLowDoseArea("R")
    _sem.SetImageShift(0, 0)
    _sem.SetImageShift(is_x, is_y)
    doRonchigramCorrection(set_track_fn=None, pos=pos, pn=pn)
    if _beam_tilt_comp:
        _echo("Ronchigram: We should not resotre the value here?.")


def recordWithRonchi(set_track_fn=None, run_ronchi=True, acquire_label="Record", pos=None, pn=None):
    if ronchiPerPositionC3 and pos is not None and pn is not None and not (run_ronchi and doRonchigram):
        _apply_stored_c3_offset(pos, pn)
    if run_ronchi and doRonchigram:
        _echo(f"NOTE: {acquire_label} - ronchigram Trial then Record stack frame.")
        doRonchigramCorrection(set_track_fn=set_track_fn, pos=pos, pn=pn)
    elif doRonchigram and not run_ronchi:
        _echo(f"NOTE: {acquire_label} - Record only (ronchigram skipped for this shot).")
    else:
        _echo(f"NOTE: {acquire_label} - Record acquire.")
    _sem.Delay(_delay_is, "s")
    lafis.add_lpp_meta_to_next_mdoc()
    _sem.R()
    _sem.S()
