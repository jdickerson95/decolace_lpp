#!Python
# ===================================================================
# ScriptName     decolace_LPP.acquire
# Purpose:       Per-hole two-shot Record: defocus guess, then ronchi,
#                science image, then +1 um CTF probe (_defocus_img).
# ===================================================================

import os

import numpy as np

import PACEtomo_ctf_calibrations as ctfcal
import PACEtomo_lafis as lafis
import PACEtomo_ronchi as ronchi


def _echo(log, text):
    if log is not None:
        log(text)
    else:
        print(text)


def _frame_saving_path(sem):
    path = sem.ReportFrameSavingPath()
    if isinstance(path, (list, tuple)):
        path = path[0] if path else "NONE"
    return str(path)


def _set_frame_folder(sem, folder_name, log=None):
    original = _frame_saving_path(sem)
    if original.upper() == "NONE":
        return None
    target = os.path.join(original, folder_name)
    try:
        os.makedirs(target, exist_ok=True)
    except Exception:
        pass
    try:
        sem.SetFolderForFrames(target)
        _echo(log, f"SetFolderForFrames -> {target}")
        return original
    except Exception:
        try:
            sem.SetFolderForFrames(folder_name)
            _echo(log, f"SetFolderForFrames -> {folder_name}")
            return original
        except Exception:
            _echo(log, "WARNING: Could not set defocus_img frame folder.")
            return None


def _restore_frame_folder(sem, original, log=None):
    if original is None:
        return
    try:
        sem.SetFolderForFrames(original)
        _echo(log, f"Restored frame folder -> {original}")
    except Exception:
        _echo(log, f"WARNING: Could not restore frame folder to {original}")


def set_frame_basename(sem, name):
    sem.SetFrameNameFormat(0, 0, 0x40)
    sem.SetFrameNameFormat(0, 1, 0x400)
    sem.SetFrameBaseName(0, 1, 0, name)


def apply_predicted_defocus(focus_base, delta_z, target_defocus, log=None):
    """SetDefocus so science CTF aims at target_defocus. Uses defocus-error slope via ChangeFocus from current."""
    desired = float(focus_base) + float(delta_z) + float(target_defocus)
    current = float(ctfcal._sem.ReportDefocus())
    commanded = ctfcal.change_focus_for_desired_delta(desired - current)
    _echo(log, f"Defocus guess: Set/Change toward {desired:.4f} um (commanded delta {commanded:.4f} um)")
    return desired, commanded


def run_ctffind(sem, defocus_lo, defocus_hi, working_xt, acquire_xt, back_project=True, log=None, label="CTF"):
    """CtfFind buffer A; optionally back-project to working X-tilt."""
    sem.NoMessageBoxOnError(1)
    try:
        cfind = sem.CtfFind("A", float(defocus_lo), float(defocus_hi))
    finally:
        sem.NoMessageBoxOnError(0)
    parsed = ctfcal.parse_ctffind(cfind)
    ax, ay = ctfcal.astig_components(parsed["astig_um"], parsed["astig_angle_deg"])
    if back_project:
        d_w, ax_w, ay_w = ctfcal.correct_ctf_to_working_xtilt(
            parsed["defocus_um"], ax, ay, working_xt, acquire_xt
        )
    else:
        d_w, ax_w, ay_w = parsed["defocus_um"], ax, ay
    mag, ang = ctfcal.astig_from_components(ax_w, ay_w)
    ctf = dict(parsed)
    ctf["defocus_raw_um"] = parsed["defocus_um"]
    ctf["defocus_um"] = d_w
    ctf["astig_x_um"] = ax_w
    ctf["astig_y_um"] = ay_w
    if np.isfinite(mag):
        ctf["astig_um"] = mag
        ctf["astig_angle_deg"] = ang
    raw = ctf.get("defocus_raw_um", float("nan"))
    res = ctf.get("resolution_A", float("nan"))
    _echo(
        log,
        f"{label} CTF {raw:.4f} -> working {d_w:.4f} um ({res:.1f} A)"
        if np.isfinite(d_w)
        else f"WARNING: {label} CtfFind failed.",
    )
    return ctf


def maybe_correct_astig(ctf, hole_count, correct_astig, every_n, log=None):
    if not correct_astig or every_n <= 0:
        return False
    if hole_count % int(every_n) != 0:
        return False
    if not ctfcal.astig_calibration_file or not ctfcal.get_astig_calibration():
        _echo(log, "WARNING: correctAstig skipped (no astig calibration).")
        return False
    if ctfcal.should_use_ctf_xtilt() and not ctfcal.get_xtilt_calibration():
        _echo(log, "WARNING: correctAstig skipped (no xtilt_defoc_astig calibration).")
        return False
    ax, ay = ctf.get("astig_x_um", np.nan), ctf.get("astig_y_um", np.nan)
    if not np.isfinite(ax) or not np.isfinite(ay):
        _echo(log, "WARNING: correctAstig skipped (invalid CTF astig).")
        return False
    dsx, dsy = ctfcal.stig_delta_to_cancel_astig(ax, ay)
    if not np.isfinite(dsx) or not np.isfinite(dsy):
        _echo(log, "WARNING: correctAstig skipped (invalid stig delta).")
        return False
    sx, sy, *_ = ctfcal._sem.ReportObjectiveStigmator()
    new_sx = float(np.clip(float(sx) + dsx, -1.0, 1.0))
    new_sy = float(np.clip(float(sy) + dsy, -1.0, 1.0))
    ctfcal._sem.SetObjectiveStigmator(new_sx, new_sy)
    _echo(
        log,
        f"Astig realign after hole {hole_count}: working astig "
        f"({ax:.4f}, {ay:.4f}) um -> stig ({new_sx:.6f}, {new_sy:.6f})",
    )
    return True


def acquire_hole(
    sem,
    hole,
    stem,
    directory,
    target_defocus,
    ctf_probe_underfocus,
    focus_base,
    delta_z,
    delay_is,
    ctf_defocus_lo,
    ctf_defocus_hi,
    count_threshold=0,
    beam_tilt_comp=True,
    hole_count=1,
    correct_astig=True,
    astig_every_n=10,
    test_mode=False,
    test_science_ctf_xtilt=True,
    log=None,
):
    """ImageShift, defocus guess, ronchi, science Record, then CTF probe."""
    idx = int(hole["index"])
    ssx, ssy = float(hole["ssx"]), float(hole["ssy"])
    science_name = f"{stem}_{idx:04d}"
    probe_name = f"{stem}_{idx:04d}_defocus_img"
    science_mrc = os.path.join(directory, science_name + ".mrc")
    probe_mrc = os.path.join(directory, probe_name + ".mrc")

    sem.SetImageShift(0, 0)
    sem.ImageShiftByMicrons(ssx, ssy)
    apply_predicted_defocus(focus_base, delta_z, target_defocus, log=log)
    is_x, is_y, *_ = sem.ReportImageShift()
    is_mag = float(np.hypot(float(is_x), float(is_y)))
    if beam_tilt_comp:
        lafis.doLafis(is_x, is_y)

    set_frame_basename(sem, science_name)
    ronchi.doRonchigramCorrection(pos=idx, pn=0)
    sem.Delay(float(delay_is), "s")

    working_xt = ctfcal.report_xtilt()
    ctf_xt = working_xt
    at_ctf_xtilt = False
    if test_mode and test_science_ctf_xtilt and ctfcal.should_use_ctf_xtilt():
        ctf_xt = ctfcal.ctf_xtilt_target()
        ctfcal.set_xtilt(ctf_xt[0], ctf_xt[1])
        at_ctf_xtilt = True
        _echo(log, f"test_mode: science at CTF X-tilt ({ctf_xt[0]:.6f}, {ctf_xt[1]:.6f})")

    sem.OpenNewFile(science_mrc)
    try:
        ronchi.recordWithRonchi(run_ronchi=False, acquire_label=f"science hole {idx}", pos=idx, pn=0)
    finally:
        sem.CloseFile()

    science_ctf = None
    if test_mode:
        science_ctf = run_ctffind(
            sem,
            ctf_defocus_lo,
            ctf_defocus_hi,
            working_xt,
            ctf_xt if at_ctf_xtilt else working_xt,
            back_project=at_ctf_xtilt,
            log=log,
            label=f"[{idx}] science",
        )

    if count_threshold > 0:
        counts = sem.ReportMeanCounts()
        if counts > count_threshold:
            sem.CenterBeamFromImage(0, 0.4)

    commanded = ctfcal.change_focus_for_desired_delta(-float(ctf_probe_underfocus))
    _echo(log, f"Probe underfocus commanded {commanded:.4f} um (want {ctf_probe_underfocus} um physical)")
    probe_nominal = float(sem.ReportDefocus())
    frame_folder = _set_frame_folder(sem, "defocus_img", log=log)
    set_frame_basename(sem, probe_name)
    ctf = None
    try:
        if not at_ctf_xtilt and ctfcal.should_use_ctf_xtilt():
            ctf_xt = ctfcal.ctf_xtilt_target()
            ctfcal.set_xtilt(ctf_xt[0], ctf_xt[1])
            at_ctf_xtilt = True
        elif not ctfcal.should_use_ctf_xtilt():
            ctf_xt = working_xt
        sem.OpenNewFile(probe_mrc)
        try:
            ronchi.recordWithRonchi(
                run_ronchi=False,
                acquire_label=f"defocus_img hole {idx}",
                pos=idx,
                pn=0,
            )
        finally:
            sem.CloseFile()
        ctf = run_ctffind(
            sem,
            ctf_defocus_lo,
            ctf_defocus_hi,
            working_xt,
            ctf_xt,
            back_project=True,
            log=log,
            label=f"[{idx}] probe",
        )
    finally:
        if at_ctf_xtilt:
            ctfcal.set_xtilt(working_xt[0], working_xt[1])
        _restore_frame_folder(sem, frame_folder, log=log)
        ctfcal.change_focus_for_desired_delta(float(ctf_probe_underfocus))
        if beam_tilt_comp:
            lafis.restoreLafis()
        sem.SetImageShift(0, 0)

    maybe_correct_astig(ctf or {}, hole_count, correct_astig, astig_every_n, log=log)
    return {
        "science_mrc": science_mrc,
        "probe_mrc": probe_mrc,
        "probe_nominal_defocus": probe_nominal,
        "ctf": ctf,
        "science_ctf": science_ctf,
        "is_x": float(is_x),
        "is_y": float(is_y),
        "is_mag": is_mag,
    }
