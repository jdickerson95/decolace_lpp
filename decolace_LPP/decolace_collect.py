#!Python
# ===================================================================
# ScriptName     decolace_LPP
# Purpose:       Hexagonal 2D montage with ronchi, LAFIS, spline
#                defocus, 6-neighbor fill-in, and a CTF probe image.
# ===================================================================

############ SETTINGS ############

# Area: inclusive navigator range [start, end], or a list of 3+ clockwise items.
nav_item_range = [10, 17]
realign_nav_item = 0            # 0: skip RealignToItem; else map/center nav item (stage only)
beam_radius = 0.168             # um; hex packing radius (FOWL)
add_overlap = 0.05              # fraction; hex step = (1-overlap)*2*beam_radius*cos(30)
directory = r""                 # empty: SerialEM working directory
stem = "decolace"               # MRC / frame basename prefix
dry_run = True                  # True: build lattice, plot, write state, do not acquire

# Test mode: fittable science defocus, CtfFind both shots, write stats.
# Independent of dry_run. Use a distinct stem (e.g. decolace_test).
test_mode = False
test_target_defocus = -0.5      # um; science CTF target used only when test_mode
test_science_ctf_xtilt = True   # True: science Record at ctfXtilt; False: working X-tilt

# Defocus
target_defocus = -0.02          # um; science-image CTF target (too close to fit reliably)
ctf_probe_underfocus = 1.0      # um physical underfocus for the _defocus_img probe
geometry_file = r""             # JSON from PACEtomo_measureGeometry.py; empty = start flat
spline_resolution = (10, 10)
spline_n_iter = 800

# Calibrations (empty = skip that correction)
xtilt_calibration_file = r""    # check_xtilt_defoc_astig.py
defocus_error_file = r""        # calibrate_defocus_error.py
astig_calibration_file = r""    # calibrate_astigmatism.py

# CTF / X-tilt
hasXLens = True
useCtfXtilt = True
ctfXtiltX = 0.002836
ctfXtiltY = 0.003867
ctfDefocusLo = -12.0
ctfDefocusHi = -0.2

# Astig realign from shot-2 CTF (working X-tilt)
correctAstig = True
astigEveryN = 10                # apply after every N completed holes

# Ronchi (laser must realign every hole)
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
ronchiCorrMatrix = [[0.212, 1.28], [1.22, -0.243]]
ronchiCorrectC3 = True
ronchiC3CorrectionFactor = 20 / 9.1
ronchiMinErrForC3Correction = 0.3
ronchiMinErrForC3CorrectionRedo = 0.5
redo_ronchi_after_C3 = True
ronchiPerPositionC3 = True
ronchiXLensTolerance = 0.000125

# LAFIS
beamTiltComp = True
xt_is_matrix = [[0.000324, -0.000347], [0.001100, 0.00028125]]
df_is_matrix = [[0.041381, 0.012342], [0.041381, 0.012342]]

delayIS = 2.0
count_threshold = 1000          # CenterBeamFromImage if mean counts exceed this; 0 disables
checkDewar = True

########## END SETTINGS ##########

import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
import json
import os
from datetime import datetime

try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _here = os.getcwd()
_repo = os.path.dirname(_here)
for _p in (_here, _repo, os.getcwd()):
    if _p and _p not in sys.path:
        sys.path.append(_p)

import numpy as np
import serialem as sem

import PACEtomo_ctf_calibrations as ctfcal
import PACEtomo_geometry as pacegeo
import PACEtomo_lafis as lafis
import PACEtomo_ronchi as ronchi

try:
    from decolace_LPP import acquire as dacq
    from decolace_LPP import area as darea
    from decolace_LPP import predict as dpred
    from decolace_LPP import schedule as dsched
    from decolace_LPP import test_stats as dtest
except ImportError:
    import acquire as dacq
    import area as darea
    import predict as dpred
    import schedule as dsched
    import test_stats as dtest


def log(text):
    if text.startswith("WARNING:") or text.startswith("ERROR:"):
        if sem.IsVersionAtLeast("40200", "20240205"):
            sem.SetNextLogOutputStyle(1 if text.startswith("ERROR:") else 0, 5 if text.startswith("WARNING:") else 2)
    sem.EchoBreakLines(text)


def state_path(out_dir):
    return os.path.join(out_dir, f"{stem}_state.json")


def _json_safe(obj):
    """json.dump cannot serialize numpy scalars or NaN/Inf."""
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def _finite_or_none(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def save_state(out_dir, area, extra=None):
    payload = dict(area)
    if extra:
        payload.update(extra)
    payload["updated"] = datetime.now().isoformat(timespec="seconds")
    path = state_path(out_dir)
    with open(path, "w") as fh:
        json.dump(_json_safe(payload), fh, indent=2)
    return path


def load_state(out_dir):
    path = state_path(out_dir)
    if not os.path.exists(path):
        return None
    with open(path, "r") as fh:
        return json.load(fh)


def configure_shared():
    ctfcal.configure(
        sem_module=sem,
        logger=log,
        has_x_lens=hasXLens,
        use_ctf_xtilt=useCtfXtilt,
        ctf_xtilt_x=ctfXtiltX,
        ctf_xtilt_y=ctfXtiltY,
        xtilt_calibration_path=xtilt_calibration_file,
        defocus_error_path=defocus_error_file,
        astig_calibration_path=astig_calibration_file,
    )
    lafis.configure(
        sem_module=sem,
        logger=log,
        has_x_lens=hasXLens,
        xt_matrix=xt_is_matrix,
        df_matrix=df_is_matrix,
    )
    ronchi.configure(
        sem_module=sem,
        logger=log,
        has_x_lens=hasXLens,
        do_ronchigram=doRonchigram,
        delay_is=delayIS,
        beam_tilt_comp=beamTiltComp,
        ronchiBaseSuffix=ronchiBaseSuffix,
        ronchiC3Offset=ronchiC3Offset,
        ronchiDelay=ronchiDelay,
        ronchiBinning=ronchiBinning,
        ronchiPixelSize=ronchiPixelSize,
        ronchiTargetPhaseA=ronchiTargetPhaseA,
        ronchiTargetPhaseB=ronchiTargetPhaseB,
        ronchiCorrectKs=ronchiCorrectKs,
        ronchiPeakRadius=ronchiPeakRadius,
        ronchiCorrMatrix=ronchiCorrMatrix,
        ronchiCorrectC3=ronchiCorrectC3,
        ronchiC3CorrectionFactor=ronchiC3CorrectionFactor,
        ronchiMinErrForC3Correction=ronchiMinErrForC3Correction,
        ronchiMinErrForC3CorrectionRedo=ronchiMinErrForC3CorrectionRedo,
        redo_ronchi_after_C3=redo_ronchi_after_C3,
        ronchiPerPositionC3=ronchiPerPositionC3,
        ronchiXLensTolerance=ronchiXLensTolerance,
    )


def load_initial_geometry(area):
    if not geometry_file:
        log("NOTE: No geometry_file; starting with a flat z=0 map at the lattice centroid.")
        return None, [], []
    geo = pacegeo.load_geometry(geometry_file)
    origin_ss = None
    origin_nav = geo.get("origin_nav_item")
    if origin_nav is not None:
        origin_ss = dpred.origin_ss_from_nav(
            sem,
            origin_nav,
            area["polygon"]["reference_nav_item"],
            area["polygon"]["reference_stage"],
            log=log,
        )
    recentered = dpred.recenter_geometry(geo, area["centroid_ss_ref"], origin_ss_ref=origin_ss)
    extra_ss = list(zip(recentered["points"]["ssx"], recentered["points"]["ssy"]))
    extra_z = list(recentered["points"]["z_um"])
    log(
        f"Loaded geometry {os.path.abspath(geometry_file)} "
        f"({recentered['n_points']} points), re-centered to hex centroid."
    )
    return recentered, extra_ss, extra_z


def main():
    sem.SuppressReports()
    out_dir = directory or sem.ReportDirectory()
    os.makedirs(out_dir, exist_ok=True)
    configure_shared()
    if doRonchigram:
        ronchi.checkRonchigramSetup()

    science_target = test_target_defocus if test_mode else target_defocus
    probe_target = science_target - float(ctf_probe_underfocus)

    log("##### decolace_LPP #####")
    log(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    log(f"Output: {out_dir}")
    log("SS origin = mean of all hex holes (not the first geo / nav point).")
    if test_mode:
        log(
            f"test_mode=True: science target={science_target:.3f} um, "
            f"probe target={probe_target:.3f} um, "
            f"science_ctf_xtilt={test_science_ctf_xtilt}"
        )

    if realign_nav_item:
        log(f"RealignToItem nav {int(realign_nav_item)} (stage only)")
        sem.RealignToOtherItem(int(realign_nav_item), 0, 0, 0.05, 4, 1)

    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)
    focus_base = float(sem.ReportDefocus())
    log(f"focus_base (ReportDefocus at centroid setup): {focus_base:.4f} um")

    area = load_state(out_dir)
    if area and area.get("holes"):
        log(f"Resuming from {state_path(out_dir)} ({sum(1 for h in area['holes'] if h.get('acquired'))} acquired)")
    else:
        area = darea.build_area(
            sem,
            nav_item_range,
            beam_radius,
            add_overlap=add_overlap,
            log=log,
        )
        area["focus_base"] = focus_base
        save_state(out_dir, area)
        plot_path = os.path.join(out_dir, f"{stem}_lattice.png")
        try:
            darea.plot_area(area, path=plot_path, show=False)
            log(f"Wrote lattice plot {plot_path}")
        except Exception as exc:
            log(f"WARNING: Could not plot lattice: {exc}")

    focus_base = float(area.get("focus_base", focus_base))
    holes = area["holes"]
    neighbors = dsched.neighbor_graph(holes, area["hex_step_um"])
    geometry, extra_ss, extra_z = load_initial_geometry(area)
    if extra_ss:
        live = dpred.geometry_from_holes(
            holes, extra_ss=extra_ss, extra_z=extra_z, requested_resolution=spline_resolution, n_iter=spline_n_iter
        )
        if live is not None:
            geometry = live

    c3_store = [[{"c3_offset": None}] for _ in holes]
    ronchi.configure(c3_position_store=lambda: c3_store)

    n_done = sum(1 for h in holes if h.get("acquired"))
    dtest.mark_center_holes(holes)
    log(f"{len(holes)} holes, {n_done} already acquired, {len(holes) - n_done} remaining")

    if dry_run:
        log("dry_run=True: lattice/state written, not acquiring.")
        sem.SuppressReports(0)
        return

    while True:
        idx, mode = dsched.next_hole(holes, neighbors, area["vertices_ss"])
        if idx is None:
            break
        hole = holes[idx]
        delta_z, used_mode = dpred.predict_delta_z(holes, neighbors, idx, geometry, mode)
        hole["predict_mode"] = used_mode
        log(
            f"\nHole {idx + 1}/{len(holes)} SS ({hole['ssx']:.3f}, {hole['ssy']:.3f}) "
            f"mode={used_mode} delta_z={delta_z:.4f} um"
        )
        if checkDewar:
            try:
                sem.ManageDewarsAndPumps(1)
            except Exception:
                pass
        result = dacq.acquire_hole(
            sem,
            hole,
            stem=stem,
            directory=out_dir,
            target_defocus=science_target,
            ctf_probe_underfocus=ctf_probe_underfocus,
            focus_base=focus_base,
            delta_z=delta_z,
            delay_is=delayIS,
            ctf_defocus_lo=ctfDefocusLo,
            ctf_defocus_hi=ctfDefocusHi,
            count_threshold=count_threshold,
            beam_tilt_comp=beamTiltComp,
            hole_count=n_done + 1,
            correct_astig=correctAstig,
            astig_every_n=astigEveryN,
            test_mode=test_mode,
            test_science_ctf_xtilt=test_science_ctf_xtilt,
            log=log,
        )
        ctf = result.get("ctf") or {}
        measured = ctf.get("defocus_um", float("nan"))
        if np.isfinite(measured):
            hole["z_focus"] = dpred.z_focus_from_probe(
                delta_z, measured, science_target, ctf_probe_underfocus
            )
            hole["defocus_um"] = _finite_or_none(measured)
            hole["astig_x_um"] = _finite_or_none(ctf.get("astig_x_um"))
            hole["astig_y_um"] = _finite_or_none(ctf.get("astig_y_um"))
            live = dpred.geometry_from_holes(
                holes,
                extra_ss=extra_ss,
                extra_z=extra_z,
                requested_resolution=spline_resolution,
                n_iter=spline_n_iter,
            )
            if live is not None:
                geometry = live
                log(f"Spline updated ({live['n_points']} points, res {live['spline']['resolution']})")
        else:
            log("WARNING: Probe CTF invalid; not updating spline for this hole.")
        hole["acquired"] = True
        if test_mode:
            dtest.attach_hole_test(hole, result, science_target, probe_target)
        n_done += 1
        save_state(out_dir, area, extra={"focus_base": focus_base, "n_acquired": n_done})
        if test_mode:
            dtest.write_outputs(holes, out_dir, stem, log=log, echo_table=False)

    log(f"##### decolace_LPP finished: {n_done}/{len(holes)} holes #####")
    try:
        darea.plot_area(area, path=os.path.join(out_dir, f"{stem}_lattice.png"), show=False)
    except Exception:
        pass
    if test_mode:
        dtest.write_outputs(holes, out_dir, stem, log=log)
    sem.SuppressReports(0)


if __name__ == "__main__" or True:
    main()
