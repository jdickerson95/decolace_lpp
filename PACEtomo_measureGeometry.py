#!Python
# ===================================================================
# ScriptName     PACEtomo_measureGeometry
# Purpose:       Measure lamella geometry from navigator points: plane
#                (pretilt/rotation) plus a cubic B-spline residual so
#                PACEtomo can use z(SSX, SSY) that is not purely planar.
# ===================================================================

############ SETTINGS ############

# Navigator item indices for geo points (at least 3). First item is the
# specimen origin (SSX=SSY=0) used by PACEtomo target coordinates.
geo_nav_item_list = [1, 2, 3]

# beam_tilt | beam_tilt_sem | ctf | beam_tilt_then_ctf
# beam_tilt_then_ctf: large-range beam-tilt, then CtfFind after bringing
# defocus into the CTF search range.
geo_defocus_method = "beam_tilt_then_ctf"

settle_delay_s = 5.0

# Spline residual control-point grid (SSX, SSY). 3x3 is a gentle bowl/saddle.
spline_resolution = (3, 3)
spline_n_iter = 800

# CTF search / refine
useCtfXtilt = True              # set ctfXtilt for CtfFind (off laser); False = measure at working X-tilt
ctfXtiltX = 0.002836
ctfXtiltY = 0.003867
xtilt_calibration_file = ""     # JSON from check_xtilt_defoc_astig.py; empty = no back-project
defocus_error_file = ""         # JSON from calibrate_defocus_error.py; scales ChangeFocus
ctfDefocusLo = -12.0
ctfDefocusHi = -0.2
ctf_target_defocus = -4.0
ctf_resolution_max_A = 20.0
ctf_max_attempts = 3
ctf_retry_delay_s = 5

# Beam-tilt (match PACEtomo / calibrate_beam_tilt_scaling.py)
tilt_angle_mrad = 10.0
beam_tilt_correction = 1.73
defocus_tilt_correction = beam_tilt_correction
beam_tilt_xtilt_x = 0.0
beam_tilt_xtilt_y = 0.0
hasXLens = True
spherical_aberration_mm = 2.7
measure_cycles = 1

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""
csv_measurements = "geometry_measurements.csv"
geometry_json = "geometry.json"
plot_file = "geometry_fit.png"

########## END SETTINGS ##########

import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
import csv
import os
from datetime import datetime

try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _here = os.getcwd()
for _p in (_here, os.path.dirname(_here), os.getcwd()):
    if _p and _p not in sys.path:
        sys.path.append(_p)

import numpy as np
import serialem as sem

import PACEtomo_beamTiltDefocus as btdef
import PACEtomo_geometry as pacegeo
import PACEtomo_ctf_calibrations as ctfcal

btdef.configure(sem_module=sem, has_x_lens=hasXLens)


def echo(text):
    sem.Echo(text)


ctfcal.configure(
    sem_module=sem,
    logger=echo,
    has_x_lens=hasXLens,
    use_ctf_xtilt=useCtfXtilt,
    ctf_xtilt_x=ctfXtiltX,
    ctf_xtilt_y=ctfXtiltY,
    xtilt_calibration_path=xtilt_calibration_file,
    defocus_error_path=defocus_error_file,
)


def prepare_output_paths():
    global csv_measurements, geometry_json, plot_file
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        csv_measurements = os.path.join(save_dir, os.path.basename(csv_measurements))
        geometry_json = os.path.join(save_dir, os.path.basename(geometry_json))
        plot_file = os.path.join(save_dir, os.path.basename(plot_file))
    csv_measurements = os.path.abspath(csv_measurements)
    geometry_json = os.path.abspath(geometry_json)
    plot_file = os.path.abspath(plot_file)


def beam_tilt_measure_defocus():
    xt_x = beam_tilt_xtilt_x if hasXLens else None
    xt_y = beam_tilt_xtilt_y if hasXLens else None
    return btdef.measure_defocus(
        tilt_angle_mrad=tilt_angle_mrad,
        beam_tilt_correction=beam_tilt_correction,
        defocus_tilt_correction=defocus_tilt_correction,
        xtilt_x=xt_x,
        xtilt_y=xt_y,
        cs_mm=spherical_aberration_mm,
    )


def serialEM_measure_defocus():
    xt_x = beam_tilt_xtilt_x if hasXLens else None
    xt_y = beam_tilt_xtilt_y if hasXLens else None
    return btdef.measure_serialEM_defocus(
        xtilt_x=xt_x,
        xtilt_y=xt_y,
    )


def ctf_measure_defocus():
    """CtfFind defocus [um] and resolution [A], back-projected to working X-tilt."""
    ctf = ctfcal.acquire_ctf(
        ctfDefocusLo,
        ctfDefocusHi,
        shot="F",
        max_attempts=ctf_max_attempts,
        resolution_max_A=ctf_resolution_max_A,
        retry_delay_s=ctf_retry_delay_s,
    )
    return float(ctf["defocus_um"]), float(ctf["resolution_A"])


def _measure_beam_tilt(kind):
    defocus = np.nan
    speed_x = speed_y = 0.0
    for _ in range(measure_cycles):
        if kind == "beam_tilt_sem":
            defocus, speed_x, speed_y = serialEM_measure_defocus()
        else:
            defocus, speed_x, speed_y = beam_tilt_measure_defocus()
    return float(defocus), np.array([speed_x, speed_y])


def _in_ctf_range(defocus):
    return np.isfinite(defocus) and ctfDefocusLo <= defocus <= ctfDefocusHi


def measure_point_defocus():
    """Return (defocus_um, drift_xy, method_used, ctf_resolution_A)."""
    method = geo_defocus_method
    if method not in ("beam_tilt", "beam_tilt_sem", "ctf", "beam_tilt_then_ctf"):
        sem.OKBox(
            f"ERROR: Unknown geo_defocus_method '{method}'. "
            "Use 'beam_tilt', 'beam_tilt_sem', 'ctf', or 'beam_tilt_then_ctf'."
        )
        sem.Exit()
    if method == "ctf":
        defocus, resolution = ctf_measure_defocus()
        return float(defocus), np.array([0.0, 0.0]), "ctf", float(resolution)
    if method in ("beam_tilt", "beam_tilt_sem"):
        defocus, drift = _measure_beam_tilt(method)
        echo(f"Defocus: {defocus:.4f} um, drift=({drift[0]:.3f}, {drift[1]:.3f}) nm/s")
        return defocus, drift, method, np.nan

    bt_kind = "beam_tilt"
    d_bt, drift = _measure_beam_tilt(bt_kind)
    echo(f"Beam-tilt defocus: {d_bt:.4f} um, drift=({drift[0]:.3f}, {drift[1]:.3f}) nm/s")
    if not np.isfinite(d_bt):
        return d_bt, drift, "beam_tilt", np.nan

    delta = 0.0
    if not _in_ctf_range(d_bt):
        target = float(ctf_target_defocus)
        if not _in_ctf_range(target):
            target = 0.5 * (ctfDefocusLo + ctfDefocusHi)
        delta = target - d_bt
        echo(f"ChangeFocus desired {delta:.4f} um to bring CTF into range (target {target:.3f} um)")
        commanded = ctfcal.change_focus_for_desired_delta(delta)
        echo(f"  commanded ChangeFocus({commanded:.4f})")
    try:
        d_ctf, resolution = ctf_measure_defocus()
    finally:
        if delta != 0.0:
            undo = ctfcal.change_focus_for_desired_delta(-delta)
            echo(f"Restored focus with ChangeFocus({undo:.4f})")
    if np.isfinite(d_ctf) and np.isfinite(resolution) and resolution <= ctf_resolution_max_A:
        d_orig = d_ctf - delta
        echo(f"Using CTF defocus in original focus state: {d_orig:.4f} um")
        return float(d_orig), drift, "beam_tilt_then_ctf", float(resolution)
    echo("WARNING: CTF refine failed or resolution too poor; using beam-tilt defocus.")
    return float(d_bt), drift, "beam_tilt", float(resolution) if np.isfinite(resolution) else np.nan


def nav_stage_xy(nav_idx):
    item = sem.ReportOtherItem(int(nav_idx))
    return float(item[1]), float(item[2])


def main():
    sem.SuppressReports()
    prepare_output_paths()
    nav_ids = [int(x) for x in geo_nav_item_list]
    if len(nav_ids) < 3:
        sem.OKBox("You need at least 3 navigator points in geo_nav_item_list.")
        sem.Exit()

    echo("##### PACEtomo measure geometry #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")
    echo(f"Nav items: {nav_ids}")
    echo(f"Defocus method: {geo_defocus_method}")
    echo(f"Spline resolution: {spline_resolution}")

    origin_idx = nav_ids[0]
    origin_x, origin_y = nav_stage_xy(origin_idx)
    echo(f"Origin nav item {origin_idx}: stage ({origin_x:.3f}, {origin_y:.3f})")

    sem.GoToLowDoseArea("R")
    s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
    sem.SetImageShift(0, 0)

    rows = []
    try:
        for i, nav_idx in enumerate(nav_ids, start=1):
            echo("------------------------------------------------")
            echo(f"Geo point {i}/{len(nav_ids)}: nav item {nav_idx}")
            sem.SetImageShift(0, 0)
            sem.MoveToNavItem(int(nav_idx))
            stage_x, stage_y = nav_stage_xy(nav_idx)
            ssx, ssy = s2ssMatrix @ np.array([stage_x - origin_x, stage_y - origin_y])
            echo(f"  SSX={ssx:.3f} um, SSY={ssy:.3f} um")
            sem.GoToLowDoseArea("R")
            sem.SetImageShift(0, 0)
            sem.Delay(settle_delay_s, "s")
            defocus, drift, method_used, resolution = measure_point_defocus()
            row = {
                "nav_item": int(nav_idx),
                "ssx": float(ssx),
                "ssy": float(ssy),
                "stage_x": float(stage_x),
                "stage_y": float(stage_y),
                "defocus_um": float(defocus) if np.isfinite(defocus) else np.nan,
                "method_used": method_used,
                "ctf_resolution_A": float(resolution) if np.isfinite(resolution) else np.nan,
                "drift_x_nm_per_s": float(drift[0]),
                "drift_y_nm_per_s": float(drift[1]),
            }
            rows.append(row)
            echo(
                f"  defocus={row['defocus_um']:.4f} um ({method_used}), "
                f"res={row['ctf_resolution_A']:.2f} A"
            )
            sem.SetImageShift(0, 0)
    finally:
        sem.SetImageShift(0, 0)
        sem.MoveToNavItem(int(origin_idx))
        echo(f"Returned to origin nav item {origin_idx}")

    kept = []
    for row in rows:
        if abs(float(row["defocus_um"])) >= 0.01 and np.isfinite(row["defocus_um"]):
            kept.append(row)
        else:
            echo(
                f"WARNING: Discarding nav item {row['nav_item']}: "
                "measured defocus is 0 or invalid."
            )
    if len(kept) < 3:
        sem.OKBox(
            f"Not enough valid geo points ({len(kept)}). Need at least 3."
        )
        sem.Exit()

    origin_defocus = np.nan
    for row in rows:
        if row["nav_item"] == origin_idx:
            origin_defocus = float(row["defocus_um"])
            break
    if not np.isfinite(origin_defocus) or abs(origin_defocus) < 0.01:
        sem.OKBox(
            "Origin nav item defocus is invalid. The first geo_nav_item_list "
            "entry must yield a usable measurement (z=0 at that point)."
        )
        sem.Exit()
    echo(f"Origin defocus (subtracted for z): {origin_defocus:.4f} um")

    ssx = np.array([r["ssx"] for r in kept], dtype=float)
    ssy = np.array([r["ssy"] for r in kept], dtype=float)
    z = np.array([r["defocus_um"] - origin_defocus for r in kept], dtype=float)
    for row in rows:
        if np.isfinite(row["defocus_um"]):
            row["z_um"] = float(row["defocus_um"] - origin_defocus)
        else:
            row["z_um"] = np.nan

    try:
        geometry = pacegeo.build_geometry(
            ssx, ssy, z, spline_resolution=spline_resolution, spline_n_iter=spline_n_iter
        )
    except ImportError as exc:
        echo(f"ERROR: {exc}")
        sem.OKBox(str(exc))
        sem.Exit()

    geometry["geo_defocus_method"] = geo_defocus_method
    geometry["origin_nav_item"] = int(origin_idx)
    geometry["origin_defocus_um"] = float(origin_defocus)
    geometry["geo_nav_item_list"] = nav_ids

    with open(csv_measurements, "w", newline="") as fh:
        fieldnames = list(rows[0].keys())
        if "z_um" not in fieldnames:
            fieldnames.append("z_um")
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    pacegeo.save_geometry(geometry_json, geometry)
    pacegeo.make_plot(plot_file, geometry)

    echo("================================================")
    echo(f"Fitted plane + spline into {geometry['n_points']} points "
         f"({len(rows) - len(kept)} discarded).")
    echo(f"  pretilt:  {geometry['pretilt']:.2f} deg")
    echo(f"  rotation: {geometry['rotation']:.2f} deg")
    echo(f"  plane RMS:  {geometry['plane_rms_um']:.4f} um")
    echo(f"  spline RMS: {geometry['spline_rms_um']:.4f} um (residual)")
    echo(f"  total RMS:  {geometry['fit_rms_um']:.4f} um")
    if geometry["plane_mean_error"] > 1:
        echo("WARNING: Plane fit shows large error. Check autofocus/CTF and geo points.")
    if abs(geometry["pretilt"]) >= 30 and abs(180 - abs(geometry["pretilt"])) >= 30:
        echo("WARNING: Pretilt value seems abnormally high.")
    echo(f"Set PACEtomo geometryMode='spline' and geometry_file to:")
    echo(f"  {geometry_json}")
    echo("================================================")
    echo(f"Saved measurements: {csv_measurements}")
    echo(f"Saved geometry JSON: {geometry_json}")
    echo(f"Saved plot: {plot_file}")

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
