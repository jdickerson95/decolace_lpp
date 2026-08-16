#!Python
# ===================================================================
# ScriptName     auto_on_peak
# Purpose:       Measure working-xtilt defocus at CTF X-tilt, then
#                5x5 X-tilt grid within 1 fringe (phase-only CtfFind),
#                fit the maximum-phase X-tilt, and take a +20
#                ronchigram for ronchiTargetPhaseA / B.
# ===================================================================
import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
import csv
import json
import os
from datetime import datetime

try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _here = os.getcwd()
for _p in (_here, os.path.dirname(_here), os.getcwd()):
    if _p and _p not in sys.path:
        sys.path.append(_p)

import matplotlib.pyplot as plt
import numpy as np
import serialem as sem

import PACEtomo_ctf_calibrations as ctfcal
import PACEtomo_ronchi as ronchi

############ SETTINGS ############

# CTF X-tilt (off the laser) for the initial defocus measurement.
# Match decolace_collect.py / PACEtomo.py.
ctfXtiltX = 0.002836
ctfXtiltY = 0.003867
xtilt_lens_index = 2
xtilt_calibration_file = r""  # check_xtilt_defoc_astig JSON; empty = no back-project

# X-tilt change equivalent to 1 fringe. Required. See laser_helper/README.md.
fringe_xtilt_x = 0.0
fringe_xtilt_y = 0.0

# Odd grid size centered on the starting (working) X-tilt, spanning 1 fringe.
grid_size = 5
settle_delay_s = 1.0

# Initial CtfFind search range [um] at CTF X-tilt (underfocus negative).
ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2

# Phase-only CtfFind on the grid. SerialEM requires (hi - lo) < 90 deg.
phase_search_lo = 0.0
phase_search_hi = 80.0

# After the peak fit, Trial ronchigram at C3 = current + working_offset.
working_offset = 20.0
ronchiDelay = 1.0
ronchiBinning = 32
ronchiPixelSize = 0.98e-4 * 2
ronchiPeakRadius = 100

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""

csv_measurements = "auto_on_peak_measurements.csv"
result_json = "auto_on_peak.json"
plot_file = "auto_on_peak_fit.png"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


def prepare_output_paths():
    global csv_measurements, result_json, plot_file
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        csv_measurements = os.path.join(save_dir, os.path.basename(csv_measurements))
        result_json = os.path.join(save_dir, os.path.basename(result_json))
        plot_file = os.path.join(save_dir, os.path.basename(plot_file))
    csv_measurements = os.path.abspath(csv_measurements)
    result_json = os.path.abspath(result_json)
    plot_file = os.path.abspath(plot_file)


def report_xtilt():
    xt = sem.ReportXLensDeflector(xtilt_lens_index)
    return float(xt[0]), float(xt[1])


def set_xtilt(x, y):
    sem.SetXLensDeflector(xtilt_lens_index, float(x), float(y))


def grid_offsets(size, fringe_x, fringe_y):
    size = int(size)
    if size < 3 or size % 2 == 0:
        raise ValueError("grid_size must be an odd integer >= 3")
    half = size // 2
    ticks = [i / half for i in range(-half, half + 1)]
    return [(0.5 * t * fringe_x, 0.5 * s * fringe_y) for s in ticks for t in ticks]


def parse_ctffind(cfind):
    return ctfcal.parse_ctffind(cfind)


def run_ctffind(defocus_lo, defocus_hi, phase_lo=None, phase_hi=None):
    sem.NoMessageBoxOnError(1)
    try:
        if phase_lo is None or phase_hi is None:
            cfind = sem.CtfFind("A", float(defocus_lo), float(defocus_hi))
        else:
            cfind = sem.CtfFind(
                "A",
                float(defocus_lo),
                float(defocus_hi),
                0,
                0,
                float(phase_lo),
                float(phase_hi),
            )
    finally:
        sem.NoMessageBoxOnError(0)
    if len(cfind) == 0:
        echo("WARNING: CtfFind failed.")
        return parse_ctffind([])
    return parse_ctffind(cfind)


def acquire_record():
    sem.GoToLowDoseArea("R")
    sem.R()


def analyze_ronchi(image):
    return ronchi.analyze_ronchigram(
        image,
        ronchiPixelSize,
        ronchiBinning,
        target_phase_a=0.0,
        target_phase_b=0.0,
        correct_ks=[[0.0, 0.0], [0.0, 0.0]],
        peak_radius=ronchiPeakRadius,
    )


def acquire_ronchi(c3):
    sem.GoToLowDoseArea("T")
    sem.SetImageDistanceOffset(c3)
    sem.Delay(ronchiDelay, "s")
    sem.T()
    return analyze_ronchi(np.asarray(sem.bufferImage("A")))


def format_ks(ks):
    ks = np.asarray(ks, dtype=float)
    return (
        f"[[{ks[0, 0]:.6f}, {ks[0, 1]:.6f}], "
        f"[{ks[1, 0]:.6f}, {ks[1, 1]:.6f}]]"
    )


def fit_quadratic_peak(x, y, z):
    """z = a0 + ax x + ay y + axx x^2 + axy x y + ayy y^2. Peak of the paraboloid."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(np.sum(mask))
    empty = {
        "a0": np.nan,
        "ax": np.nan,
        "ay": np.nan,
        "axx": np.nan,
        "axy": np.nan,
        "ayy": np.nan,
        "peak_x": np.nan,
        "peak_y": np.nan,
        "peak_z": np.nan,
        "is_maximum": False,
        "rms": np.nan,
        "n": n,
    }
    if n < 6:
        return empty
    A = np.column_stack((
        np.ones(n),
        x[mask],
        y[mask],
        x[mask] ** 2,
        x[mask] * y[mask],
        y[mask] ** 2,
    ))
    coeffs, _, _, _ = np.linalg.lstsq(A, z[mask], rcond=None)
    a0, ax, ay, axx, axy, ayy = [float(c) for c in coeffs]
    pred = A @ coeffs
    rms = float(np.sqrt(np.mean((pred - z[mask]) ** 2)))
    hess = np.array([[2.0 * axx, axy], [axy, 2.0 * ayy]], dtype=float)
    try:
        peak = np.linalg.solve(hess, np.array([-ax, -ay], dtype=float))
        peak_x, peak_y = float(peak[0]), float(peak[1])
        peak_z = float(
            a0 + ax * peak_x + ay * peak_y
            + axx * peak_x ** 2 + axy * peak_x * peak_y + ayy * peak_y ** 2
        )
    except np.linalg.LinAlgError:
        peak_x, peak_y, peak_z = np.nan, np.nan, np.nan
    det = float(np.linalg.det(hess)) if np.all(np.isfinite(hess)) else np.nan
    is_maximum = bool(np.isfinite(det) and axx < 0 and det > 0)
    return {
        "a0": a0,
        "ax": ax,
        "ay": ay,
        "axx": axx,
        "axy": axy,
        "ayy": ayy,
        "peak_x": peak_x,
        "peak_y": peak_y,
        "peak_z": peak_z,
        "is_maximum": is_maximum,
        "rms": rms,
        "n": n,
    }


def in_grid(x, y, xs, ys):
    return (
        np.nanmin(xs) - 1e-12 <= x <= np.nanmax(xs) + 1e-12
        and np.nanmin(ys) - 1e-12 <= y <= np.nanmax(ys) + 1e-12
    )


def make_plot(path, rows, peak_x, peak_y, discrete_x, discrete_y):
    xs = np.array([r["xtilt_x"] for r in rows], dtype=float)
    ys = np.array([r["xtilt_y"] for r in rows], dtype=float)
    phase = np.array([r["phase_shift_deg"] for r in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(7, 6), tight_layout=True)
    sc = ax.scatter(xs, ys, c=phase, cmap="viridis", s=90)
    fig.colorbar(sc, ax=ax, shrink=0.8, label="CtfFind phase (deg)")
    if np.isfinite(discrete_x) and np.isfinite(discrete_y):
        ax.scatter(
            [discrete_x],
            [discrete_y],
            marker="o",
            facecolors="none",
            edgecolors="white",
            s=160,
            linewidths=1.5,
            label="grid max",
        )
    if np.isfinite(peak_x) and np.isfinite(peak_y):
        ax.scatter(
            [peak_x],
            [peak_y],
            color="C3",
            marker="*",
            s=160,
            zorder=5,
            label="quadratic peak",
        )
    ax.set_xlabel("X-tilt X")
    ax.set_ylabel("X-tilt Y")
    ax.set_title("CTF phase vs X-tilt (1-fringe grid)")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    prepare_output_paths()
    echo("##### Laser helper: auto on peak #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")

    if fringe_xtilt_x == 0.0 or fringe_xtilt_y == 0.0:
        echo(
            "ERROR: set fringe_xtilt_x and fringe_xtilt_y to the X-tilt "
            "change equivalent to 1 fringe. See laser_helper/README.md."
        )
        sem.Exit()
    if phase_search_hi - phase_search_lo >= 90.0:
        echo("ERROR: SerialEM requires phase_search_hi - phase_search_lo < 90 deg.")
        sem.Exit()

    ctfcal.configure(
        sem_module=sem,
        use_ctf_xtilt=True,
        ctf_xtilt_x=ctfXtiltX,
        ctf_xtilt_y=ctfXtiltY,
        xtilt_lens_index_value=xtilt_lens_index,
        xtilt_calibration_path=xtilt_calibration_file,
    )

    start_xtilt = report_xtilt()
    start_c3 = float(sem.ReportImageDistanceOffset())
    echo(
        f"Starting X-tilt (working): ({start_xtilt[0]:.6f}, {start_xtilt[1]:.6f})"
    )
    echo(f"Starting C3: {start_c3:.6f}")
    echo(f"CTF X-tilt: ({ctfXtiltX:.6f}, {ctfXtiltY:.6f})")
    echo(
        f"1-fringe X-tilt: X={fringe_xtilt_x:g}, Y={fringe_xtilt_y:g}; "
        f"grid={grid_size}x{grid_size}"
    )

    offsets = grid_offsets(grid_size, fringe_xtilt_x, fringe_xtilt_y)
    rows = []
    working_defocus = np.nan
    measured_defocus = np.nan
    peak_xtilt = (np.nan, np.nan)
    ronchi_result = None

    try:
        echo("------------------------------------------------")
        echo("Measuring defocus at CTF X-tilt")
        set_xtilt(ctfXtiltX, ctfXtiltY)
        sem.Delay(settle_delay_s, "s")
        acquire_record()
        ctf0 = run_ctffind(ctf_defocus_lo, ctf_defocus_hi)
        measured_defocus = ctf0["defocus_um"]
        echo(
            f"  CtfFind at CTF X-tilt: defocus={measured_defocus:.4f} um, "
            f"phase={ctf0['phase_shift_deg']:.2f} deg, "
            f"res={ctf0['resolution_A']:.2f} A"
        )
        if not np.isfinite(measured_defocus):
            echo("ERROR: initial CtfFind defocus is invalid.")
            sem.Exit()

        working_defocus, _, _ = ctfcal.correct_ctf_to_working_xtilt(
            measured_defocus,
            np.nan,
            np.nan,
            start_xtilt,
            (ctfXtiltX, ctfXtiltY),
        )
        if xtilt_calibration_file:
            echo(
                f"  Back-projected defocus at starting X-tilt: "
                f"{working_defocus:.4f} um"
            )
        else:
            echo(
                "  No xtilt_calibration_file; using CTF-X-tilt defocus as the "
                "fixed value (not back-projected)."
            )
            working_defocus = measured_defocus
        echo(f"Saved working defocus (fixed for phase grid): {working_defocus:.4f} um")

        echo("------------------------------------------------")
        echo("Restoring starting X-tilt for the 1-fringe phase grid")
        set_xtilt(start_xtilt[0], start_xtilt[1])
        sem.Delay(settle_delay_s, "s")

        for i, (dxtilt_x, dxtilt_y) in enumerate(offsets, start=1):
            xtilt_x = start_xtilt[0] + dxtilt_x
            xtilt_y = start_xtilt[1] + dxtilt_y
            echo("------------------------------------------------")
            echo(
                f"{i}/{len(offsets)} X-tilt ({xtilt_x:.6f}, {xtilt_y:.6f})  "
                f"d=({dxtilt_x:.6f}, {dxtilt_y:.6f})"
            )
            set_xtilt(xtilt_x, xtilt_y)
            sem.Delay(settle_delay_s, "s")
            acquire_record()
            ctf = run_ctffind(
                working_defocus,
                working_defocus,
                phase_search_lo,
                phase_search_hi,
            )
            row = {
                "dxtilt_x": float(dxtilt_x),
                "dxtilt_y": float(dxtilt_y),
                "xtilt_x": float(xtilt_x),
                "xtilt_y": float(xtilt_y),
                "defocus_um": ctf["defocus_um"],
                "phase_shift_deg": ctf["phase_shift_deg"],
                "fit_score": ctf["fit_score"],
                "resolution_A": ctf["resolution_A"],
            }
            rows.append(row)
            echo(
                f"  phase={ctf['phase_shift_deg']:.2f} deg, "
                f"score={ctf['fit_score']:.4f}, res={ctf['resolution_A']:.2f} A"
            )

        xs = np.array([r["xtilt_x"] for r in rows], dtype=float)
        ys = np.array([r["xtilt_y"] for r in rows], dtype=float)
        phase = np.array([r["phase_shift_deg"] for r in rows], dtype=float)
        fit = fit_quadratic_peak(xs, ys, phase)

        finite = np.isfinite(phase)
        if not np.any(finite):
            echo("ERROR: no valid phase measurements.")
            sem.Exit()
        imax = int(np.nanargmax(phase))
        discrete_x = float(xs[imax])
        discrete_y = float(ys[imax])
        discrete_z = float(phase[imax])

        echo("================================================")
        echo(
            f"QUADRATIC  phase = a0 + ax X + ay Y + axx X^2 + axy X Y + ayy Y^2  "
            f"(RMS {fit['rms']:.3f} deg)"
        )
        echo(
            f"  peak X-tilt=({fit['peak_x']:.6f}, {fit['peak_y']:.6f}), "
            f"phase={fit['peak_z']:.2f} deg, is_maximum={fit['is_maximum']}"
        )
        echo(
            f"  grid max X-tilt=({discrete_x:.6f}, {discrete_y:.6f}), "
            f"phase={discrete_z:.2f} deg"
        )

        use_fit = (
            fit["is_maximum"]
            and np.isfinite(fit["peak_x"])
            and np.isfinite(fit["peak_y"])
            and in_grid(fit["peak_x"], fit["peak_y"], xs, ys)
        )
        if use_fit:
            peak_xtilt = (fit["peak_x"], fit["peak_y"])
            echo("Using quadratic peak (inside the 1-fringe grid).")
        else:
            peak_xtilt = (discrete_x, discrete_y)
            echo(
                "WARNING: quadratic peak is not a maximum inside the grid; "
                "using the discrete grid maximum."
            )

        echo("------------------------------------------------")
        echo(
            f"Ronchigram at peak X-tilt ({peak_xtilt[0]:.6f}, {peak_xtilt[1]:.6f}), "
            f"C3 = start + {working_offset:g}"
        )
        set_xtilt(peak_xtilt[0], peak_xtilt[1])
        sem.Delay(settle_delay_s, "s")
        ronchi_c3 = start_c3 + working_offset
        ronchi_result = acquire_ronchi(ronchi_c3)
        phases = np.asarray(ronchi_result["phases"], dtype=float)
        ks = np.asarray(ronchi_result["ks"], dtype=float)
        echo(f"  ronchiTargetPhaseA (vertical) = {phases[0]:.8f}")
        echo(f"  ronchiTargetPhaseB (horizontal) = {phases[1]:.8f}")
        echo(f"  ks = {format_ks(ks)}")
        echo("================================================")
        echo("COPY INTO decolace_collect.py AND PACEtomo.py:")
        echo(f"  ronchiTargetPhaseA = {phases[0]:.8f}")
        echo(f"  ronchiTargetPhaseB = {phases[1]:.8f}")
        echo(
            f"Set working XLensDeflector({xtilt_lens_index}) to "
            f"({peak_xtilt[0]:.6f}, {peak_xtilt[1]:.6f})"
        )
        echo("================================================")

        with open(csv_measurements, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

        result = {
            "model": "auto_on_peak",
            "created": datetime.now().isoformat(timespec="seconds"),
            "start_xtilt_x": start_xtilt[0],
            "start_xtilt_y": start_xtilt[1],
            "start_c3": start_c3,
            "ctf_xtilt_x": float(ctfXtiltX),
            "ctf_xtilt_y": float(ctfXtiltY),
            "measured_defocus_um_at_ctf_xtilt": float(measured_defocus),
            "working_defocus_um": float(working_defocus),
            "xtilt_calibration_file": xtilt_calibration_file or None,
            "fringe_xtilt_x": float(fringe_xtilt_x),
            "fringe_xtilt_y": float(fringe_xtilt_y),
            "grid_size": int(grid_size),
            "quadratic_fit": fit,
            "discrete_max_xtilt_x": discrete_x,
            "discrete_max_xtilt_y": discrete_y,
            "discrete_max_phase_deg": discrete_z,
            "peak_xtilt_x": float(peak_xtilt[0]),
            "peak_xtilt_y": float(peak_xtilt[1]),
            "used_quadratic_peak": bool(use_fit),
            "ronchi_c3": float(ronchi_c3),
            "ronchiTargetPhaseA": float(phases[0]),
            "ronchiTargetPhaseB": float(phases[1]),
            "ronchiCorrectKs": ks.tolist(),
        }
        with open(result_json, "w") as fh:
            json.dump(result, fh, indent=2)

        make_plot(
            plot_file,
            rows,
            peak_xtilt[0],
            peak_xtilt[1],
            discrete_x,
            discrete_y,
        )
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved result JSON: {result_json}")
        echo(f"Saved plot: {plot_file}")

    finally:
        sem.SetImageDistanceOffset(start_c3)
        set_xtilt(start_xtilt[0], start_xtilt[1])
        echo(f"Restored ImageDistanceOffset to {start_c3:.6f}")
        echo(
            f"Restored XLensDeflector({xtilt_lens_index}) to "
            f"({start_xtilt[0]:.6f}, {start_xtilt[1]:.6f})"
        )

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
