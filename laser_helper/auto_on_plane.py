#!Python
# ===================================================================
# ScriptName     auto_on_plane
# Purpose:       Sweep C3 (ImageDistanceOffset), measure ronchi fringe
#                spacing, fit the in-plane C3, then report ronchiCorrectKs
#                at plane + working_offset (default −20, same as collect).
# ===================================================================
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

import matplotlib.pyplot as plt
import numpy as np
import serialem as sem

import PACEtomo_ronchi as ronchi

############ SETTINGS ############

# C3 sweep relative to the starting ImageDistanceOffset (assumed roughly on plane).
# Measures at first_step, then first_step+step, ... until max_offset, then the same
# negative branch. Does not measure at 0.
first_step = 20.0
step = 5.0
max_offset = 50.0

# After the fit, measure ks at C3_plane + working_offset. That ks is ronchiCorrectKs.
# Must match collect/PACEtomo ronchiC3Offset (Trial is at plane −20).
working_offset = -20.0

# Trial / FFT (must match collect/PACEtomo ronchi settings).
ronchiDelay = 1.0
ronchiBinning = 32
ronchiPixelSize = 0.98e-4 * 2
ronchiPeakRadius = 100

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""

csv_measurements = "auto_on_plane_measurements.csv"
plot_file = "auto_on_plane_fit.png"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


def prepare_output_paths():
    global csv_measurements, plot_file
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        csv_measurements = os.path.join(save_dir, os.path.basename(csv_measurements))
        plot_file = os.path.join(save_dir, os.path.basename(plot_file))
    csv_measurements = os.path.abspath(csv_measurements)
    plot_file = os.path.abspath(plot_file)


def sweep_deltas():
    n = int(round((max_offset - first_step) / step)) + 1
    if n < 1:
        raise ValueError("max_offset must be >= first_step")
    pos = [first_step + i * step for i in range(n)]
    return pos + [-d for d in pos]


def fringe_magnitude(ks):
    ks = np.asarray(ks, dtype=float)
    return 0.5 * (float(np.linalg.norm(ks[0])) + float(np.linalg.norm(ks[1])))


def analyze_ks(image):
    result = ronchi.analyze_ronchigram(
        image,
        ronchiPixelSize,
        ronchiBinning,
        target_phase_a=0.0,
        target_phase_b=0.0,
        correct_ks=[[0.0, 0.0], [0.0, 0.0]],
        peak_radius=ronchiPeakRadius,
    )
    return np.asarray(result["ks"], dtype=float)


def acquire_ks(c3):
    sem.GoToLowDoseArea("T")
    sem.SetImageDistanceOffset(c3)
    sem.Delay(ronchiDelay, "s")
    sem.T()
    image = np.asarray(sem.bufferImage("A"))
    return analyze_ks(image)


def format_ks(ks):
    ks = np.asarray(ks, dtype=float)
    return (
        f"[[{ks[0, 0]:.6f}, {ks[0, 1]:.6f}], "
        f"[{ks[1, 0]:.6f}, {ks[1, 1]:.6f}]]"
    )


def fit_line(c3, signed_mag):
    x = np.asarray(c3, dtype=float)
    y = np.asarray(signed_mag, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.sum(mask) < 2:
        return {"slope": np.nan, "intercept": np.nan, "r2": np.nan, "rms": np.nan}
    slope, intercept = np.polyfit(x[mask], y[mask], 1)
    pred = slope * x[mask] + intercept
    ss_res = float(np.sum((y[mask] - pred) ** 2))
    ss_tot = float(np.sum((y[mask] - np.mean(y[mask])) ** 2))
    r2 = np.nan if ss_tot == 0 else 1.0 - ss_res / ss_tot
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(r2),
        "rms": float(np.sqrt(np.mean((pred - y[mask]) ** 2))),
    }


def make_plot(path, rows, fit, c3_plane, verify_c3, verify_signed):
    c3 = np.array([float(r["c3"]) for r in rows], dtype=float)
    signed = np.array([float(r["signed_magnitude"]) for r in rows], dtype=float)
    x_line = np.linspace(np.min(c3), np.max(c3), 200)
    y_line = fit["slope"] * x_line + fit["intercept"]

    fig, ax = plt.subplots(figsize=(7, 5), tight_layout=True)
    ax.scatter(c3, signed, label="sweep")
    ax.plot(x_line, y_line, color="C1", label="linear fit")
    if np.isfinite(c3_plane):
        ax.axvline(c3_plane, color="gray", linestyle="--", linewidth=1, label="C3 plane")
    ax.axhline(0, color="gray", linestyle=":", linewidth=1)
    if np.isfinite(verify_c3) and np.isfinite(verify_signed):
        ax.scatter(
            [verify_c3],
            [verify_signed],
            color="C3",
            marker="*",
            s=120,
            zorder=5,
            label="plane + working_offset",
        )
    ax.set_xlabel("C3 ImageDistanceOffset")
    ax.set_ylabel("Signed fringe-spacing magnitude (1/um)")
    ax.set_title("Ronchi fringe spacing vs C3")
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    prepare_output_paths()
    echo("##### Laser helper: auto on plane #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")
    echo(
        f"sweep first={first_step:g}, step={step:g}, max={max_offset:g}; "
        f"working_offset={working_offset:g}"
    )

    start_c3 = float(sem.ReportImageDistanceOffset())
    echo(f"Starting C3 (assumed roughly on plane): {start_c3:.6f}")
    deltas = sweep_deltas()
    echo(f"Deltas: {', '.join(f'{d:+g}' for d in deltas)}")

    rows = []
    verify_c3 = np.nan
    verify_ks = None
    verify_signed = np.nan
    c3_plane = np.nan
    try:
        for i, delta in enumerate(deltas):
            c3 = start_c3 + delta
            echo("------------------------------------------------")
            echo(f"C3 = {c3:.6f} (start {delta:+g})")
            ks = acquire_ks(c3)
            mag = fringe_magnitude(ks)
            signed = mag if delta > 0 else -mag
            row = {
                "step": i,
                "delta": float(delta),
                "c3": float(c3),
                "magnitude": mag,
                "signed_magnitude": signed,
                "ks_00": float(ks[0, 0]),
                "ks_01": float(ks[0, 1]),
                "ks_10": float(ks[1, 0]),
                "ks_11": float(ks[1, 1]),
            }
            rows.append(row)
            echo(
                f"  |ks|={mag:.4f}, signed={signed:.4f}, ks={format_ks(ks)}"
            )

        fit = fit_line(
            [r["c3"] for r in rows],
            [r["signed_magnitude"] for r in rows],
        )
        slope = fit["slope"]
        intercept = fit["intercept"]
        if not np.isfinite(slope) or slope == 0:
            echo("ERROR: linear fit failed (slope is zero or invalid).")
            sem.Exit()
        c3_plane = -intercept / slope

        echo("================================================")
        echo(
            f"FIT  signed_mag = {slope:.6f} * C3 + {intercept:.6f}  "
            f"(R^2 {fit['r2']:.4f}, RMS {fit['rms']:.4f})"
        )
        echo(f"C3 on plane (signed |ks| = 0): {c3_plane:.6f}")

        verify_c3 = c3_plane + working_offset
        echo("------------------------------------------------")
        echo(f"Verification C3 = plane + {working_offset:g} = {verify_c3:.6f}")
        verify_ks = acquire_ks(verify_c3)
        verify_mag = fringe_magnitude(verify_ks)
        verify_signed = verify_mag if working_offset > 0 else -verify_mag
        echo(f"  |ks|={verify_mag:.4f}, ks={format_ks(verify_ks)}")
        echo("================================================")
        echo("COPY INTO decolace_collect.py AND PACEtomo.py:")
        echo(f"  ronchiCorrectKs = {format_ks(verify_ks)}")
        echo(f"  ronchiC3Offset = {working_offset:g}")
        echo(f"Set session ImageDistanceOffset (C3 baseline) to {c3_plane:.6f}")
        echo("================================================")

        with open(csv_measurements, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

        make_plot(plot_file, rows, fit, c3_plane, verify_c3, verify_signed)
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved plot: {plot_file}")

    finally:
        sem.SetImageDistanceOffset(start_c3)
        echo(f"Restored ImageDistanceOffset to {start_c3:.6f}")

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
