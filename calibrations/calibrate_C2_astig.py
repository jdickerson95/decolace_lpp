#!Python
# ===================================================================
# ScriptName     calibrate_C2_astig
# Purpose:       Sweep C2 (condenser) stigmator on a grid, collect a
#                ronchigram at each point, measure fringe spacing in x
#                and y, and fit
#                  [ks_x, ks_y] = M @ [dC2_x, dC2_y] + b
#                so auto_c2_stig.py can equalize the two spacings with
#                  dC2 = inv(M) @ ([mean, mean] - [ks_x, ks_y])
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

import PACEtomo_ronchi as ronchi

############ SETTINGS ############

# Odd grid size: 3, 5, or 7.
grid_size = 5

# Step between neighboring C2 stigmator points (SerialEM units, typically -1..1).
stig_step = 0.02

# Settle time after SetCondenserStigmator [s].
settle_delay_s = 1.0

# Trial ronchigram at C3 = current + working_offset.
# Must match collect/PACEtomo ronchiC3Offset (Trial is at plane −20).
working_offset = -20.0
ronchiDelay = 1.0
ronchiBinning = 32
ronchiPixelSize = 0.98e-4 * 2
ronchiPeakRadius = 100

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""

csv_measurements = "c2_astig_measurements.csv"
calibration_json = "c2_astig_calibration.json"
plot_file = "c2_astig_fit.png"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


def prepare_output_paths():
    global csv_measurements, calibration_json, plot_file
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        csv_measurements = os.path.join(save_dir, os.path.basename(csv_measurements))
        calibration_json = os.path.join(save_dir, os.path.basename(calibration_json))
        plot_file = os.path.join(save_dir, os.path.basename(plot_file))
    csv_measurements = os.path.abspath(csv_measurements)
    calibration_json = os.path.abspath(calibration_json)
    plot_file = os.path.abspath(plot_file)


def grid_offsets(size, step):
    size = int(size)
    if size not in (3, 5, 7):
        raise ValueError("grid_size must be 3, 5, or 7")
    half = size // 2
    ticks = [i * step for i in range(-half, half + 1)]
    return [(dx, dy) for dy in ticks for dx in ticks]


def clip_stig(value):
    return float(np.clip(value, -1.0, 1.0))


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
    return analyze_ks(np.asarray(sem.bufferImage("A")))


def fringe_xy(ks):
    """Fringe spatial frequencies (1/µm) of the two ronchi axes.

    ks[0] is the vertical-fringe (mostly x) peak; ks[1] is the
    horizontal-fringe (mostly y) peak. Magnitudes match auto_on_plane.
    """
    ks = np.asarray(ks, dtype=float)
    return float(np.linalg.norm(ks[0])), float(np.linalg.norm(ks[1]))


def format_ks(ks):
    ks = np.asarray(ks, dtype=float)
    return (
        f"[[{ks[0, 0]:.6f}, {ks[0, 1]:.6f}], "
        f"[{ks[1, 0]:.6f}, {ks[1, 1]:.6f}]]"
    )


def fit_c2_map(dstig_x, dstig_y, ks_x, ks_y):
    """Least-squares [ks_x, ks_y] = M @ [dC2_x, dC2_y] + b."""
    dsx = np.asarray(dstig_x, dtype=float).ravel()
    dsy = np.asarray(dstig_y, dtype=float).ravel()
    ax = np.asarray(ks_x, dtype=float).ravel()
    ay = np.asarray(ks_y, dtype=float).ravel()
    mask = np.isfinite(dsx) & np.isfinite(dsy) & np.isfinite(ax) & np.isfinite(ay)
    n = int(np.sum(mask))
    if n < 3:
        raise ValueError(f"Need at least 3 finite C2/ks points (got {n}).")
    A = np.column_stack((dsx[mask], dsy[mask], np.ones(n)))
    cx, _, _, _ = np.linalg.lstsq(A, ax[mask], rcond=None)
    cy, _, _, _ = np.linalg.lstsq(A, ay[mask], rcond=None)
    M = np.array([[cx[0], cx[1]], [cy[0], cy[1]]], dtype=float)
    b = np.array([cx[2], cy[2]], dtype=float)
    det = float(np.linalg.det(M))
    if abs(det) < 1e-12:
        raise ValueError("C2-to-fringe-spacing matrix is singular; increase stig_step.")
    inv_M = np.linalg.inv(M)
    pred = (M @ np.vstack((dsx[mask], dsy[mask]))).T + b
    meas = np.column_stack((ax[mask], ay[mask]))
    rms = float(np.sqrt(np.mean(np.sum((pred - meas) ** 2, axis=1))))
    return {
        "M": M.tolist(),
        "inv_M": inv_M.tolist(),
        "b": b.tolist(),
        "rms_ks": rms,
        "n": n,
        "det": det,
        "formula": (
            "[ks_x, ks_y] = M @ [dC2_x, dC2_y] + b; "
            "dC2 = inv(M) @ ([mean, mean] - [ks_x, ks_y]) to equalize spacings"
        ),
    }


def make_plot(path, rows, fit):
    dsx = np.array([r["dc2_x"] for r in rows], dtype=float)
    dsy = np.array([r["dc2_y"] for r in rows], dtype=float)
    kx = np.array([r["ks_x"] for r in rows], dtype=float)
    ky = np.array([r["ks_y"] for r in rows], dtype=float)
    M = np.array(fit["M"], dtype=float)
    b = np.array(fit["b"], dtype=float)
    pred = (M @ np.vstack((dsx, dsy))).T + b
    pred_x, pred_y = pred[:, 0], pred[:, 1]
    delta = kx - ky

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)
    axes[0].scatter(pred_x, kx, label="ks_x")
    axes[0].scatter(pred_y, ky, marker="x", label="ks_y")
    lo = np.nanmin([pred_x, pred_y, kx, ky])
    hi = np.nanmax([pred_x, pred_y, kx, ky])
    axes[0].plot([lo, hi], [lo, hi], "--", color="gray")
    axes[0].set_xlabel("Fitted fringe spacing (1/um)")
    axes[0].set_ylabel("Measured fringe spacing (1/um)")
    axes[0].set_title("C2 stigmator map fit")
    axes[0].legend(fontsize=8)

    sc = axes[1].scatter(dsx, dsy, c=delta, cmap="coolwarm", s=80)
    fig.colorbar(sc, ax=axes[1], shrink=0.8, label="ks_x - ks_y (1/um)")
    axes[1].axhline(0, color="gray", linestyle="--", linewidth=1)
    axes[1].axvline(0, color="gray", linestyle="--", linewidth=1)
    axes[1].set_xlabel("dC2_x (SerialEM units)")
    axes[1].set_ylabel("dC2_y (SerialEM units)")
    axes[1].set_title("Fringe-spacing difference vs C2")
    axes[1].set_aspect("equal", adjustable="box")

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    prepare_output_paths()
    offsets = grid_offsets(grid_size, stig_step)
    echo("##### C2 (condenser) stigmator fringe-spacing calibration #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")
    echo(
        f"grid={grid_size}x{grid_size}, step={stig_step}, n={len(offsets)}, "
        f"working_offset={working_offset:g}"
    )

    orig_c2 = sem.ReportCondenserStigmator()
    start_sx, start_sy = float(orig_c2[0]), float(orig_c2[1])
    start_c3 = float(sem.ReportImageDistanceOffset())
    ronchi_c3 = start_c3 + working_offset
    echo(f"Saved CondenserStigmator: ({start_sx:.6f}, {start_sy:.6f})")
    echo(f"Starting C3: {start_c3:.6f}; Trial C3 = start + {working_offset:g} = {ronchi_c3:.6f}")

    rows = []
    try:
        for i, (dsx, dsy) in enumerate(offsets, start=1):
            sx = clip_stig(start_sx + dsx)
            sy = clip_stig(start_sy + dsy)
            echo("------------------------------------------------")
            echo(
                f"{i}/{len(offsets)} C2 ({sx:.4f}, {sy:.4f})  "
                f"d=({dsx:.4f}, {dsy:.4f})"
            )
            sem.SetCondenserStigmator(sx, sy)
            sem.Delay(settle_delay_s, "s")
            ks = acquire_ks(ronchi_c3)
            ks_x, ks_y = fringe_xy(ks)
            row = {
                "dc2_x": float(dsx),
                "dc2_y": float(dsy),
                "c2_x": float(sx),
                "c2_y": float(sy),
                "ks_x": ks_x,
                "ks_y": ks_y,
                "ks_diff": ks_x - ks_y,
                "ks_00": float(ks[0, 0]),
                "ks_01": float(ks[0, 1]),
                "ks_10": float(ks[1, 0]),
                "ks_11": float(ks[1, 1]),
            }
            rows.append(row)
            echo(
                f"  ks_x={ks_x:.4f}, ks_y={ks_y:.4f}, "
                f"diff={ks_x - ks_y:.4f}, ks={format_ks(ks)}"
            )
        fit = fit_c2_map(
            [r["dc2_x"] for r in rows],
            [r["dc2_y"] for r in rows],
            [r["ks_x"] for r in rows],
            [r["ks_y"] for r in rows],
        )
        calibration = {
            "model": "c2_stigmator_fringe_spacing",
            "created": datetime.now().isoformat(timespec="seconds"),
            "grid_size": int(grid_size),
            "stig_step": float(stig_step),
            "settle_delay_s": float(settle_delay_s),
            "working_offset": float(working_offset),
            "ronchiBinning": int(ronchiBinning),
            "ronchiPixelSize": float(ronchiPixelSize),
            "ronchiPeakRadius": float(ronchiPeakRadius),
            "start_c2_x": start_sx,
            "start_c2_y": start_sy,
            "start_c3": start_c3,
            **fit,
        }
        with open(csv_measurements, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        with open(calibration_json, "w") as fh:
            json.dump(calibration, fh, indent=2)
        make_plot(plot_file, rows, fit)
        echo("================================================")
        echo("CALIBRATION")
        echo(f"  M = {fit['M']}")
        echo(f"  inv(M) = {fit['inv_M']}")
        echo(f"  b = {fit['b']}")
        echo(f"  RMS {fit['rms_ks']:.6f} 1/um  (n={fit['n']})")
        echo("  To equalize ks_x and ks_y: dC2 = inv(M) @ ([mean, mean] - [ks_x, ks_y])")
        echo("================================================")
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved calibration JSON: {calibration_json}")
        echo(f"Saved plot: {plot_file}")
        echo("Paste the JSON path into laser_helper/auto_c2_stig.py → c2_astig_calibration_file")
    finally:
        sem.SetImageDistanceOffset(start_c3)
        sem.SetCondenserStigmator(start_sx, start_sy)
        echo(f"Restored ImageDistanceOffset to {start_c3:.6f}")
        echo(f"Restored CondenserStigmator to ({start_sx:.6f}, {start_sy:.6f})")

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
