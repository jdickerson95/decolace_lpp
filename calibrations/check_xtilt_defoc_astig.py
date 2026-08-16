#!Python
# ===================================================================
# ScriptName     check_xtilt_defoc_astig
# Purpose:       Measure how X-tilt changes CtfFind defocus and
#                astigmatism on a 3x3, 5x5, or 7x7 grid around a
#                starting X-tilt. Fit planes
#                  z = a0 + ax * dX + ay * dY
#                so a given X-tilt change (dX, dY) predicts
#                  delta_defocus_um = ax * dX + ay * dY
#                  delta_astig_um   = bx * dX + by * dY
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

############ SETTINGS ############

# Center of the X-tilt grid (XLensDeflector values).
start_xtilt_x = 0.0
start_xtilt_y = 0.0

# Odd grid size: 3, 5, or 7.
grid_size = 5

# Step between neighboring grid points (same units as XLensDeflector).
xtilt_step = 0.0005

# Settle time after setting X-tilt [s].
settle_delay_s = 1.0

xtilt_lens_index = 2

# CTF search range [microns].
ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""

csv_measurements = "xtilt_defoc_astig_measurements.csv"
calibration_json = "xtilt_defoc_astig_calibration.json"
plot_file = "xtilt_defoc_astig_fits.png"

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


def parse_ctffind(cfind):
    return ctfcal.parse_ctffind(cfind)


def run_ctffind():
    sem.NoMessageBoxOnError(1)
    try:
        cfind = sem.CtfFind("A", ctf_defocus_lo, ctf_defocus_hi)
    finally:
        sem.NoMessageBoxOnError(0)
    if len(cfind) == 0:
        echo("WARNING: CtfFind failed.")
        return parse_ctffind([])
    return parse_ctffind(cfind)


def astig_components(astig_um, angle_deg):
    return ctfcal.astig_components(astig_um, angle_deg)


def fit_plane(dx, dy, z):
    """Least-squares z = a0 + ax * dx + ay * dy."""
    x = np.array(dx, dtype=float)
    y = np.array(dy, dtype=float)
    zz = np.array(z, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(zz)
    n = int(np.sum(mask))
    if n < 3:
        return {
            "a0": np.nan,
            "ax": np.nan,
            "ay": np.nan,
            "rms": np.nan,
            "n": n,
            "grad": np.nan,
        }
    A = np.column_stack((np.ones(n), x[mask], y[mask]))
    coeffs, _, _, _ = np.linalg.lstsq(A, zz[mask], rcond=None)
    a0, ax, ay = [float(c) for c in coeffs]
    pred = a0 + ax * x[mask] + ay * y[mask]
    rms = float(np.sqrt(np.mean((pred - zz[mask]) ** 2)))
    return {
        "a0": a0,
        "ax": ax,
        "ay": ay,
        "rms": rms,
        "n": n,
        "grad": float(np.hypot(ax, ay)),
    }


def predict_delta(fit, dxtilt_x, dxtilt_y):
    ax, ay = fit["ax"], fit["ay"]
    if not np.isfinite(ax) or not np.isfinite(ay):
        return np.nan
    return float(ax * dxtilt_x + ay * dxtilt_y)


def make_plot(path, rows, defocus_fit, astig_fit):
    dx = np.array([r["dxtilt_x"] for r in rows], dtype=float)
    dy = np.array([r["dxtilt_y"] for r in rows], dtype=float)
    defocus = np.array([r["defocus_um"] for r in rows], dtype=float)
    astig = np.array([r["astig_um"] for r in rows], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), tight_layout=True)

    for ax, values, title, cmap in (
        (axes[0, 0], defocus, "Defocus (um)", "viridis"),
        (axes[0, 1], astig, "Astigmatism (um)", "magma"),
    ):
        sc = ax.scatter(dx, dy, c=values, cmap=cmap, s=80)
        fig.colorbar(sc, ax=ax, shrink=0.8)
        ax.set_xlabel("dX-tilt")
        ax.set_ylabel("dY-tilt")
        ax.set_title(title)
        ax.set_aspect("equal", adjustable="box")

    defocus_pred = defocus_fit["a0"] + defocus_fit["ax"] * dx + defocus_fit["ay"] * dy
    astig_pred = astig_fit["a0"] + astig_fit["ax"] * dx + astig_fit["ay"] * dy
    axes[1, 0].scatter(defocus_pred, defocus)
    lo = np.nanmin([defocus_pred, defocus])
    hi = np.nanmax([defocus_pred, defocus])
    axes[1, 0].plot([lo, hi], [lo, hi], "--", color="gray")
    axes[1, 0].set_xlabel("Plane-fit defocus (um)")
    axes[1, 0].set_ylabel("Measured defocus (um)")
    axes[1, 0].set_title("Defocus plane fit")

    axes[1, 1].scatter(astig_pred, astig)
    lo = np.nanmin([astig_pred, astig])
    hi = np.nanmax([astig_pred, astig])
    axes[1, 1].plot([lo, hi], [lo, hi], "--", color="gray")
    axes[1, 1].set_xlabel("Plane-fit astigmatism (um)")
    axes[1, 1].set_ylabel("Measured astigmatism (um)")
    axes[1, 1].set_title("Astigmatism plane fit")

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    prepare_output_paths()
    offsets = grid_offsets(grid_size, xtilt_step)
    echo("##### X-tilt defocus / astigmatism check #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")
    echo(
        f"start X-tilt=({start_xtilt_x:.6f}, {start_xtilt_y:.6f}), "
        f"grid={grid_size}x{grid_size}, step={xtilt_step}, "
        f"n={len(offsets)}, settle={settle_delay_s:.3f} s"
    )

    original_xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
    echo(
        f"Saved XLensDeflector({xtilt_lens_index}): "
        f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
    )

    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)

    rows = []
    try:
        for i, (dxtilt_x, dxtilt_y) in enumerate(offsets, start=1):
            xtilt_x = start_xtilt_x + dxtilt_x
            xtilt_y = start_xtilt_y + dxtilt_y
            echo("------------------------------------------------")
            echo(
                f"{i}/{len(offsets)} X-tilt ({xtilt_x:.6f}, {xtilt_y:.6f})  "
                f"d=({dxtilt_x:.6f}, {dxtilt_y:.6f})"
            )
            sem.SetXLensDeflector(xtilt_lens_index, float(xtilt_x), float(xtilt_y))
            sem.Delay(settle_delay_s, "s")
            sem.L()
            ctf = run_ctffind()
            astig_x, astig_y = astig_components(ctf["astig_um"], ctf["astig_angle_deg"])
            row = {
                "dxtilt_x": float(dxtilt_x),
                "dxtilt_y": float(dxtilt_y),
                "xtilt_x": float(xtilt_x),
                "xtilt_y": float(xtilt_y),
                "defocus_um": ctf["defocus_um"],
                "astig_um": ctf["astig_um"],
                "astig_angle_deg": ctf["astig_angle_deg"],
                "astig_x_um": float(astig_x) if np.isfinite(astig_x) else np.nan,
                "astig_y_um": float(astig_y) if np.isfinite(astig_y) else np.nan,
                "phase_shift_deg": ctf["phase_shift_deg"],
                "fit_score": ctf["fit_score"],
                "resolution_A": ctf["resolution_A"],
            }
            rows.append(row)
            echo(
                f"  defocus={ctf['defocus_um']:.4f} um, "
                f"astig={ctf['astig_um']:.4f} um @ {ctf['astig_angle_deg']:.1f} deg, "
                f"res={ctf['resolution_A']:.2f} A"
            )

        defocus_fit = fit_plane(
            [r["dxtilt_x"] for r in rows],
            [r["dxtilt_y"] for r in rows],
            [r["defocus_um"] for r in rows],
        )
        astig_fit = fit_plane(
            [r["dxtilt_x"] for r in rows],
            [r["dxtilt_y"] for r in rows],
            [r["astig_um"] for r in rows],
        )
        astig_x_fit = fit_plane(
            [r["dxtilt_x"] for r in rows],
            [r["dxtilt_y"] for r in rows],
            [r["astig_x_um"] for r in rows],
        )
        astig_y_fit = fit_plane(
            [r["dxtilt_x"] for r in rows],
            [r["dxtilt_y"] for r in rows],
            [r["astig_y_um"] for r in rows],
        )

        example_dx = xtilt_step
        example_dy = 0.0
        calibration = {
            "model": "xtilt_defocus_astigmatism_plane",
            "created": datetime.now().isoformat(timespec="seconds"),
            "start_xtilt_x": float(start_xtilt_x),
            "start_xtilt_y": float(start_xtilt_y),
            "grid_size": int(grid_size),
            "xtilt_step": float(xtilt_step),
            "settle_delay_s": float(settle_delay_s),
            "defocus_um": defocus_fit,
            "astig_um": astig_fit,
            "astig_x_um": astig_x_fit,
            "astig_y_um": astig_y_fit,
            "formula": (
                "delta_z = ax * dxtilt_x + ay * dxtilt_y  "
                "(dxtilt relative to start_xtilt)"
            ),
            "example": {
                "dxtilt_x": example_dx,
                "dxtilt_y": example_dy,
                "delta_defocus_um": predict_delta(defocus_fit, example_dx, example_dy),
                "delta_astig_um": predict_delta(astig_fit, example_dx, example_dy),
            },
        }

        with open(csv_measurements, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

        with open(calibration_json, "w") as fh:
            json.dump(calibration, fh, indent=2)

        make_plot(plot_file, rows, defocus_fit, astig_fit)

        echo("================================================")
        echo("FITS  (z = a0 + ax * dX + ay * dY)")
        echo(
            f"  defocus: a0={defocus_fit['a0']:.6f} um, "
            f"ax={defocus_fit['ax']:.6f} um/xtilt, "
            f"ay={defocus_fit['ay']:.6f} um/xtilt, "
            f"|grad|={defocus_fit['grad']:.6f}, RMS={defocus_fit['rms']:.6f} um"
        )
        echo(
            f"  astig:   a0={astig_fit['a0']:.6f} um, "
            f"ax={astig_fit['ax']:.6f} um/xtilt, "
            f"ay={astig_fit['ay']:.6f} um/xtilt, "
            f"|grad|={astig_fit['grad']:.6f}, RMS={astig_fit['rms']:.6f} um"
        )
        echo(
            f"  For dX-tilt={example_dx:g}, dY-tilt=0: "
            f"delta defocus={calibration['example']['delta_defocus_um']:.6f} um, "
            f"delta astig={calibration['example']['delta_astig_um']:.6f} um"
        )
        echo("================================================")
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved calibration JSON: {calibration_json}")
        echo(f"Saved plot: {plot_file}")

    finally:
        sem.SetXLensDeflector(
            xtilt_lens_index, float(original_xtilt[0]), float(original_xtilt[1])
        )
        echo(
            f"Restored XLensDeflector({xtilt_lens_index}) to "
            f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
        )

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
