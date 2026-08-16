#!Python
# ===================================================================
# ScriptName     calibrate_astigmatism
# Purpose:       Sweep objective stigmator on a grid, measure astigmatism
#                with CtfFind (optional CTF X-tilt, back-projected to the
#                working X-tilt), and fit
#                  [astig_x, astig_y] = M @ [dStig_x, dStig_y] + b
#                so a measured astig vector can be cancelled with
#                  dStig = -inv(M) @ astig_vec
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

# Odd grid size: 3, 5, or 7.
grid_size = 5

# Step between neighboring stigmator points (SerialEM units, typically -1..1).
stig_step = 0.02

# Settle time after SetObjectiveStigmator [s].
settle_delay_s = 1.0

# Optional CTF X-tilt (move off laser). Disable for a simpler on-axis measure.
useCtfXtilt = True
hasXLens = True
ctfXtiltX = 0.002836
ctfXtiltY = 0.003867
xtilt_lens_index = 2
xtilt_calibration_file = ""  # JSON from check_xtilt_defoc_astig.py; empty = skip back-project

# CTF search range [microns].
ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2
ctf_resolution_max_A = 20.0
ctf_max_attempts = 3
ctf_retry_delay_s = 5.0

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""

csv_measurements = "astigmatism_measurements.csv"
calibration_json = "astigmatism_calibration.json"
plot_file = "astigmatism_fit.png"

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


def make_plot(path, rows, fit):
    dsx = np.array([r["dstig_x"] for r in rows], dtype=float)
    dsy = np.array([r["dstig_y"] for r in rows], dtype=float)
    ax = np.array([r["astig_x_um"] for r in rows], dtype=float)
    ay = np.array([r["astig_y_um"] for r in rows], dtype=float)
    M = np.array(fit["M"], dtype=float)
    b = np.array(fit["b"], dtype=float)
    pred = (M @ np.vstack((dsx, dsy))).T + b
    pred_x, pred_y = pred[:, 0], pred[:, 1]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)
    axes[0].scatter(pred_x, ax, label="astig_x")
    axes[0].scatter(pred_y, ay, marker="x", label="astig_y")
    lo = np.nanmin([pred_x, pred_y, ax, ay])
    hi = np.nanmax([pred_x, pred_y, ax, ay])
    axes[0].plot([lo, hi], [lo, hi], "--", color="gray")
    axes[0].set_xlabel("Fitted astig component (um)")
    axes[0].set_ylabel("Measured astig component (um)")
    axes[0].set_title("Stigmator map fit")
    axes[0].legend(fontsize=8)

    res_x = ax - pred_x
    res_y = ay - pred_y
    axes[1].scatter(dsx, res_x, label="astig_x residual")
    axes[1].scatter(dsy, res_y, marker="x", label="astig_y residual")
    axes[1].axhline(0, color="gray", linestyle="--", linewidth=1)
    axes[1].set_xlabel("dStig (SerialEM units)")
    axes[1].set_ylabel("Residual (um)")
    axes[1].set_title("Fit residual")
    axes[1].legend(fontsize=8)

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    prepare_output_paths()
    ctfcal.configure(
        sem_module=sem,
        logger=echo,
        has_x_lens=hasXLens,
        use_ctf_xtilt=useCtfXtilt,
        ctf_xtilt_x=ctfXtiltX,
        ctf_xtilt_y=ctfXtiltY,
        xtilt_lens_index_value=xtilt_lens_index,
        xtilt_calibration_path=xtilt_calibration_file,
    )
    offsets = grid_offsets(grid_size, stig_step)
    echo("##### Objective stigmator astigmatism calibration #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")
    echo(
        f"grid={grid_size}x{grid_size}, step={stig_step}, n={len(offsets)}, "
        f"useCtfXtilt={useCtfXtilt}"
    )

    orig_stig = sem.ReportObjectiveStigmator()
    start_sx, start_sy = float(orig_stig[0]), float(orig_stig[1])
    echo(f"Saved ObjectiveStigmator: ({start_sx:.6f}, {start_sy:.6f})")

    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)

    rows = []
    try:
        for i, (dsx, dsy) in enumerate(offsets, start=1):
            sx = clip_stig(start_sx + dsx)
            sy = clip_stig(start_sy + dsy)
            echo("------------------------------------------------")
            echo(
                f"{i}/{len(offsets)} stig ({sx:.4f}, {sy:.4f})  "
                f"d=({dsx:.4f}, {dsy:.4f})"
            )
            sem.SetObjectiveStigmator(sx, sy)
            sem.Delay(settle_delay_s, "s")
            ctf = ctfcal.acquire_ctf(
                ctf_defocus_lo,
                ctf_defocus_hi,
                shot="L",
                max_attempts=ctf_max_attempts,
                resolution_max_A=ctf_resolution_max_A,
                retry_delay_s=ctf_retry_delay_s,
            )
            row = {
                "dstig_x": float(dsx),
                "dstig_y": float(dsy),
                "stig_x": float(sx),
                "stig_y": float(sy),
                "defocus_um": ctf["defocus_um"],
                "astig_um": ctf["astig_um"],
                "astig_angle_deg": ctf["astig_angle_deg"],
                "astig_x_um": ctf.get("astig_x_um", np.nan),
                "astig_y_um": ctf.get("astig_y_um", np.nan),
                "resolution_A": ctf["resolution_A"],
                "used_ctf_xtilt": bool(ctf.get("used_ctf_xtilt", False)),
            }
            rows.append(row)
            echo(
                f"  astig={row['astig_um']:.4f} um "
                f"(x={row['astig_x_um']:.4f}, y={row['astig_y_um']:.4f})"
            )
        fit = ctfcal.fit_stig_map(
            [r["dstig_x"] for r in rows],
            [r["dstig_y"] for r in rows],
            [r["astig_x_um"] for r in rows],
            [r["astig_y_um"] for r in rows],
        )
        calibration = {
            "model": "objective_stigmator_astigmatism",
            "created": datetime.now().isoformat(timespec="seconds"),
            "grid_size": int(grid_size),
            "stig_step": float(stig_step),
            "settle_delay_s": float(settle_delay_s),
            "start_stig_x": start_sx,
            "start_stig_y": start_sy,
            "useCtfXtilt": bool(useCtfXtilt),
            "ctfXtiltX": float(ctfXtiltX) if useCtfXtilt else None,
            "ctfXtiltY": float(ctfXtiltY) if useCtfXtilt else None,
            "xtilt_calibration_file": xtilt_calibration_file or None,
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
        echo(f"  RMS {fit['rms_um']:.6f} um  (n={fit['n']})")
        echo("  To cancel astig at working X-tilt: dStig = -inv(M) @ [astig_x, astig_y]")
        echo("================================================")
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved calibration JSON: {calibration_json}")
        echo(f"Saved plot: {plot_file}")
    finally:
        sem.SetObjectiveStigmator(start_sx, start_sy)
        echo(f"Restored ObjectiveStigmator to ({start_sx:.6f}, {start_sy:.6f})")

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
