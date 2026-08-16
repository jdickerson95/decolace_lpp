#!Python
# ===================================================================
# ScriptName     calibrate_defocus_error
# Purpose:       Sweep ChangeFocus in equal increments, wait for settle,
#                measure defocus with CtfFind, and fit
#                  measured_um = intercept_um + slope * commanded_change_um
#                Use slope to command the ChangeFocus needed for a desired
#                CTF defocus change:
#                  ChangeFocus(desired_delta_um / slope)
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

############ SETTINGS ############

# Commanded ChangeFocus increment [um] at each step (negative = more underfocus).
increment_um = -0.5

# Number of increments after the starting measurement (total points = n_increments + 1).
n_increments = 8

# Settle time after ChangeFocus [s].
settle_delay_s = 1.0

# CTF search range [microns].
ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2

# X-tilt for CtfFind (match PACEtomo ctfXtiltX/Y). Restored after run.
use_ctf_xtilt = True
xtilt_lens_index = 2
ctf_xtilt_x = 0.002836
ctf_xtilt_y = 0.003867

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""

csv_measurements = "defocus_error_measurements.csv"
calibration_json = "defocus_error_calibration.json"
plot_file = "defocus_error_fit.png"

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


def run_ctffind():
    sem.NoMessageBoxOnError(1)
    try:
        cfind = sem.CtfFind("A", ctf_defocus_lo, ctf_defocus_hi)
    finally:
        sem.NoMessageBoxOnError(0)
    if len(cfind) == 0:
        echo("ERROR: CtfFind failed.")
        sem.Exit()
    return float(cfind[0]), float(cfind[-1])


def acquire_ctf():
    sem.L()
    return run_ctffind()


def report_defocus_um():
    return float(sem.ReportDefocus())


def fit_line(commanded_um, measured_um):
    x = np.array(commanded_um, dtype=float)
    y = np.array(measured_um, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.sum(mask) < 2:
        return {
            "intercept_um": np.nan,
            "slope": np.nan,
            "rms_um": np.nan,
        }
    slope, intercept = np.polyfit(x[mask], y[mask], 1)
    pred = slope * x[mask] + intercept
    return {
        "intercept_um": float(intercept),
        "slope": float(slope),
        "rms_um": float(np.sqrt(np.mean((pred - y[mask]) ** 2))),
    }


def make_plot(path, rows, fit):
    commanded = np.array([float(r["commanded_change_um"]) for r in rows], dtype=float)
    measured = np.array([float(r["ctf_defocus_um"]) for r in rows], dtype=float)
    x_line = np.linspace(np.min(commanded), np.max(commanded), 100)
    y_line = fit["intercept_um"] + fit["slope"] * x_line
    residual = measured - (fit["intercept_um"] + fit["slope"] * commanded)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)

    ax = axes[0]
    ax.scatter(commanded, measured, label="CtfFind")
    ax.plot(x_line, y_line, color="C1", label="fit")
    ax.plot(x_line, fit["intercept_um"] + x_line, "--", color="gray", label="ideal slope=1")
    ax.set_xlabel("Commanded ChangeFocus from start (um)")
    ax.set_ylabel("CTF defocus (um)")
    ax.set_title("Measured vs commanded defocus")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.scatter(commanded, residual, label="residual")
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("Commanded ChangeFocus from start (um)")
    ax.set_ylabel("Residual (um)")
    ax.set_title("Fit residual")
    ax.legend(fontsize=8)

    fig.savefig(path, dpi=150)
    plt.show()


def main():
    sem.SuppressReports()
    prepare_output_paths()
    echo("##### Defocus error calibration #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Output directory: {os.path.dirname(csv_measurements)}")
    echo(
        f"increment={increment_um:.4f} um, n_increments={n_increments}, "
        f"settle={settle_delay_s:.3f} s"
    )

    original_xtilt = None
    if use_ctf_xtilt:
        original_xtilt = sem.ReportXLensDeflector(xtilt_lens_index)
        echo(
            f"Saved XLensDeflector({xtilt_lens_index}): "
            f"({float(original_xtilt[0]):.6f}, {float(original_xtilt[1]):.6f})"
        )
        sem.SetXLensDeflector(xtilt_lens_index, float(ctf_xtilt_x), float(ctf_xtilt_y))
        echo(f"CTF X-tilt ({ctf_xtilt_x:.6f}, {ctf_xtilt_y:.6f})")

    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)

    commanded_change = 0.0
    rows = []
    try:
        echo("------------------------------------------------")
        echo("Starting measurement (no ChangeFocus yet)")
        start_report = report_defocus_um()
        ctf_defocus, ctf_res = acquire_ctf()
        row = {
            "step": 0,
            "commanded_change_um": commanded_change,
            "report_defocus_um": start_report,
            "ctf_defocus_um": ctf_defocus,
            "ctf_resolution_A": ctf_res,
        }
        rows.append(row)
        echo(
            f"step 0: commanded=0.000 um, ReportDefocus={start_report:.4f} um, "
            f"CTF={ctf_defocus:.4f} um"
        )

        for step in range(1, n_increments + 1):
            echo("------------------------------------------------")
            echo(f"ChangeFocus({increment_um:.4f})")
            sem.ChangeFocus(increment_um)
            commanded_change += increment_um
            sem.Delay(settle_delay_s, "s")
            report = report_defocus_um()
            ctf_defocus, ctf_res = acquire_ctf()
            row = {
                "step": step,
                "commanded_change_um": commanded_change,
                "report_defocus_um": report,
                "ctf_defocus_um": ctf_defocus,
                "ctf_resolution_A": ctf_res,
            }
            rows.append(row)
            echo(
                f"step {step}: commanded={commanded_change:.4f} um, "
                f"ReportDefocus={report:.4f} um, CTF={ctf_defocus:.4f} um"
            )

        fit = fit_line(
            [r["commanded_change_um"] for r in rows],
            [r["ctf_defocus_um"] for r in rows],
        )
        slope = fit["slope"]
        intercept = fit["intercept_um"]
        command_scale = float("nan") if not np.isfinite(slope) or slope == 0 else 1.0 / slope

        calibration = {
            "model": "defocus_error",
            "created": datetime.now().isoformat(timespec="seconds"),
            "increment_um": float(increment_um),
            "n_increments": int(n_increments),
            "settle_delay_s": float(settle_delay_s),
            "use_ctf_xtilt": bool(use_ctf_xtilt),
            "ctf_xtilt_x": float(ctf_xtilt_x) if use_ctf_xtilt else None,
            "ctf_xtilt_y": float(ctf_xtilt_y) if use_ctf_xtilt else None,
            "intercept_um": intercept,
            "slope": slope,
            "rms_um": fit["rms_um"],
            "command_scale": command_scale,
            "formula": (
                "measured_um = intercept_um + slope * commanded_change_um; "
                "ChangeFocus(desired_delta_um / slope) to change CTF defocus by "
                "desired_delta_um"
            ),
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
        echo(
            f"  measured = {intercept:.6f} + {slope:.6f} * commanded_change  "
            f"(RMS {fit['rms_um']:.6f} um)"
        )
        echo(
            f"  To change CTF defocus by D um: ChangeFocus(D / {slope:.6f}) "
            f"= {command_scale:.6f} * D"
        )
        echo("================================================")
        echo(f"Saved measurements: {csv_measurements}")
        echo(f"Saved calibration JSON: {calibration_json}")
        echo(f"Saved plot: {plot_file}")

    finally:
        if commanded_change != 0.0:
            sem.ChangeFocus(-commanded_change)
            echo(f"Restored focus with ChangeFocus({-commanded_change:.4f})")
        if original_xtilt is not None:
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
