#!Python
# ===================================================================
# ScriptName     auto_c2_stig
# Purpose:       Measure a ronchigram, then move C2 (condenser)
#                stigmators so fringe spacing in x and y is the same,
#                using the map from calibrate_C2_astig.py:
#                  dC2 = inv(M) @ ([mean, mean] - [ks_x, ks_y])
# ===================================================================
import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
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

import numpy as np
import serialem as sem

import PACEtomo_ronchi as ronchi

############ SETTINGS ############

# JSON from calibrations/calibrate_C2_astig.py (required).
c2_astig_calibration_file = r""

# Trial ronchigram at C3 = current + working_offset.
# Must match collect/PACEtomo ronchiC3Offset (Trial is at plane −20).
working_offset = -20.0
ronchiDelay = 1.0
ronchiBinning = 32
ronchiPixelSize = 0.98e-4 * 2
ronchiPeakRadius = 100

settle_delay_s = 1.0

# Stop when |ks_x - ks_y| is below this (1/µm), or after max_iterations.
tolerance = 0.05
max_iterations = 2

# Clip each iteration's C2 change (SerialEM units).
max_stig_delta = 0.2

# Output directory. Empty string writes to the current SerialEM working directory.
save_dir = r""

result_json = "auto_c2_stig.json"

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


def prepare_output_paths():
    global result_json
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        result_json = os.path.join(save_dir, os.path.basename(result_json))
    result_json = os.path.abspath(result_json)


def clip_stig(value):
    return float(np.clip(value, -1.0, 1.0))


def report_c2():
    c2 = sem.ReportCondenserStigmator()
    return float(c2[0]), float(c2[1])


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
    """Fringe spatial frequencies (1/µm) of the two ronchi axes."""
    ks = np.asarray(ks, dtype=float)
    return float(np.linalg.norm(ks[0])), float(np.linalg.norm(ks[1]))


def format_ks(ks):
    ks = np.asarray(ks, dtype=float)
    return (
        f"[[{ks[0, 0]:.6f}, {ks[0, 1]:.6f}], "
        f"[{ks[1, 0]:.6f}, {ks[1, 1]:.6f}]]"
    )


def load_calibration(path):
    with open(path, "r") as fh:
        cal = json.load(fh)
    if "inv_M" not in cal:
        raise ValueError("Calibration JSON is missing inv_M.")
    inv_m = np.array(cal["inv_M"], dtype=float)
    if inv_m.shape != (2, 2) or not np.all(np.isfinite(inv_m)):
        raise ValueError("Calibration inv_M is not a finite 2x2 matrix.")
    return cal, inv_m


def stig_delta_to_equalize(inv_m, ks_x, ks_y):
    """C2 delta that moves both spacings to their mean, preserving mean |k|."""
    mean = 0.5 * (float(ks_x) + float(ks_y))
    delta_ks = np.array([mean - float(ks_x), mean - float(ks_y)], dtype=float)
    delta = inv_m @ delta_ks
    return float(delta[0]), float(delta[1]), float(mean)


def main():
    sem.SuppressReports()
    prepare_output_paths()
    if not c2_astig_calibration_file:
        sem.OKBox("Set c2_astig_calibration_file to the JSON from calibrate_C2_astig.py.")
        sem.Exit()
    cal_path = os.path.abspath(c2_astig_calibration_file)
    try:
        cal, inv_m = load_calibration(cal_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        echo(f"ERROR: could not load C2 calibration: {exc}")
        sem.OKBox(f"Could not load C2 calibration:\n{exc}")
        sem.Exit()

    echo("##### Laser helper: auto C2 stig #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"C2 calibration: {cal_path}")
    echo(f"working_offset={working_offset:g}, tolerance={tolerance:g} 1/um")

    start_c2 = report_c2()
    start_c3 = float(sem.ReportImageDistanceOffset())
    ronchi_c3 = start_c3 + working_offset
    echo(f"Starting CondenserStigmator: ({start_c2[0]:.6f}, {start_c2[1]:.6f})")
    echo(f"Starting C3: {start_c3:.6f}; Trial C3 = {ronchi_c3:.6f}")

    history = []
    applied = False
    try:
        sx, sy = start_c2
        final_ks_x = np.nan
        final_ks_y = np.nan
        final_ks = None
        for i in range(max_iterations + 1):
            echo("------------------------------------------------")
            echo(f"Measure {i + 1}: C2 ({sx:.6f}, {sy:.6f})")
            ks = acquire_ks(ronchi_c3)
            ks_x, ks_y = fringe_xy(ks)
            diff = ks_x - ks_y
            echo(
                f"  ks_x={ks_x:.4f}, ks_y={ks_y:.4f}, "
                f"|diff|={abs(diff):.4f}, ks={format_ks(ks)}"
            )
            history.append({
                "step": i,
                "c2_x": sx,
                "c2_y": sy,
                "ks_x": ks_x,
                "ks_y": ks_y,
                "ks_diff": diff,
                "ks": ks.tolist(),
            })
            final_ks_x, final_ks_y, final_ks = ks_x, ks_y, ks
            if abs(diff) <= tolerance:
                echo(f"|ks_x - ks_y| <= {tolerance:g}; C2 stig is equalized.")
                break
            if i >= max_iterations:
                echo(
                    f"WARNING: still |ks_x - ks_y|={abs(diff):.4f} "
                    f"after {max_iterations} correction(s)."
                )
                break

            dcx, dcy, mean = stig_delta_to_equalize(inv_m, ks_x, ks_y)
            echo(
                f"  target ks={mean:.4f} (both axes); "
                f"raw dC2=({dcx:.6f}, {dcy:.6f})"
            )
            step_norm = float(np.hypot(dcx, dcy))
            if step_norm > max_stig_delta:
                scale = max_stig_delta / step_norm
                dcx *= scale
                dcy *= scale
                echo(
                    f"  WARNING: clipped dC2 to max_stig_delta={max_stig_delta:g} "
                    f"→ ({dcx:.6f}, {dcy:.6f})"
                )
            new_sx = clip_stig(sx + dcx)
            new_sy = clip_stig(sy + dcy)
            if new_sx != sx + dcx or new_sy != sy + dcy:
                echo("  WARNING: C2 change was clipped to [-1, 1].")
            sem.SetCondenserStigmator(new_sx, new_sy)
            if settle_delay_s > 0:
                sem.Delay(settle_delay_s, "s")
            echo(f"  Set CondenserStigmator to ({new_sx:.6f}, {new_sy:.6f})")
            sx, sy = new_sx, new_sy
            applied = True

        echo("================================================")
        echo(
            f"Final C2 ({sx:.6f}, {sy:.6f}); "
            f"ks_x={final_ks_x:.4f}, ks_y={final_ks_y:.4f}, "
            f"|diff|={abs(final_ks_x - final_ks_y):.4f}"
        )
        if not applied:
            echo("Left CondenserStigmator unchanged.")
        echo("================================================")

        result = {
            "model": "auto_c2_stig",
            "created": datetime.now().isoformat(timespec="seconds"),
            "c2_astig_calibration_file": cal_path,
            "start_c2_x": start_c2[0],
            "start_c2_y": start_c2[1],
            "start_c3": start_c3,
            "ronchi_c3": float(ronchi_c3),
            "working_offset": float(working_offset),
            "tolerance": float(tolerance),
            "max_iterations": int(max_iterations),
            "applied": bool(applied),
            "final_c2_x": sx,
            "final_c2_y": sy,
            "final_ks_x": float(final_ks_x),
            "final_ks_y": float(final_ks_y),
            "final_ks_diff": float(final_ks_x - final_ks_y),
            "final_ks": None if final_ks is None else np.asarray(final_ks).tolist(),
            "history": history,
        }
        with open(result_json, "w") as fh:
            json.dump(result, fh, indent=2)
        echo(f"Saved result JSON: {result_json}")

    finally:
        sem.SetImageDistanceOffset(start_c3)
        left_c2 = report_c2()
        echo(f"Restored ImageDistanceOffset to {start_c3:.6f}")
        echo(f"Left CondenserStigmator at ({left_c2[0]:.6f}, {left_c2[1]:.6f})")

    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
