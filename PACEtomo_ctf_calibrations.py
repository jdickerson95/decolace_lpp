#!Python
# ===================================================================
# ScriptName     PACEtomo_ctf_calibrations
# Purpose:       Shared CTF X-tilt back-projection, defocus-error
#                ChangeFocus scaling, and stigmator-astigmatism map.
# ===================================================================

import json
import os

import numpy as np
import serialem as sem

_sem = sem
_logger = None

hasXLens = True
useCtfXtilt = True
ctfXtiltX = 0.0
ctfXtiltY = 0.0
xtilt_lens_index = 2

xtilt_calibration_file = ""
defocus_error_file = ""
astig_calibration_file = ""

_xtilt_cal = None
_defocus_error_cal = None
_astig_cal = None


def _echo(text):
    if _logger is not None:
        _logger(text)
    else:
        _sem.Echo(text)


def configure(
    sem_module=None,
    logger=None,
    has_x_lens=None,
    use_ctf_xtilt=None,
    ctf_xtilt_x=None,
    ctf_xtilt_y=None,
    xtilt_lens_index_value=None,
    xtilt_calibration_path=None,
    defocus_error_path=None,
    astig_calibration_path=None,
):
    """Set SerialEM module, flags, and optional calibration JSON paths."""
    global _sem, _logger, hasXLens, useCtfXtilt, ctfXtiltX, ctfXtiltY
    global xtilt_lens_index, xtilt_calibration_file, defocus_error_file
    global astig_calibration_file, _xtilt_cal, _defocus_error_cal, _astig_cal
    if sem_module is not None:
        _sem = sem_module
    _logger = logger
    if has_x_lens is not None:
        hasXLens = bool(has_x_lens)
    if use_ctf_xtilt is not None:
        useCtfXtilt = bool(use_ctf_xtilt)
    if ctf_xtilt_x is not None:
        ctfXtiltX = float(ctf_xtilt_x)
    if ctf_xtilt_y is not None:
        ctfXtiltY = float(ctf_xtilt_y)
    if xtilt_lens_index_value is not None:
        xtilt_lens_index = int(xtilt_lens_index_value)
    if xtilt_calibration_path is not None:
        xtilt_calibration_file = xtilt_calibration_path or ""
        _xtilt_cal = None
    if defocus_error_path is not None:
        defocus_error_file = defocus_error_path or ""
        _defocus_error_cal = None
    if astig_calibration_path is not None:
        astig_calibration_file = astig_calibration_path or ""
        _astig_cal = None


def _load_json(path):
    if not path:
        return None
    path = os.path.abspath(path)
    with open(path, "r") as fh:
        return json.load(fh)


def get_xtilt_calibration():
    global _xtilt_cal
    if _xtilt_cal is None and xtilt_calibration_file:
        _xtilt_cal = _load_json(xtilt_calibration_file)
    return _xtilt_cal


def get_defocus_error_calibration():
    global _defocus_error_cal
    if _defocus_error_cal is None and defocus_error_file:
        _defocus_error_cal = _load_json(defocus_error_file)
    return _defocus_error_cal


def get_astig_calibration():
    global _astig_cal
    if _astig_cal is None and astig_calibration_file:
        _astig_cal = _load_json(astig_calibration_file)
    return _astig_cal


def parse_ctffind(cfind):
    """
    CtfFind reportedValue1-6:
      mean defocus [um], astigmatism [um], angle [deg],
      extra phase shift [deg], fitting score, resolution [A].
    """
    vals = [float(v) for v in cfind] if cfind is not None else []
    while len(vals) < 6:
        vals.append(np.nan)
    return {
        "defocus_um": vals[0],
        "astig_um": vals[1],
        "astig_angle_deg": vals[2],
        "phase_shift_deg": vals[3],
        "fit_score": vals[4],
        "resolution_A": vals[5],
    }


def astig_components(astig_um, angle_deg):
    """2-theta astigmatism components [um]."""
    astig_um = float(astig_um)
    angle_deg = float(angle_deg)
    if not np.isfinite(astig_um) or not np.isfinite(angle_deg):
        return np.nan, np.nan
    theta = np.deg2rad(2.0 * angle_deg)
    return astig_um * np.cos(theta), astig_um * np.sin(theta)


def astig_from_components(astig_x_um, astig_y_um):
    """Magnitude [um] and angle [deg] from 2-theta components."""
    ax = float(astig_x_um)
    ay = float(astig_y_um)
    mag = float(np.hypot(ax, ay))
    if mag < 1e-12:
        return 0.0, 0.0
    angle = 0.5 * np.degrees(np.arctan2(ay, ax))
    return mag, float(angle)


def _plane_ax_ay(cal, key):
    if not cal or key not in cal:
        return np.nan, np.nan
    block = cal[key]
    if not isinstance(block, dict):
        return np.nan, np.nan
    return float(block.get("ax", np.nan)), float(block.get("ay", np.nan))


def should_use_ctf_xtilt():
    return bool(useCtfXtilt) and bool(hasXLens)


def report_xtilt():
    if not hasXLens:
        return float("nan"), float("nan")
    xtilt = _sem.ReportXLensDeflector(xtilt_lens_index)
    return float(xtilt[0]), float(xtilt[1])


def set_xtilt(x, y):
    if not hasXLens:
        return
    _sem.SetXLensDeflector(xtilt_lens_index, float(x), float(y))


def ctf_xtilt_target():
    return float(ctfXtiltX), float(ctfXtiltY)


def correct_ctf_to_working_xtilt(
    defocus_um,
    astig_x_um,
    astig_y_um,
    working_xt,
    ctf_xt,
):
    """Map CtfFind values from CTF X-tilt back to the working X-tilt.

    value_working = value_measured + ax * (working_x - ctf_x)
                                 + ay * (working_y - ctf_y)
    No-op if CTF X-tilt is disabled or no xtilt JSON is loaded.
    """
    defocus_um = float(defocus_um) if defocus_um is not None else np.nan
    astig_x_um = float(astig_x_um) if astig_x_um is not None else np.nan
    astig_y_um = float(astig_y_um) if astig_y_um is not None else np.nan
    if not should_use_ctf_xtilt():
        return defocus_um, astig_x_um, astig_y_um
    cal = get_xtilt_calibration()
    if not cal:
        return defocus_um, astig_x_um, astig_y_um
    wx, wy = float(working_xt[0]), float(working_xt[1])
    cx, cy = float(ctf_xt[0]), float(ctf_xt[1])
    dx = wx - cx
    dy = wy - cy
    if abs(dx) < 1e-12 and abs(dy) < 1e-12:
        return defocus_um, astig_x_um, astig_y_um

    def _shift(value, key):
        ax, ay = _plane_ax_ay(cal, key)
        if not np.isfinite(value) or not np.isfinite(ax) or not np.isfinite(ay):
            return value
        return float(value + ax * dx + ay * dy)

    defocus_um = _shift(defocus_um, "defocus_um")
    astig_x_um = _shift(astig_x_um, "astig_x_um")
    astig_y_um = _shift(astig_y_um, "astig_y_um")
    return defocus_um, astig_x_um, astig_y_um


def change_focus_command_um(desired_delta_um):
    """Commanded ChangeFocus [um] for a desired CTF defocus change."""
    delta = float(desired_delta_um)
    cal = get_defocus_error_calibration()
    slope = 1.0
    if cal is not None:
        if "slope" in cal and np.isfinite(float(cal["slope"])) and abs(float(cal["slope"])) > 1e-8:
            slope = float(cal["slope"])
        elif "command_scale" in cal and abs(float(cal.get("command_scale", 0))) > 1e-8:
            slope = 1.0 / float(cal["command_scale"])
    if not np.isfinite(slope) or abs(slope) < 1e-8:
        slope = 1.0
    return float(delta / slope)


def change_focus_for_desired_delta(desired_delta_um):
    """ChangeFocus scaled by defocus-error slope (1 if no calibration)."""
    command = change_focus_command_um(desired_delta_um)
    _sem.ChangeFocus(command)
    return command


def stig_delta_to_cancel_astig(astig_x_um, astig_y_um):
    """Objective-stigmator delta (SerialEM units) to cancel 2θ astig [um]."""
    cal = get_astig_calibration()
    if not cal:
        raise ValueError("No astigmatism calibration loaded (astig_calibration_file).")
    inv_m = np.array(cal["inv_M"], dtype=float)
    vec = np.array([float(astig_x_um), float(astig_y_um)], dtype=float)
    if not np.all(np.isfinite(vec)):
        return np.nan, np.nan
    delta = -inv_m @ vec
    return float(delta[0]), float(delta[1])


def fit_stig_map(dstig_x, dstig_y, astig_x, astig_y):
    """Least-squares [astig_x, astig_y] = M @ [dStig_x, dStig_y] + b."""
    dsx = np.asarray(dstig_x, dtype=float).ravel()
    dsy = np.asarray(dstig_y, dtype=float).ravel()
    ax = np.asarray(astig_x, dtype=float).ravel()
    ay = np.asarray(astig_y, dtype=float).ravel()
    mask = np.isfinite(dsx) & np.isfinite(dsy) & np.isfinite(ax) & np.isfinite(ay)
    n = int(np.sum(mask))
    if n < 3:
        raise ValueError(f"Need at least 3 finite stig/astig points (got {n}).")
    A = np.column_stack((dsx[mask], dsy[mask], np.ones(n)))
    cx, _, _, _ = np.linalg.lstsq(A, ax[mask], rcond=None)
    cy, _, _, _ = np.linalg.lstsq(A, ay[mask], rcond=None)
    M = np.array([[cx[0], cx[1]], [cy[0], cy[1]]], dtype=float)
    b = np.array([cx[2], cy[2]], dtype=float)
    det = float(np.linalg.det(M))
    if abs(det) < 1e-12:
        raise ValueError("Stigmator-to-astigmatism matrix is singular; increase stig_step.")
    inv_M = np.linalg.inv(M)
    pred = (M @ np.vstack((dsx[mask], dsy[mask]))).T + b
    meas = np.column_stack((ax[mask], ay[mask]))
    rms = float(np.sqrt(np.mean(np.sum((pred - meas) ** 2, axis=1))))
    return {
        "M": M.tolist(),
        "inv_M": inv_M.tolist(),
        "b": b.tolist(),
        "rms_um": rms,
        "n": n,
        "det": det,
        "formula": (
            "[astig_x, astig_y] = M @ [dStig_x, dStig_y] + b; "
            "dStig = -inv(M) @ astig_vec to cancel measured astig"
        ),
    }


def acquire_ctf(
    defocus_lo,
    defocus_hi,
    shot="F",
    max_attempts=3,
    resolution_max_A=20.0,
    retry_delay_s=5.0,
):
    """Optional CTF X-tilt, CtfFind, restore X-tilt, back-project to working X-tilt.

    Returns dict with raw and working-X-tilt defocus/astig, plus resolution.
    """
    working_xt = report_xtilt()
    moved = False
    ctf_xt = working_xt
    try:
        if should_use_ctf_xtilt():
            ctf_xt = ctf_xtilt_target()
            set_xtilt(ctf_xt[0], ctf_xt[1])
            moved = True
        cfind = []
        _sem.NoMessageBoxOnError(1)
        try:
            for attempt in range(1, int(max_attempts) + 1):
                if attempt > 1:
                    _sem.Delay(float(retry_delay_s), "s")
                if str(shot).upper() == "L":
                    _sem.L()
                else:
                    _sem.F()
                cfind = _sem.CtfFind("A", float(defocus_lo), float(defocus_hi))
                if len(cfind) == 0:
                    _echo(f"ERROR: CtfFind failed on attempt {attempt}/{max_attempts}.")
                    if attempt < int(max_attempts):
                        continue
                    parsed = parse_ctffind([])
                    parsed.update(
                        {
                            "astig_x_um": np.nan,
                            "astig_y_um": np.nan,
                            "defocus_raw_um": np.nan,
                            "astig_x_raw_um": np.nan,
                            "astig_y_raw_um": np.nan,
                            "working_xtilt_x": working_xt[0],
                            "working_xtilt_y": working_xt[1],
                            "ctf_xtilt_x": ctf_xt[0],
                            "ctf_xtilt_y": ctf_xt[1],
                            "used_ctf_xtilt": moved,
                        }
                    )
                    return parsed
                parsed = parse_ctffind(cfind)
                if np.isfinite(parsed["resolution_A"]) and parsed["resolution_A"] <= float(resolution_max_A):
                    break
                _echo(
                    f"WARNING: CtfFind resolution {parsed['resolution_A']:.2f} A > "
                    f"{resolution_max_A} A (attempt {attempt}/{max_attempts})"
                )
            else:
                _echo(
                    f"WARNING: CtfFind resolution still > {resolution_max_A} A after "
                    f"{max_attempts} attempts; using last result."
                )
        finally:
            _sem.NoMessageBoxOnError(0)
        parsed = parse_ctffind(cfind)
        ax, ay = astig_components(parsed["astig_um"], parsed["astig_angle_deg"])
        parsed["astig_x_um"] = float(ax) if np.isfinite(ax) else np.nan
        parsed["astig_y_um"] = float(ay) if np.isfinite(ay) else np.nan
        parsed["defocus_raw_um"] = parsed["defocus_um"]
        parsed["astig_x_raw_um"] = parsed["astig_x_um"]
        parsed["astig_y_raw_um"] = parsed["astig_y_um"]
        d_w, ax_w, ay_w = correct_ctf_to_working_xtilt(
            parsed["defocus_um"],
            parsed["astig_x_um"],
            parsed["astig_y_um"],
            working_xt,
            ctf_xt,
        )
        parsed["defocus_um"] = d_w
        parsed["astig_x_um"] = ax_w
        parsed["astig_y_um"] = ay_w
        mag, ang = astig_from_components(ax_w, ay_w)
        if np.isfinite(mag):
            parsed["astig_um"] = mag
            parsed["astig_angle_deg"] = ang
        parsed["working_xtilt_x"] = working_xt[0]
        parsed["working_xtilt_y"] = working_xt[1]
        parsed["ctf_xtilt_x"] = ctf_xt[0]
        parsed["ctf_xtilt_y"] = ctf_xt[1]
        parsed["used_ctf_xtilt"] = moved
        _echo(
            f"CtfFind: {parsed['defocus_raw_um']:.4f} um"
            + (
                f" -> {parsed['defocus_um']:.4f} um at working X-tilt"
                if moved and np.isfinite(parsed["defocus_um"])
                else ""
            )
            + f" ({parsed['resolution_A']:.2f} A)"
        )
        return parsed
    finally:
        if moved:
            set_xtilt(working_xt[0], working_xt[1])
