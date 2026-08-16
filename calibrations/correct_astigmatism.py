#!Python
# ===================================================================
# ScriptName     correct_astigmatism
# Purpose:       Measure astigmatism with CtfFind (optional CTF X-tilt
#                off the laser, back-projected to the working X-tilt)
#                and apply the stigmator delta from
#                calibrate_astigmatism.py:
#                  dStig = -inv(M) @ [astig_x, astig_y]
# ===================================================================
import sys
sys.path.append(r"C:\Program Files\SerialEM\PythonModules")
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

import PACEtomo_ctf_calibrations as ctfcal

############ SETTINGS ############

# JSON from calibrate_astigmatism.py (required).
astig_calibration_file = r""

# JSON from check_xtilt_defoc_astig.py (used when useCtfXtilt is True).
xtilt_calibration_file = r""

# Measure CTF at ctfXtilt (off laser), then map astig back to working X-tilt.
useCtfXtilt = True
hasXLens = True
ctfXtiltX = 0.002836
ctfXtiltY = 0.003867
xtilt_lens_index = 2

settle_delay_s = 1.0
confirm_after = True

ctf_defocus_lo = -10.0
ctf_defocus_hi = -0.2
ctf_resolution_max_A = 20.0
ctf_max_attempts = 3
ctf_retry_delay_s = 5.0

########## END SETTINGS ##########


def echo(text):
    sem.Echo(text)


def measure():
    if settle_delay_s > 0:
        sem.Delay(settle_delay_s, "s")
    return ctfcal.acquire_ctf(
        ctf_defocus_lo,
        ctf_defocus_hi,
        shot="L",
        max_attempts=ctf_max_attempts,
        resolution_max_A=ctf_resolution_max_A,
        retry_delay_s=ctf_retry_delay_s,
    )


def main():
    sem.SuppressReports()
    if not astig_calibration_file:
        sem.OKBox("Set astig_calibration_file to the JSON from calibrate_astigmatism.py.")
        sem.Exit()
    ctfcal.configure(
        sem_module=sem,
        logger=echo,
        has_x_lens=hasXLens,
        use_ctf_xtilt=useCtfXtilt,
        ctf_xtilt_x=ctfXtiltX,
        ctf_xtilt_y=ctfXtiltY,
        xtilt_lens_index_value=xtilt_lens_index,
        xtilt_calibration_path=xtilt_calibration_file,
        astig_calibration_path=astig_calibration_file,
    )
    echo("##### Correct objective astigmatism #####")
    echo(f"Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    echo(f"Astig calibration: {os.path.abspath(astig_calibration_file)}")
    echo(f"useCtfXtilt={useCtfXtilt}")

    orig_stig = sem.ReportObjectiveStigmator()
    sx, sy = float(orig_stig[0]), float(orig_stig[1])
    echo(f"Current ObjectiveStigmator: ({sx:.6f}, {sy:.6f})")

    sem.GoToLowDoseArea("R")
    sem.SetImageShift(0, 0)

    ctf = measure()
    ax, ay = ctf.get("astig_x_um", np.nan), ctf.get("astig_y_um", np.nan)
    echo(
        f"Astig at working X-tilt: {ctf['astig_um']:.4f} um "
        f"(x={ax:.4f}, y={ay:.4f}) @ {ctf['astig_angle_deg']:.1f} deg"
    )
    if not np.isfinite(ax) or not np.isfinite(ay):
        echo("ERROR: CtfFind astigmatism is invalid; not changing stigmator.")
        sem.SuppressReports(0)
        sem.Exit()

    dsx, dsy = ctfcal.stig_delta_to_cancel_astig(ax, ay)
    echo(f"Stigmator delta: ({dsx:.6f}, {dsy:.6f})")
    new_sx = float(np.clip(sx + dsx, -1.0, 1.0))
    new_sy = float(np.clip(sy + dsy, -1.0, 1.0))
    if new_sx != sx + dsx or new_sy != sy + dsy:
        echo("WARNING: Stigmator change was clipped to [-1, 1].")
    sem.SetObjectiveStigmator(new_sx, new_sy)
    echo(f"Set ObjectiveStigmator to ({new_sx:.6f}, {new_sy:.6f})")

    if confirm_after:
        ctf2 = measure()
        echo(
            f"Astig after correction: {ctf2['astig_um']:.4f} um "
            f"(x={ctf2.get('astig_x_um', np.nan):.4f}, "
            f"y={ctf2.get('astig_y_um', np.nan):.4f})"
        )

    echo("Working X-tilt restored (CTF X-tilt not left applied).")
    sem.SuppressReports(0)
    sem.Exit()


if __name__ == "__main__":
    main()
