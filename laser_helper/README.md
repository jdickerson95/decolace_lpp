# laser_helper

SerialEM scripts for the session ronchi numbers you paste into collect: in-plane C3, `ronchiCorrectKs`, and `ronchiTargetPhaseA` / `B`.

Run **after** the CTF calibrations (`calibrate_defocus_error`, `check_xtilt_defoc_astig`). Trial must match Record (same beam position). Copy `ronchiPixelSize` / `ronchiBinning` / `ronchiPeakRadius` and `ctfXtiltX` / `ctfXtiltY` from [`decolace_LPP/decolace_collect.py`](../decolace_LPP/decolace_collect.py) if you changed them.

Order: **on plane, then on peak.**

## `auto_on_plane.py`

Start **roughly on plane**. The script:

1. Reads the current `ImageDistanceOffset` (C3).
2. Moves `+20`, then `+5` steps to `+50`, measuring ronchi fringe spacing at each point.
3. Repeats the same sweep in the negative direction.
4. Takes the magnitude of each fringe spacing, then assigns a negative sign to the negative-C3 branch.
5. Fits a line and reports the C3 where signed spacing is 0 (in plane).
6. Goes to that C3 **−20**, measures `ks`, and prints `ronchiCorrectKs` to copy.

### Where to paste

From the SerialEM log:

- `ronchiCorrectKs = [[...], [...]]`
- `ronchiC3Offset = -20`

into **both** [`decolace_LPP/decolace_collect.py`](../decolace_LPP/decolace_collect.py) and [`PACEtomo.py`](../PACEtomo.py) (`########## Ronchigram / laser alignment ##########`).

Set the session `ImageDistanceOffset` to the reported in-plane C3. Collect then adds `ronchiC3Offset` only for the Trial shot.

Writes `auto_on_plane_measurements.csv` and `auto_on_plane_fit.png`.

## How to measure 1-fringe X-tilt (needed for `auto_on_peak.py`)

`auto_on_peak.py` does **not** guess the X-tilt period. You type `fringe_xtilt_x` and `fringe_xtilt_y`: the XLensDeflector change that shifts the ronchigram by **one fringe** (one full period) along X and along Y.

Do this after `auto_on_plane.py`, with C3 on plane and a Trial at **−20**:

1. `GoToLowDoseArea T`, `SetImageDistanceOffset` to plane−20, take a Trial. You should see a clear two-direction fringe pattern (or two FFT spots).
2. Note the current X-tilt (`ReportXLensDeflector 2`) and the look of the pattern (or the FFT peak phases in the log if you already ran a ronchi analyze).
3. Change **only X-tilt X** in small steps. Watch one set of fringes walk across the field. Stop when the pattern looks the same as the start — that is one period.  
   `fringe_xtilt_x = (X after) − (X before)`.
4. Restore X-tilt X. Repeat changing **only X-tilt Y**.  
   `fringe_xtilt_y = (Y after) − (Y before)`.
5. Use the **absolute** deltas (sign does not matter; the 5×5 grid is centered and spans ±½ fringe).

FFT check: one fringe is a **2π** wrap of that fringe’s phase. If the vertical-fringe phase returns to the same value after a 2π jump, that X-tilt delta is `fringe_xtilt_x`.

Typical size is small (often a few 1e-4). If the 5×5 later looks flat or the quadratic peak sits on the edge, the 1-fringe values are too large or too small; remeasure.

Paste both numbers into the `auto_on_peak.py` settings block before you run it. The script exits if either is still `0`.

## `auto_on_peak.py`

Requires: CTF X-tilt values, 1-fringe X-tilt (above), and preferably the `check_xtilt_defoc_astig.py` JSON so defocus can be back-projected to the starting X-tilt.

Start at the **working** (laser) X-tilt, C3 on plane.

1. Moves to `ctfXtiltX` / `ctfXtiltY`, Records, CtfFind (full defocus range). That is a clean CTF off the laser.
2. Back-projects that defocus to the starting X-tilt if `xtilt_calibration_file` is set. Saves that number as the **fixed** defocus.
3. Restores the starting X-tilt.
4. 5×5 X-tilt grid centered on start, spanning **one fringe** in X and in Y (`±0.5 * fringe_xtilt_*`).
5. At each point: Record, CtfFind with min=max=saved defocus and a phase search (`phase_search_lo` to `phase_search_hi`, range less than 90 deg). Only phase is fitted.
6. Fits a quadratic to phase vs X-tilt and takes the **maximum**. If that peak is not a maximum inside the grid, uses the discrete grid max instead.
7. Sets that X-tilt, takes a Trial at C3 **−20**, and prints `ronchiTargetPhaseA` / `ronchiTargetPhaseB`.

Restores starting C3 and X-tilt at the end.

### Where to paste

From the SerialEM log:

- `ronchiTargetPhaseA = ...`
- `ronchiTargetPhaseB = ...`

into the same two files as `ronchiCorrectKs`. Also set the working `XLensDeflector` to the printed peak X-tilt.

Writes `auto_on_peak_measurements.csv`, `auto_on_peak.json`, and `auto_on_peak_fit.png`.
