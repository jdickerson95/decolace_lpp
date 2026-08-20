# Ronchi / laser alignment — parameters you must update

decolace **realigns the laser at every hole** (`doRonchigram = True`). That only works if the FFT analysis knows **your** fringe spacing and **your** target phases.

The numbers checked into the repo (`ronchiCorrectKs`, `ronchiTargetPhaseA`, `ronchiTargetPhaseB`) are from **one** microscope session. They are not universal. **You must replace them with values from your own ronchigram calibration** before collecting.

## Where to edit (this is the important part)

Edit the **settings block of the script you actually run**. `PACEtomo_ronchi.py` has fallback defaults only; at startup the running script **overwrites** them via `ronchi.configure(...)`.

| If you are running | Edit these variables here |
| --- | --- |
| **decolace** | [`decolace_LPP/decolace_collect.py`](../decolace_collect.py) settings, around `ronchiTargetPhaseA` / `ronchiTargetPhaseB` / `ronchiCorrectKs` |
| **PACEtomo tilt series** | [`PACEtomo.py`](../../PACEtomo.py) settings, section `########## Ronchigram / laser alignment ##########` |
| **selectTargets** (if it takes ronchi images) | [`PACEtomo_selectTargets.py`](../../PACEtomo_selectTargets.py) settings, same names |

Keep decolace and PACEtomo **in sync**. Copy the same three values (plus `ronchiCorrMatrix`, `ronchiPixelSize`, `ronchiC3Offset` if you changed those) into both files.

Do **not** treat `PACEtomo_ronchi.py` module-level defaults as the place you maintain the calibration. They will be ignored once `decolace_collect.py` or `PACEtomo.py` calls `configure`.

## What each parameter is

| Setting | Meaning |
| --- | --- |
| `ronchiTargetPhaseA` | Target FFT phase [rad] of the **vertical** laser fringe |
| `ronchiTargetPhaseB` | Target FFT phase [rad] of the **horizontal** laser fringe |
| `ronchiCorrectKs` | 2×2 reference fringe-spacing vectors in 1/µm (`[[kx_v, ky_v], [kx_h, ky_h]]`). Used as the “aligned” ks; the script drives C3 so measured ks match this, then drives X-tilt so phases match A/B |
| `ronchiCorrMatrix` | Phase error → XLensDeflector coupling (scaled by 1e-5 inside the analyzer) |
| `ronchiPixelSize` | Unbinned pixel size [µm] at Trial/Record mag |
| `ronchiBinning` | Extra binning used only in the FFT (must match how you measured ks) |
| `ronchiC3Offset` | Temporary `SetImageDistanceOffset` added for the Trial shot only |
| `ronchiC3CorrectionFactor` | µm C3 per 1/µm mean diagonal ks error |
| `doRonchigram` | Must stay `True` for decolace |

## SerialEM Trial setup (required)

Trial low-dose area must match **Record** in position, beam, mag, and image shift. The only intended difference is a **very short Trial exposure** so ronchi dose is negligible.

If Trial is not at the Record beam position, the laser correction is applied at the wrong place and decolace will not stay aligned.

## How to get `correctKs` and target phases

Use the scripts in [`laser_helper/`](../../laser_helper/README.md). Run them **after** the CTF calibrations ([02_calibrations.md](02_calibrations.md)), with Trial matched to Record.

Order: **C2 stig, on plane, then on peak.**

### 0. C2 stig → equal x/y fringe spacing

C2 (condenser) astigmatism makes the two ronchi fringe spacings different. Equalize them before measuring `ronchiCorrectKs`.

Once per setup (or when C2 coupling changes), run [`calibrations/calibrate_C2_astig.py`](../../calibrations/calibrate_C2_astig.py) starting **roughly on plane**. It grids the condenser stigmator, takes a Trial at −20 at each point, measures `ks_x` / `ks_y` (the two FFT-peak magnitudes in 1/µm), and writes `c2_astig_calibration.json`. Paste that path into [`laser_helper/auto_c2_stig.py`](../../laser_helper/auto_c2_stig.py) → `c2_astig_calibration_file`.

Each session, still roughly on plane, run [`laser_helper/auto_c2_stig.py`](../../laser_helper/auto_c2_stig.py). It measures one ronchigram, then moves C2 so `ks_x ≈ ks_y` (mean spacing unchanged). C3 is restored; **C2 is left at the corrected value**.

If `|ks_x − ks_y|` is already below `tolerance` (default 0.05 1/µm), it does not move C2.

### 1. On plane → `ronchiCorrectKs`

Run [`laser_helper/auto_on_plane.py`](../../laser_helper/auto_on_plane.py) starting **roughly on plane**.

It sweeps C3 (`ImageDistanceOffset`) from +20 to +50 and −20 to −50, measures fringe-spacing magnitude, signs the negative branch, and fits the C3 where spacing is 0. Then it goes to that C3 **−20**, prints `ronchiCorrectKs`, and **leaves `ImageDistanceOffset` at the in-plane C3**.

Copy from the log into `decolace_collect.py` and `PACEtomo.py`:

- `ronchiCorrectKs = [[...], [...]]`
- `ronchiC3Offset = -20`

### 2. One-fringe X-tilt (manual)

`auto_on_peak.py` needs the X-tilt change that moves the ronchigram by **one fringe** along X and along Y. You type those into `fringe_xtilt_x` and `fringe_xtilt_y`.

With C3 on plane and a Trial at −20: change only X-tilt X until the pattern repeats (one period); that delta is `fringe_xtilt_x`. Restore X, repeat for Y. One fringe is also a 2π wrap of that fringe’s FFT phase.

Full steps: [`laser_helper/README.md`](../../laser_helper/README.md).

### 3. On peak → `ronchiTargetPhaseA` / `B`

Run [`laser_helper/auto_on_peak.py`](../../laser_helper/auto_on_peak.py) from the **working** (laser) X-tilt.

It measures defocus at `ctfXtilt` (off the laser), back-projects to the starting X-tilt if you set `xtilt_calibration_file`, then does a 5×5 X-tilt grid spanning one fringe. At each point CtfFind is run with min=max=that defocus so only phase is fitted. A quadratic fit gives the X-tilt of **maximum** phase. A −20 Trial there prints the target phases.

Copy from the log into the same two files:

- `ronchiTargetPhaseA = ...`
- `ronchiTargetPhaseB = ...`

Set the working `XLensDeflector` to the printed peak X-tilt.

After that, residual phase error in the collect log should be small, and C3 corrections should only fire when ks drift.

If `||ks error||` stays huge or phases wrap randomly, `ronchiPixelSize` / `ronchiBinning` / `ronchiPeakRadius` are wrong, or Trial is not at Record.

## Related CTF X-tilt

`ctfXtiltX` / `ctfXtiltY` in `decolace_collect.py` must be far enough from the **ronchi-aligned** working X-tilt that the probe image is off the laser. Use the same pair in the calibration scripts ([02_calibrations.md](02_calibrations.md)) and in `auto_on_peak.py`.
