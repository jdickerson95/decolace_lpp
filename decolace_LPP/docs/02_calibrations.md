# Calibration scripts: what they write and where to paste it

Run these from SerialEM like any other script. Each writes a **JSON** (the calibration), plus a CSV and a PNG plot. The JSON is what decolace and PACEtomo actually read.

**You must copy the JSON path into the consumer script.** Nothing is applied automatically.

## Where to paste JSON paths

| You ran | JSON it writes (default name) | Paste the full path into |
| --- | --- | --- |
| `calibrations/calibrate_defocus_error.py` | `defocus_error_calibration.json` | `decolace_collect.py` → `defocus_error_file`<br>`PACEtomo.py` → `defocus_error_file`<br>`PACEtomo_measureGeometry.py` → `defocus_error_file` |
| `calibrations/check_xtilt_defoc_astig.py` | `xtilt_defoc_astig_calibration.json` | `decolace_collect.py` → `xtilt_calibration_file`<br>`PACEtomo.py` → `xtilt_calibration_file`<br>`calibrate_astigmatism.py` → `xtilt_calibration_file`<br>`PACEtomo_measureGeometry.py` → `xtilt_calibration_file` |
| `calibrations/calibrate_astigmatism.py` | `astigmatism_calibration.json` | `decolace_collect.py` → `astig_calibration_file`<br>`calibrations/correct_astigmatism.py` → `astig_calibration_file` |
| `PACEtomo_measureGeometry.py` | `geometry.json` | `decolace_collect.py` → `geometry_file`<br>`PACEtomo.py` → `geometry_file` (and `geometryMode = "spline"`) |

Example in `decolace_LPP/decolace_collect.py`:

```python
xtilt_calibration_file = r"Z:\26aug15\cals\xtilt_defoc_astig_calibration.json"
defocus_error_file     = r"Z:\26aug15\cals\defocus_error_calibration.json"
astig_calibration_file = r"Z:\26aug15\cals\astigmatism_calibration.json"
geometry_file          = r"Z:\26aug15\lamella1\geometry.json"
```

If a path is `r""`, that correction is skipped. `correctAstig = True` needs **both** the xtilt JSON and the astig JSON.

Also copy **`ctfXtiltX` / `ctfXtiltY`** (and `useCtfXtilt`) so they are the **same numbers** in every script that uses CTF X-tilt (PACEtomo, decolace, measureGeometry, calibrate_defocus_error, calibrate_astigmatism).

## Output directory

Each calibration script has:

```python
save_dir = r""   # empty = SerialEM current directory
```

Set `save_dir` to a calibrations folder you will keep, e.g. `r"Z:\cals"`. The log prints the absolute JSON path when the run finishes. Use that path.

## 1. `calibrate_defocus_error.py`

**What it does:** commanded `ChangeFocus` steps vs CtfFind. Fits `measured = intercept + slope * commanded`.

**Why decolace needs it:** science defocus and the +1 µm probe are applied with `ChangeFocus(desired / slope)`. Without the JSON, slope is treated as 1.

**How to run:**

1. Record mag, carbon or a region with a reliable CTF (~−1 to −4 µm).
2. Match `use_ctf_xtilt` / `ctf_xtilt_x` / `ctf_xtilt_y` to PACEtomo/decolace.
3. Set `save_dir`, run.
4. Check `defocus_error_fit.png`. Slope should be a finite number near 1 (not 0).
5. Paste `defocus_error_calibration.json` into `defocus_error_file` as above.

## 2. `check_xtilt_defoc_astig.py`

**What it does:** X-tilt grid around `start_xtilt_x/y`, CtfFind at each point, plane fits for defocus and 2θ astig vs X-tilt.

**Why decolace needs it:** the `_defocus_img` is taken at `ctfXtilt` (off the laser). Defocus and astig are mapped back to the ronchi/working X-tilt with this plane.

**How to run:**

1. Same mag/CTF as collection. Set `start_xtilt_x/y` to the **working** (laser-aligned) X-tilt you use at Record, or 0 if that is your origin.
2. `grid_size` 5 and `xtilt_step` 0.0005 is a typical start; the grid should span as far as `ctfXtilt` is from working X-tilt.
3. Set `save_dir`, run.
4. Check `xtilt_defoc_astig_fits.png`.
5. Paste `xtilt_defoc_astig_calibration.json` into `xtilt_calibration_file`.

## 3. `calibrate_astigmatism.py`

**What it does:** objective-stigmator grid, CtfFind astig, fit `[astig_x, astig_y] = M @ dStig + b`. Cancel with `dStig = -inv(M) @ astig`.

**Why decolace needs it:** every `astigEveryN` holes, shot 2’s astig (back-projected to working X-tilt) is cancelled with the stigmator. Goal is ~zero astig.

**How to run:**

1. Run **after** the xtilt JSON exists if `useCtfXtilt = True`. Set `xtilt_calibration_file` in **this** script too.
2. Match `ctfXtiltX/Y` to decolace.
3. Set `save_dir`, run. Stigmator is restored at the end.
4. Check `astigmatism_fit.png` / that `inv_M` exists in the JSON.
5. Paste `astigmatism_calibration.json` into `decolace_collect.py` → `astig_calibration_file`.

Optional one-shot test: `calibrations/correct_astigmatism.py` with the same two JSON paths.

## 4. Beam-tilt scripts (optional for decolace)

These feed **PACEtomo autofocus**, not the decolace two-shot CTF loop. Run them if you also collect tilt series.

| Script | Typical result you copy |
| --- | --- |
| `calibrate_beam_tilt_scaling.py` | `beam_tilt_correction` / `defocus_tilt_correction` numbers into `PACEtomo.py` |
| `calibrate_beam_tilt_correction.py` | same family of scale factors |
| `calibrate_beam_tilt_serialEM.py` | SerialEM beam-tilt vs defocus check |
| `calibrate_beam_tilt_xtilt_matrix.py` | beam-tilt / X-tilt coupling JSON; also used when diagnosing autofocus vs X-tilt |

decolace image shift uses **LAFIS** matrices instead; see [06_lafis.md](06_lafis.md).

## Keep a copy

Do not leave JSON only in SerialEM’s current directory if that folder changes per session. A stable `Z:\cals\` (or similar) path in every settings block is the intended workflow.
