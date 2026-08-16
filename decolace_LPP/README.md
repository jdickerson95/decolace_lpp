# decolace_LPP

Hexagonal 2D montage on SerialEM / LPP. Each hole: predicted defocus, ronchi (laser), science Record near focus, then a +1 µm `_defocus_img` for CTF. Collection order is “six hex neighbors, then fill the center.”

This README is the map. Step-by-step for each program is under [`docs/`](docs/).

## Read these in order

1. [docs/01_serialem_setup.md](docs/01_serialem_setup.md) — copy the repo, Python path, how settings blocks work
2. [docs/02_calibrations.md](docs/02_calibrations.md) — **what each calibration writes and which setting gets the JSON path**
3. [docs/03_ronchi.md](docs/03_ronchi.md) — **you must set `ronchiCorrectKs` and `ronchiTargetPhaseA` / `B` in `decolace_collect.py` (and PACEtomo.py)**
4. [docs/04_geometry.md](docs/04_geometry.md) — measureGeometry; decolace origin is the hex centroid, not the first geo point
5. [docs/05_area_and_collect.md](docs/05_area_and_collect.md) — navigator polygon, dry run, collect, resume, **test mode**
6. [docs/06_lafis.md](docs/06_lafis.md) — `xt_is_matrix` / `df_is_matrix`

## What to do with calibration results (short)

Each calibration script writes a JSON (plus CSV/PNG). **Paste the full JSON path into `decolace_LPP/decolace_collect.py`:**

| JSON from | Setting in `decolace_collect.py` |
| --- | --- |
| `calibrations/calibrate_defocus_error.py` → `defocus_error_calibration.json` | `defocus_error_file = r"..."` |
| `calibrations/check_xtilt_defoc_astig.py` → `xtilt_defoc_astig_calibration.json` | `xtilt_calibration_file = r"..."` |
| `calibrations/calibrate_astigmatism.py` → `astigmatism_calibration.json` | `astig_calibration_file = r"..."` |
| `PACEtomo_measureGeometry.py` → `geometry.json` | `geometry_file = r"..."` |

Empty `r""` skips that file. `correctAstig = True` needs the xtilt JSON **and** the astig JSON.

Use the same `ctfXtiltX` / `ctfXtiltY` in decolace and in every calibration that takes a CTF off the laser.

Details: [docs/02_calibrations.md](docs/02_calibrations.md).

## Ronchi parameters you must update

Do **not** leave the checked-in `ronchiCorrectKs` / `ronchiTargetPhaseA` / `ronchiTargetPhaseB`. They are microscope-session specific.

**Edit them in [`decolace_collect.py`](decolace_collect.py)** (settings block, `ronchiTargetPhaseA`, `ronchiTargetPhaseB`, `ronchiCorrectKs`).

If you also run tilt series, copy the **same three values** into [`PACEtomo.py`](../PACEtomo.py) under `########## Ronchigram / laser alignment ##########`.

`PACEtomo_ronchi.py` is only fallback defaults; the collect/PACEtomo settings overwrite them at run time.

How to measure new values: run [`laser_helper/auto_on_plane.py`](../laser_helper/auto_on_plane.py) then [`laser_helper/auto_on_peak.py`](../laser_helper/auto_on_peak.py). Steps: [docs/03_ronchi.md](docs/03_ronchi.md) and [`laser_helper/README.md`](../laser_helper/README.md).

## Code map

| File | Role |
| --- | --- |
| [`decolace_collect.py`](decolace_collect.py) | SerialEM entry: **all user settings live here** |
| [`area.py`](area.py) | Clockwise nav polygon → hex lattice; SS origin = mean of holes |
| [`schedule.py`](schedule.py) | 6-neighbor graph; seed then fill |
| [`predict.py`](predict.py) | Spline vs 6-neighbor defocus; live refit |
| [`acquire.py`](acquire.py) | One hole: defocus → LAFIS → ronchi → science → probe |
| [`test_stats.py`](test_stats.py) | Test-mode CTF stats (center vs others, image shift vs defocus/astig) |

Shared with PACEtomo (do not put session calibrations only in these): `PACEtomo_ronchi.py`, `PACEtomo_lafis.py`, `PACEtomo_ctf_calibrations.py`, `PACEtomo_geometry.py`.

## Per-hole sequence (once `dry_run = False`)

1. Image shift
2. Predicted defocus (`ChangeFocus` scaled by defocus-error JSON)
3. LAFIS
4. Ronchi (Trial)
5. Science `{stem}_{idx:04d}.mrc` at `target_defocus` (default −0.02 µm)
6. Probe `{stem}_{idx:04d}_defocus_img.mrc` at +1 µm underfocus and CTF X-tilt; CtfFind; back-project
7. Update spline; every `astigEveryN` holes, cancel astig from the probe

State: `{stem}_state.json` after every hole. Same `directory` + `stem` resumes. Delete the state file to pack a new lattice.

**Test mode:** set `test_mode = True` and a distinct `stem` (e.g. `decolace_test`). Science uses `test_target_defocus` (−0.5 µm) so both shots CtfFind. Default `test_science_ctf_xtilt = True` takes science at `ctfXtilt`. Writes `{stem}_test_measurements.csv`, `{stem}_test_summary.csv`, `{stem}_test_stats.png`. Details: [docs/05_area_and_collect.md](docs/05_area_and_collect.md).
