# Navigator area and running decolace

Entry script: [`decolace_collect.py`](../decolace_collect.py).

## 1. SerialEM Record

Use the Record mag, beam, and exposure you want for the montage. Trial must match Record except for a short exposure ([03_ronchi.md](03_ronchi.md)).

`beam_radius` is the **illuminated radius in µm at Record**, used only to pack the hex lattice (FOWL spacing):

```text
step = (1 - add_overlap) * 2 * beam_radius * cos(30°)
```

Measure it from a Record image of the beam or from known illumination.

## 2. Draw the area in the Navigator

1. Take or load a View/Preview map of the region.
2. Add **3 or more Navigator points clockwise** around the area (lamella outline, holey-film patch, etc.).
3. Note their item numbers (Navigator “Item” column, not user notes).

In `decolace_collect.py`:

```python
nav_item_range = [10, 17]   # inclusive 10,11,...,17
# or an explicit list of 3+ items:
# nav_item_range = [10, 12, 15, 18, 21]
```

These points are the **polygon only**. They are not SS = (0, 0).

Optional:

```python
realign_nav_item = 28   # map/center item for RealignToItem (stage). 0 = skip
```

That item is only for getting the stage back; it is not the specimen origin.

## 3. Point the script at calibrations and geometry

See [02_calibrations.md](02_calibrations.md) and [04_geometry.md](04_geometry.md). Minimum for a serious run:

- `defocus_error_file`
- `xtilt_calibration_file`
- `astig_calibration_file` if `correctAstig = True`
- matching `ctfXtiltX` / `ctfXtiltY`
- updated ronchi `ronchiCorrectKs` / `ronchiTargetPhaseA` / `ronchiTargetPhaseB` in **this same file** ([03_ronchi.md](03_ronchi.md))
- LAFIS `xt_is_matrix` / `df_is_matrix` ([06_lafis.md](06_lafis.md))

Other collect settings:

| Setting | Typical | Meaning |
| --- | --- | --- |
| `directory` | `r"Z:\...\lamella1"` | Where MRC, frames policy, state JSON, and lattice plot go. Empty = SerialEM working directory |
| `stem` | `"decolace"` | File prefix |
| `target_defocus` | `-0.02` | Science CTF target [µm] (too close to focus to fit) |
| `ctf_probe_underfocus` | `1.0` | Physical extra underfocus of `_defocus_img` |
| `correctAstig` | `True` | Use probe astig + stig JSON every `astigEveryN` holes |
| `astigEveryN` | `10` | |
| `dry_run` | `True` first | Lattice + plot + state only |
| `test_mode` | `False` | If True, science target is `test_target_defocus` and both shots are CtfFind’d for stats |
| `test_target_defocus` | `-0.5` | Fittable science CTF target used only in test mode (production `target_defocus` stays −0.02) |
| `test_science_ctf_xtilt` | `True` | Test-mode science Record at `ctfXtilt` (off laser). `False` = working X-tilt |

## 4. Dry run, then collect

1. Set `dry_run = True`, run `decolace_collect.py`.
2. Inspect `{stem}_lattice.png` (circles = holes, red + = centroid origin) and `{stem}_state.json`.
3. If the lattice is too sparse/dense, change `beam_radius` / `add_overlap` and delete the state file (or it will resume the old lattice).
4. Set `dry_run = False`, run again.

## 5. What happens at each hole

1. Image shift to the hole.
2. Apply predicted defocus (`focus_base + delta_z + target_defocus`), using the defocus-error slope.
3. LAFIS (beam tilt / extra defocus / X-tilt vs IS).
4. Ronchi (Trial, C3 and laser X-tilt).
5. Science Record → `{stem}_{idx:04d}.mrc` (and frames with that basename).
6. Probe: +1 µm underfocus, CTF X-tilt, Record → `{stem}_{idx:04d}_defocus_img.mrc`, CtfFind, back-project defocus and astig, restore X-tilt and defocus.
7. Update spline. If six hex neighbors already have CTF, the **next** prediction for a center hole uses only those six.
8. Every `astigEveryN` completed holes: `SetObjectiveStigmator` to cancel working-X-tilt astig.

Visit order: polygon corners / poorly connected sites first (spline), then interior fill-ins once six neighbors exist.

## 6. Output

| File | What |
| --- | --- |
| `{stem}_{idx:04d}.mrc` | Science (target defocus, ronchi X-tilt) |
| `{stem}_{idx:04d}_defocus_img.mrc` | CTF probe (+1 µm, CTF X-tilt) |
| `{stem}_state.json` | Lattice, which holes are done, measured z |
| `{stem}_lattice.png` | Plot |
| camera frames | science basename; probe basename includes `_defocus_img`. Probe frames try `SetFolderForFrames` → `defocus_img` then restore |

## 7. Test mode

`test_mode` is independent of `dry_run`. Production `target_defocus` is not changed.

1. Use a distinct `stem` (e.g. `"decolace_test"`) so state and MRC files do not mix with a real collect.
2. Set `test_mode = True`, `dry_run = False`. Default science target is `test_target_defocus = -0.5` so both shots can CtfFind.
3. `test_science_ctf_xtilt = True` (default): after ronchi at the working X-tilt, the science Record is taken at `ctfXtilt` (off the laser), CtfFind, then back-projected to working X-tilt. The probe stays at CTF X-tilt. Set this `False` to take science at the normal working X-tilt (CtfFind as-taken, no back-project).
4. After each hole, and again at the end, the script writes:
   - `{stem}_test_measurements.csv` — one row per hole (defocus error, astig, image shift)
   - `{stem}_test_summary.csv` — all / center / others × science / probe
   - `{stem}_test_stats.png` — lattice colored by science defocus error; `|IS|` vs `|defocus error|` and `|astig|`; summary table
5. The SerialEM log prints the summary table at the end. Center is the hole nearest SS (0, 0). Failed CtfFind stays NaN and is counted as `fail`.

Switching back to production: set `test_mode = False` and a production `stem`. Do not resume a test-mode state file as a real collect.

## 8. Resume and restart

State is rewritten after every hole. Re-run with the same `directory` and `stem` to continue.

To start over: move or delete `{stem}_state.json` (and optionally the MRC files). Changing `nav_item_range` or `beam_radius` without deleting state will keep the old lattice.
