# Measure geometry for decolace

Script: [`PACEtomo_measureGeometry.py`](../../PACEtomo_measureGeometry.py) at the repo root.

This is **optional**. decolace can start from a flat map (`geometry_file = r""`) and build a spline from the `_defocus_img` CTFs as it collects. A pre-map helps the first holes, which otherwise only have a spline (no 6 neighbors yet).

## PACEtomo “first = origin” vs decolace

In `PACEtomo_measureGeometry.py`, `geo_nav_item_list[0]` is SS = (0, 0) because PACEtomo targets are stored relative to a tracking item.

**decolace ignores that as the montage origin.** For decolace:

- Navigator polygon points only define the **area**.
- Specimen `(0, 0)` is the **mean of all hex holes**.
- Geo points can be anywhere on the lamella. None of them has to be the center, and the center does not have to be measured.
- On load, decolace re-expresses the JSON in the hex-centroid frame and shifts z so height at the centroid is 0 (interpolated).

You can still run measureGeometry exactly as for PACEtomo. Just do not assume the first geo nav item is the decolace center.

## How to run

1. Add ≥3 navigator points on the lamella (spread them; more is better for a 10×10 live spline later).
2. Set `geo_nav_item_list` to those item numbers.
3. Paste the same calibration JSON paths as decolace:
   - `xtilt_calibration_file`
   - `defocus_error_file`
4. Match `ctfXtiltX` / `ctfXtiltY`.
5. For a pre-map that decolace will refine, `spline_resolution = (3, 3)` in measureGeometry is enough; decolace refits up to `(10, 10)` during collection.
6. Set `save_dir` if you want a stable folder. Output JSON default name: `geometry.json`.
7. In `decolace_collect.py` set:

```python
geometry_file = r"Z:\...\geometry.json"
spline_resolution = (10, 10)
```

## What decolace does with it

- Loads the spline/plane.
- Translates SS using the measureGeometry origin nav item vs the polygon reference, then subtracts the hex centroid.
- Subtracts z at (0, 0) in that frame.
- After each good probe CTF, adds that hole and refits (coarse grid until enough points, then 10×10).
- Once a hole has **six measured hex neighbors**, prediction uses **only those six** (local plane), not the spline.
