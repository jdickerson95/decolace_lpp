# How to run SerialEM scripts in this repo

Scripts are ordinary SerialEM Python files. They are **not** imported as a library: you open the `.py` file in a SerialEM script slot (or “Run Python script”) and edit the `############ SETTINGS ############` block at the top **before** running.

## Put the repo on the microscope PC

1. Copy the whole `PACEtomo_LPP` folder (not only `decolace_LPP/`).
2. SerialEM’s Python path must include the **repo root** (the folder that contains `PACEtomo.py`, `PACEtomo_ronchi.py`, `PACEtomo_ctf_calibrations.py`, and `decolace_LPP/`).
3. Typical extra packages (same conda env as PACEtomo): `torch`, `torch-cubic-spline-grids`. `shapely` is optional.

If imports fail with `No module named PACEtomo_ronchi`, the repo root is not on `sys.path`. Either set SerialEM’s Python path, or start SerialEM with the working directory set to the repo root.

## Settings block

Every program has a settings section **above** `import serialem`. Change values there, save the file, then run it.

Windows paths must be raw strings:

```python
defocus_error_file = r"Z:\calibrations\defocus_error_calibration.json"
```

Empty `r""` means “skip this file / use SerialEM’s current directory”.

## Order of programs

Do these in order (details in the other docs):

1. SerialEM low-dose (Record + Trial matched) — [03_ronchi.md](03_ronchi.md)
2. Calibrations — [02_calibrations.md](02_calibrations.md)
3. Laser on-plane then on-peak (`laser_helper/`) — [03_ronchi.md](03_ronchi.md)
4. Optional geometry map — [04_geometry.md](04_geometry.md)
5. Navigator polygon + decolace collect — [05_area_and_collect.md](05_area_and_collect.md)

LAFIS matrices: [06_lafis.md](06_lafis.md)
