# LAFIS matrices (image-shift X-tilt and defocus)

decolace applies PACEtomo’s LAFIS correction on every hole after the defocus guess and before ronchi:

- `xt_is_matrix` — extra XLensDeflector(2) vs image shift
- `df_is_matrix` — extra defocus vs image shift
- SerialEM `AdjustBeamTiltforIS()` (coma vs IS) when `beamTiltComp = True`

## Where to edit

| Script | Variables |
| --- | --- |
| [`decolace_LPP/decolace_collect.py`](../decolace_collect.py) | `xt_is_matrix`, `df_is_matrix`, `beamTiltComp` |
| [`PACEtomo.py`](../../PACEtomo.py) | same names in the `######### LAFIS` settings section |

Copy the same matrices into both files. `PACEtomo_lafis.py` only holds fallbacks; the running script’s settings win.

## How to get the numbers

Use [`LAFIS_manual_calibration.py`](../../LAFIS_manual_calibration.py) (and your usual IS vs X-tilt / defocus notes):

1. At Record, IS = 0, note X-tilt and defocus.
2. Apply a known IS (e.g. 5 µm on one axis), realign laser/focus as you would for data, note the new X-tilt and defocus.
3. Fill `xt_is_matrix` / `df_is_matrix` so `xt1 = xt0 + IS · matrix` matches what you measured.

The matrices in the repo are from one session (`#26jul23` in PACEtomo). Replace them for this microscope.

`beamTiltComp = True` also requires SerialEM’s coma-vs-image-shift calibration (`ReportComaVsISmatrix` must succeed).
