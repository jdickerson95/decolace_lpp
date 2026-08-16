#!Python
# ===================================================================
# ScriptName     PACEtomo_lafis
# Purpose:       Image-shift beam-tilt / defocus / X-tilt (LAFIS) and
#                LPP frame-stack mdoc metadata. Shared by PACEtomo and
#                decolace_LPP.
# ===================================================================

import numpy as np
import serialem as sem

_sem = sem
_logger = None

hasXLens = True
xt_is_matrix = [[0.000324, -0.000347], [0.001100, 0.00028125]]
df_is_matrix = [[0.041381, 0.012342], [0.041381, 0.012342]]

lafisZeroImageShiftDefocus = None
lafisZeroImageShiftXLens = None
lafisIsDone = False
lafisXtCorrectionX = 0.0
lafisXtCorrectionY = 0.0


def _echo(text):
    if _logger is not None:
        _logger(text)
    else:
        _sem.Echo(text)


def configure(
    sem_module=None,
    logger=None,
    has_x_lens=None,
    xt_matrix=None,
    df_matrix=None,
):
    global _sem, _logger, hasXLens, xt_is_matrix, df_is_matrix
    if sem_module is not None:
        _sem = sem_module
    if logger is not None:
        _logger = logger
    if has_x_lens is not None:
        hasXLens = bool(has_x_lens)
    if xt_matrix is not None:
        xt_is_matrix = xt_matrix
    if df_matrix is not None:
        df_is_matrix = df_matrix


def calc_xt_is(xt0, is_delta):
    xt1 = [0.0, 0.0]
    xt1[0] = xt0[0] + is_delta[0] * xt_is_matrix[0][0] + is_delta[1] * xt_is_matrix[1][0]
    xt1[1] = xt0[1] + is_delta[0] * xt_is_matrix[0][1] + is_delta[1] * xt_is_matrix[1][1]
    return xt1


def calc_df_is(df0, is_delta):
    return df0 + is_delta[0] * df_is_matrix[0][0] + is_delta[1] * df_is_matrix[1][1]


def add_lpp_meta_to_next_mdoc():
    for k, v in (("ImageDistanceOffset", _sem.ReportImageDistanceOffset()),):
        v_str = "%.12f" % (float(v))
        _sem.AddToNextFrameStackMdoc(k, v_str)
    for k, v in (
        ("BeamTilt", _sem.ReportBeamTilt()),
        ("ObjectiveStig", _sem.ReportObjectiveStigmator()),
        ("XTilt", _sem.ReportXLensDeflector(2)),
    ):
        v_str = "%.6f        %.6f" % (float(v[0]), float(v[1]))
        _sem.AddToNextFrameStackMdoc(k, v_str)


def saveZeroImageShiftDefocusXLens():
    global lafisZeroImageShiftDefocus, lafisZeroImageShiftXLens
    lafisZeroImageShiftDefocus = _sem.ReportDefocus()
    if hasXLens:
        lafisZeroImageShiftXLens = _sem.ReportXLensDeflector(2)
    else:
        lafisZeroImageShiftXLens = None


def doLafis(is_x, is_y):
    global lafisIsDone, lafisXtCorrectionX, lafisXtCorrectionY
    _echo(f"WARNING: ***********doing LAFIS for image shift {is_x:.3f}, {is_y:.3f}")
    saveZeroImageShiftDefocusXLens()
    _sem.AdjustBeamTiltforIS()
    df0 = lafisZeroImageShiftDefocus
    xt0 = lafisZeroImageShiftXLens
    is_delta = (is_x, is_y)
    df1 = calc_df_is(df0, is_delta)
    _sem.SetDefocus(df1)
    if hasXLens:
        xt1 = calc_xt_is(xt0, is_delta)
        _sem.SetXLensDeflector(2, xt1[0], xt1[1])
        lafisXtCorrectionX = xt1[0] - xt0[0]
        lafisXtCorrectionY = xt1[1] - xt0[1]
    lafisIsDone = True


def restoreLafis():
    global lafisIsDone, lafisXtCorrectionX, lafisXtCorrectionY
    if not lafisIsDone:
        _echo("WARNING: LAFIS not done, can not restore")
        return
    _sem.RestoreBeamTilt()
    _sem.SetDefocus(lafisZeroImageShiftDefocus)
    if hasXLens and lafisZeroImageShiftXLens is not None:
        xt_x, xt_y = lafisZeroImageShiftXLens
        _sem.SetXLensDeflector(2, xt_x, xt_y)
        lafisXtCorrectionX = 0.0
        lafisXtCorrectionY = 0.0
    lafisIsDone = False
    _echo("WARNING: Lafis restored")
