#!Python
# ===================================================================
# ScriptName     decolace_LPP.predict
# Purpose:       Defocus prediction from spline (seeds) or 6 hex
#                neighbors (fill-in). Live spline refit. Geometry is
#                expressed with SS origin at the hex-lattice centroid.
# ===================================================================

import numpy as np

import PACEtomo_geometry as pacegeo


def _finite_pairs(ssx, ssy, z):
    ssx = np.asarray(ssx, dtype=float).ravel()
    ssy = np.asarray(ssy, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    mask = np.isfinite(ssx) & np.isfinite(ssy) & np.isfinite(z)
    return ssx[mask], ssy[mask], z[mask]


def recenter_geometry(geometry, centroid_ss_ref, origin_ss_ref=None):
    """Move measureGeometry points into the lattice-centroid SS frame.

    origin_ss_ref: specimen coords of PACEtomo's first geo item, in the same
    frame as centroid_ss_ref (polygon-reference SS). If None, geo SS is
    assumed already in that reference frame.
    z is shifted so evaluate_z(0, 0) == 0 at the montage centroid.
    """
    ssx = np.array(geometry["points"]["ssx"], dtype=float)
    ssy = np.array(geometry["points"]["ssy"], dtype=float)
    z = np.array(geometry["points"]["z_um"], dtype=float)
    centroid = np.asarray(centroid_ss_ref, dtype=float)
    if origin_ss_ref is not None:
        origin = np.asarray(origin_ss_ref, dtype=float)
        ssx = ssx + origin[0] - centroid[0]
        ssy = ssy + origin[1] - centroid[1]
    else:
        ssx = ssx - centroid[0]
        ssy = ssy - centroid[1]
    res = tuple(geometry.get("spline", {}).get("resolution") or (3, 3))
    n_iter = int(geometry.get("spline", {}).get("n_iter") or 800)
    recentered = pacegeo.build_geometry(ssx, ssy, z, spline_resolution=res, spline_n_iter=n_iter)
    z0 = float(np.asarray(pacegeo.evaluate_z(0.0, 0.0, recentered)))
    recentered = pacegeo.build_geometry(ssx, ssy, z - z0, spline_resolution=res, spline_n_iter=n_iter)
    recentered["z_offset_removed_um"] = z0
    recentered["ss_origin"] = "hex_lattice_centroid"
    return recentered


def origin_ss_from_nav(sem, origin_nav_item, reference_nav_item, reference_stage, log=None):
    """Specimen coords of origin_nav_item in the polygon-reference SS frame."""
    try:
        from decolace_LPP.area import specimen_from_stage_diff
    except ImportError:
        from area import specimen_from_stage_diff

    echo = log or (lambda *_a, **_k: None)
    if int(origin_nav_item) == int(reference_nav_item):
        echo("Geometry origin nav item is the polygon reference; origin SS_ref = (0, 0)")
        return np.array([0.0, 0.0], dtype=float)
    _idx, stage = _report_stage(sem, origin_nav_item)
    dxy = np.array(stage[:2], dtype=float) - np.array(reference_stage[:2], dtype=float)
    ss = specimen_from_stage_diff(sem, dxy[0], dxy[1])
    echo(f"Geometry origin nav {origin_nav_item} SS_ref ({ss[0]:.3f}, {ss[1]:.3f})")
    return ss


def _report_stage(sem, nav_idx):
    index, x, y, z, *_ = sem.ReportOtherItem(int(nav_idx))
    return int(index), np.array([float(x), float(y), float(z)], dtype=float)


def adaptive_resolution(n_points, requested=(10, 10)):
    nx, ny = int(requested[0]), int(requested[1])
    n = int(n_points)
    if n < 6:
        return (3, 3)
    if n < max(6, (nx * ny) // 4):
        return (3, 3)
    return (nx, ny)


def fit_live_geometry(ssx, ssy, z, requested_resolution=(10, 10), n_iter=800):
    ssx, ssy, z = _finite_pairs(ssx, ssy, z)
    if ssx.size < 3:
        return None
    res = adaptive_resolution(ssx.size, requested_resolution)
    return pacegeo.build_geometry(ssx, ssy, z, spline_resolution=res, spline_n_iter=n_iter)


def predict_spline(ssx, ssy, geometry):
    if geometry is None:
        return 0.0
    val = pacegeo.evaluate_z(ssx, ssy, geometry)
    return float(np.asarray(val))


def fit_local_plane(ssx, ssy, z):
    ssx, ssy, z = _finite_pairs(ssx, ssy, z)
    if ssx.size < 3:
        return None
    A = np.column_stack((ssx, ssy, np.ones(ssx.size)))
    coeff, *_ = np.linalg.lstsq(A, z, rcond=None)
    return coeff


def predict_fill6(ssx, ssy, neighbor_ss, neighbor_z):
    coeff = fit_local_plane(
        [p[0] for p in neighbor_ss],
        [p[1] for p in neighbor_ss],
        neighbor_z,
    )
    if coeff is None:
        return float(np.nanmean(neighbor_z))
    return float(coeff[0] * ssx + coeff[1] * ssy + coeff[2])


def z_focus_from_probe(
    delta_z_used,
    measured_defocus_um,
    target_defocus,
    ctf_probe_underfocus,
):
    """Specimen height residual in CTF microns, same frame as predict_delta_z.

    Science is aimed at target_defocus. The probe is commanded a further
    ctf_probe_underfocus um underfocus (ChangeFocus scaled by the defocus-error
    slope), so the expected probe CTF at working X-tilt is
        expected = target_defocus - ctf_probe_underfocus
    (e.g. -0.02 - 1.0 = -1.02 um).

    Do not mix SerialEM ReportDefocus with CtfFind. ReportDefocus after LAFIS
    SetDefocus and slope-scaled ChangeFocus is not the same number as the
    back-projected CTF, so probe_nominal - measured would walk the spline.

    If the probe CTF is more underfocus than expected (measured more negative),
    the next prediction must go more positive so science lands on target:
        z_focus = delta_z_used - (measured - expected)
    """
    expected = float(target_defocus) - float(ctf_probe_underfocus)
    error = float(measured_defocus_um) - expected
    return float(delta_z_used) - error


def predict_delta_z(holes, neighbors, idx, geometry, mode):
    hole = holes[idx]
    if mode == "fill6":
        nidx = [j for j in neighbors[idx] if holes[j].get("z_focus") is not None]
        if len(nidx) >= 6:
            xy = np.array([holes[idx]["ssx"], holes[idx]["ssy"]], dtype=float)
            nidx = sorted(
                nidx,
                key=lambda j: np.hypot(holes[j]["ssx"] - xy[0], holes[j]["ssy"] - xy[1]),
            )[:6]
            nss = [(holes[j]["ssx"], holes[j]["ssy"]) for j in nidx]
            nz = [holes[j]["z_focus"] for j in nidx]
            return predict_fill6(hole["ssx"], hole["ssy"], nss, nz), "fill6"
    return predict_spline(hole["ssx"], hole["ssy"], geometry), "spline"


def geometry_from_holes(holes, extra_ss=None, extra_z=None, requested_resolution=(10, 10), n_iter=800):
    ssx = [h["ssx"] for h in holes if h.get("z_focus") is not None]
    ssy = [h["ssy"] for h in holes if h.get("z_focus") is not None]
    z = [h["z_focus"] for h in holes if h.get("z_focus") is not None]
    if extra_ss is not None and extra_z is not None:
        ssx.extend([p[0] for p in extra_ss])
        ssy.extend([p[1] for p in extra_ss])
        z.extend(list(extra_z))
    return fit_live_geometry(ssx, ssy, z, requested_resolution=requested_resolution, n_iter=n_iter)
