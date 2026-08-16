#!Python
# ===================================================================
# ScriptName     decolace_LPP.area
# Purpose:       Navigator polygon -> hexagonal hole lattice, SS origin
#                at the mean of all holes (not a geo point).
# ===================================================================

import os
from datetime import datetime

import numpy as np

try:
    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon
    _HAS_SHAPELY = True
except ImportError:
    _HAS_SHAPELY = False
    Point = None
    Polygon = None


def nav_indices_from_range(nav_item_range):
    """Inclusive [start, end] or an explicit list of navigator indices."""
    vals = [int(v) for v in nav_item_range]
    if len(vals) == 2 and vals[1] >= vals[0] and vals[1] - vals[0] < 500:
        # Treat as inclusive range when it looks like [start, end].
        # An explicit list of two items is still a valid 2-vertex degenerate
        # case; require the user to pass 3+ for a polygon, so 2 means range.
        return list(range(vals[0], vals[1] + 1))
    return vals


def _report_nav_stage(sem, nav_idx):
    index, x, y, z, *_ = sem.ReportOtherItem(int(nav_idx))
    return int(index), np.array([float(x), float(y), float(z)], dtype=float)


def specimen_from_stage_diff(sem, dx, dy):
    sem.SetImageShift(0, 0)
    sem.ImageShiftByStageDiff(float(dx), float(dy))
    ss = sem.ReportSpecimenShift()
    sem.SetImageShift(0, 0)
    return np.array([float(ss[0]), float(ss[1])], dtype=float)


def read_polygon_vertices(sem, nav_ids, log=None):
    """Stage and specimen coords of clockwise navigator points relative to first item."""
    echo = log or (lambda *_a, **_k: None)
    if len(nav_ids) < 3:
        raise ValueError("Need at least 3 navigator points for a polygon.")
    vertices = []
    ref_idx, ref_stage = _report_nav_stage(sem, nav_ids[0])
    echo(f"Polygon reference nav item {ref_idx}: stage ({ref_stage[0]:.3f}, {ref_stage[1]:.3f})")
    for nav_idx in nav_ids:
        idx, stage = _report_nav_stage(sem, nav_idx)
        dxy = stage[:2] - ref_stage[:2]
        ss = specimen_from_stage_diff(sem, dxy[0], dxy[1])
        vertices.append(
            {
                "nav_item": int(idx),
                "stage": [float(v) for v in stage],
                "ss_ref": [float(ss[0]), float(ss[1])],
            }
        )
        echo(
            f"  nav {idx}: stage ({stage[0]:.3f}, {stage[1]:.3f}) "
            f"SS_ref ({ss[0]:.3f}, {ss[1]:.3f})"
        )
    return {
        "nav_ids": [int(v) for v in nav_ids],
        "reference_nav_item": int(ref_idx),
        "reference_stage": [float(v) for v in ref_stage],
        "vertices": vertices,
    }


def hex_step(beam_radius, add_overlap=0.05):
    return (1.0 - float(add_overlap)) * 2.0 * float(beam_radius) * np.cos(np.deg2rad(30.0))


def _point_in_polygon(verts, point):
    """Ray-cast even-odd rule (numpy only)."""
    x, y = float(point[0]), float(point[1])
    verts = np.asarray(verts, dtype=float)
    x0, y0 = verts[:, 0], verts[:, 1]
    x1, y1 = np.roll(x0, -1), np.roll(y0, -1)
    cond = ((y0 > y) != (y1 > y)) & (x < (x1 - x0) * (y - y0) / (y1 - y0 + 1e-30) + x0)
    return bool(np.count_nonzero(cond) % 2)


def _polygon_contains(verts, point, buffer_um=0.001):
    """True if point is inside the polygon (shapely if available)."""
    xy = (float(point[0]), float(point[1]))
    if _HAS_SHAPELY:
        poly = Polygon([tuple(v) for v in verts]).buffer(buffer_um)
        return bool(poly.contains(Point(xy)))
    return _point_in_polygon(verts, xy)


def pack_hex_lattice(ss_vertices, beam_radius, add_overlap=0.05, expansion=1.0, direction=1, rotation_direction=1):
    """Hex packing inside the polygon using FOWL spacing (bbox grid, polygon clip)."""
    verts = np.asarray(ss_vertices, dtype=float) * float(expansion)
    if verts.shape[0] < 3:
        raise ValueError("Need at least 3 specimen vertices.")
    if _HAS_SHAPELY:
        polygon = Polygon([tuple(v) for v in verts]).buffer(0.001)
        if polygon.is_empty:
            raise ValueError("Polygon is empty; check clockwise navigator points.")
        contains = lambda p: polygon.contains(Point(tuple(p)))
    else:
        contains = lambda p: _polygon_contains(verts, p)

    step = hex_step(beam_radius, add_overlap)
    dy = step * np.sqrt(3.0) / 2.0
    xmin, ymin = verts.min(axis=0) - 0.5 * step
    xmax, ymax = verts.max(axis=0) + 0.5 * step
    positions = []
    row = 0
    y = ymin
    while y <= ymax + 1e-9:
        x = xmin + (0.5 * step if row % 2 else 0.0)
        while x <= xmax + 1e-9:
            p = np.array([x, y], dtype=float)
            if contains(p):
                positions.append(p)
            x += step
        y += dy
        row += 1
    if not positions:
        raise ValueError("Hex packing produced no holes; check beam_radius and navigator points.")
    return np.array(positions, dtype=float)


def center_on_centroid(positions):
    positions = np.asarray(positions, dtype=float)
    centroid = np.mean(positions, axis=0)
    return positions - centroid, centroid


def build_area(
    sem,
    nav_item_range,
    beam_radius,
    add_overlap=0.05,
    expansion=1.0,
    log=None,
):
    """Read nav polygon, pack hex lattice, shift SS origin to mean of all holes."""
    nav_ids = nav_indices_from_range(nav_item_range)
    poly = read_polygon_vertices(sem, nav_ids, log=log)
    ss_verts = np.array([v["ss_ref"] for v in poly["vertices"]], dtype=float)
    packed = pack_hex_lattice(ss_verts, beam_radius, add_overlap=add_overlap, expansion=expansion)
    centered, centroid_ss_ref = center_on_centroid(packed)
    verts_centered = ss_verts - centroid_ss_ref
    holes = []
    for i, ss in enumerate(centered):
        holes.append(
            {
                "index": i,
                "ssx": float(ss[0]),
                "ssy": float(ss[1]),
                "acquired": False,
                "z_focus": None,
                "defocus_um": None,
                "astig_x_um": None,
                "astig_y_um": None,
                "predict_mode": None,
            }
        )
    echo = log or (lambda *_a, **_k: None)
    echo(f"Hex lattice: {len(holes)} holes, centroid SS_ref ({centroid_ss_ref[0]:.3f}, {centroid_ss_ref[1]:.3f})")
    echo("SS origin is the mean of all holes, not a navigator or geo point.")
    return {
        "created": datetime.now().isoformat(timespec="seconds"),
        "beam_radius": float(beam_radius),
        "add_overlap": float(add_overlap),
        "hex_step_um": float(hex_step(beam_radius, add_overlap)),
        "polygon": poly,
        "centroid_ss_ref": [float(centroid_ss_ref[0]), float(centroid_ss_ref[1])],
        "vertices_ss": verts_centered.tolist(),
        "holes": holes,
    }


def plot_area(area, path=None, show=True):
    import matplotlib.pyplot as plt

    verts = np.array(area["vertices_ss"], dtype=float)
    holes = area["holes"]
    radius = float(area["beam_radius"])
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.add_patch(plt.Polygon(verts, closed=True, fill=None, edgecolor="r", linewidth=1.5))
    for hole in holes:
        color = "0.5" if hole.get("acquired") else "b"
        ax.add_patch(
            plt.Circle((hole["ssx"], hole["ssy"]), radius=radius, fill=None, edgecolor=color, linewidth=0.6)
        )
    ax.plot(0, 0, "r+", markersize=12, label="centroid (SS origin)")
    ax.set_aspect("equal")
    ax.set_xlabel("SSX (um)")
    ax.set_ylabel("SSY (um)")
    ax.set_title(f"decolace_LPP hex lattice ({len(holes)} holes)")
    ax.legend()
    fig.tight_layout()
    if path:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        fig.savefig(path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return path
