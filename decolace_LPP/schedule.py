#!Python
# ===================================================================
# ScriptName     decolace_LPP.schedule
# Purpose:       Hex-neighbor graph and 6-then-fill visit order.
# ===================================================================

import numpy as np


def hole_xy(holes):
    return np.array([[h["ssx"], h["ssy"]] for h in holes], dtype=float)


def neighbor_graph(holes, step_um, tol=0.15):
    """Undirected neighbors: distance within step * (1 +/- tol), excluding self."""
    xy = hole_xy(holes)
    step = float(step_um)
    lo = step * (1.0 - float(tol))
    hi = step * (1.0 + float(tol))
    n = len(holes)
    neighbors = [[] for _ in range(n)]
    for i in range(n):
        d = np.linalg.norm(xy - xy[i], axis=1)
        for j in np.where((d >= lo) & (d <= hi))[0]:
            j = int(j)
            if j != i:
                neighbors[i].append(j)
    return neighbors


def measured_neighbor_indices(holes, neighbors, idx):
    return [j for j in neighbors[idx] if holes[j].get("z_focus") is not None]


def is_fill_ready(holes, neighbors, idx):
    if holes[idx].get("acquired"):
        return False
    return len(measured_neighbor_indices(holes, neighbors, idx)) >= 6


def next_hole(holes, neighbors, vertex_ss):
    """Prefer a 6-neighbor fill-in; else the best spline seed (corners, then wavefront)."""
    pending = [i for i, h in enumerate(holes) if not h.get("acquired")]
    if not pending:
        return None, None
    verts = np.asarray(vertex_ss, dtype=float)
    fill = [i for i in pending if is_fill_ready(holes, neighbors, i)]
    if fill:
        def fill_key(i):
            nmeas = len(measured_neighbor_indices(holes, neighbors, i))
            dist0 = float(np.hypot(holes[i]["ssx"], holes[i]["ssy"]))
            return (-nmeas, dist0, i)
        idx = sorted(fill, key=fill_key)[0]
        return idx, "fill6"
    def seed_key(i):
        nmeas = len(measured_neighbor_indices(holes, neighbors, i))
        xy = np.array([holes[i]["ssx"], holes[i]["ssy"]], dtype=float)
        dvert = float(np.min(np.linalg.norm(verts - xy, axis=1))) if len(verts) else 0.0
        return (-nmeas, dvert, i)
    idx = sorted(pending, key=seed_key)[0]
    return idx, "spline"
