#!Python
# ===================================================================
# ScriptName     PACEtomo_geometry
# Purpose:       Plane + cubic B-spline residual height map z(SSX, SSY).
#                Torch is imported only when fitting or evaluating a spline.
# ===================================================================

import json
from datetime import datetime

import numpy as np

MODEL_NAME = "plane_plus_bspline_residual"


def plane_z(ssx, ssy, pretilt, rotation):
    """PACEtomo plane height [um]; 0 at specimen origin."""
    ssx = np.asarray(ssx, dtype=float)
    ssy = np.asarray(ssy, dtype=float)
    return np.tan(np.radians(pretilt)) * (
        np.cos(np.radians(rotation)) * ssy - np.sin(np.radians(rotation)) * ssx
    )


def fit_plane(ssx, ssy, z):
    """SVD plane through (SSX, SSY, z). Returns pretilt/rotation like selectTargets."""
    ssx = np.asarray(ssx, dtype=float).ravel()
    ssy = np.asarray(ssy, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    mask = np.isfinite(ssx) & np.isfinite(ssy) & np.isfinite(z)
    ssx, ssy, z = ssx[mask], ssy[mask], z[mask]
    n = int(ssx.size)
    if n < 3:
        raise ValueError(f"Need at least 3 finite geo points to fit a plane (got {n}).")
    points = np.vstack((ssx, ssy, z))
    center = np.mean(points, axis=1)
    svd = np.linalg.svd(points - center[:, None])
    norm = svd[0][:, -1]
    if abs(float(norm[2])) < 1e-8:
        raise ValueError("Fitted plane is nearly vertical; check geo-point z values.")
    errors = []
    for point in zip(ssx, ssy, z):
        errors.append(float(np.dot(norm, np.array(point) - center) ** 2))
    mean_error = float(np.mean(errors))
    z_plane = (
        center[2]
        - (norm[0] * (ssx - center[0]) + norm[1] * (ssy - center[1])) / norm[2]
    )
    rms = float(np.sqrt(np.mean((z - z_plane) ** 2)))
    sign = 1 if norm[1] <= 0 else -1
    pretilt = float(sign * np.degrees(np.arccos(np.clip(norm[2], -1.0, 1.0))))
    rotation = float(-np.degrees(np.arctan(norm[0] / norm[1])))
    return {
        "n": n,
        "center": [float(c) for c in center],
        "norm": [float(v) for v in norm],
        "pretilt": pretilt,
        "rotation": rotation,
        "mean_error": mean_error,
        "rms_um": rms,
        "errors": errors,
    }


def _bbox_from_points(ssx, ssy, pad_frac=0.1, min_span=1.0):
    ssx = np.asarray(ssx, dtype=float)
    ssy = np.asarray(ssy, dtype=float)
    ssx_min, ssx_max = float(np.min(ssx)), float(np.max(ssx))
    ssy_min, ssy_max = float(np.min(ssy)), float(np.max(ssy))
    ssx_span = max(ssx_max - ssx_min, min_span)
    ssy_span = max(ssy_max - ssy_min, min_span)
    ssx_pad = ssx_span * pad_frac
    ssy_pad = ssy_span * pad_frac
    ssx_mid = 0.5 * (ssx_min + ssx_max)
    ssy_mid = 0.5 * (ssy_min + ssy_max)
    return {
        "ssx_min": ssx_mid - 0.5 * ssx_span - ssx_pad,
        "ssx_max": ssx_mid + 0.5 * ssx_span + ssx_pad,
        "ssy_min": ssy_mid - 0.5 * ssy_span - ssy_pad,
        "ssy_max": ssy_mid + 0.5 * ssy_span + ssy_pad,
    }


def normalize_ss(ssx, ssy, bbox):
    ssx = np.asarray(ssx, dtype=float)
    ssy = np.asarray(ssy, dtype=float)
    dx = bbox["ssx_max"] - bbox["ssx_min"]
    dy = bbox["ssy_max"] - bbox["ssy_min"]
    u = (ssx - bbox["ssx_min"]) / dx if dx != 0 else np.full_like(ssx, 0.5)
    v = (ssy - bbox["ssy_min"]) / dy if dy != 0 else np.full_like(ssy, 0.5)
    return np.clip(u, 0.0, 1.0), np.clip(v, 0.0, 1.0)


def _require_torch():
    try:
        import torch
        from torch_cubic_spline_grids import CubicBSplineGrid2d
    except ImportError as exc:
        raise ImportError(
            "Spline geometry requires torch and torch-cubic-spline-grids. "
            "Install with: pip install torch torch-cubic-spline-grids"
        ) from exc
    return torch, CubicBSplineGrid2d


def _grid_from_data(grid_data, CubicBSplineGrid2d):
    import torch
    data = torch.as_tensor(grid_data, dtype=torch.float32)
    if data.ndim == 2:
        data = data.unsqueeze(0)
    return CubicBSplineGrid2d.from_grid_data(data)


def fit_spline_residual(
    ssx,
    ssy,
    residual,
    resolution=(3, 3),
    n_iter=800,
    lr=0.05,
    weight_decay=1e-3,
    bbox=None,
):
    """Fit CubicBSplineGrid2d to residual z after subtracting the PACEtomo plane."""
    torch, CubicBSplineGrid2d = _require_torch()
    ssx = np.asarray(ssx, dtype=float).ravel()
    ssy = np.asarray(ssy, dtype=float).ravel()
    residual = np.asarray(residual, dtype=float).ravel()
    mask = np.isfinite(ssx) & np.isfinite(ssy) & np.isfinite(residual)
    ssx, ssy, residual = ssx[mask], ssy[mask], residual[mask]
    if ssx.size < 3:
        raise ValueError("Need at least 3 finite points to fit a spline residual.")
    if bbox is None:
        bbox = _bbox_from_points(ssx, ssy)
    res = tuple(int(v) for v in resolution)
    if len(res) != 2:
        raise ValueError("spline resolution must be (n_ssx, n_ssy).")
    u, v = normalize_ss(ssx, ssy, bbox)
    coords = torch.tensor(np.column_stack((u, v)), dtype=torch.float32)
    target = torch.tensor(residual, dtype=torch.float32).reshape(-1, 1)
    grid = CubicBSplineGrid2d(resolution=res, n_channels=1)
    optimizer = torch.optim.Adam(grid.parameters(), lr=lr, weight_decay=weight_decay)
    last_loss = float("nan")
    for _ in range(int(n_iter)):
        optimizer.zero_grad()
        pred = grid(coords)
        loss = torch.mean((pred - target) ** 2)
        loss.backward()
        optimizer.step()
        last_loss = float(loss.detach())
    with torch.no_grad():
        pred = grid(coords).cpu().numpy().ravel()
    rms = float(np.sqrt(np.mean((pred - residual) ** 2)))
    return {
        "resolution": list(res),
        "bbox": {k: float(v) for k, v in bbox.items()},
        "grid_data": grid.data.detach().cpu().numpy().tolist(),
        "rms_um": rms,
        "n_iter": int(n_iter),
        "final_mse": last_loss,
        "n": int(ssx.size),
    }


def evaluate_spline_residual(ssx, ssy, geometry):
    """Spline residual [um] at specimen coordinates."""
    torch, CubicBSplineGrid2d = _require_torch()
    spline = geometry.get("spline") or geometry
    grid = _grid_from_data(spline["grid_data"], CubicBSplineGrid2d)
    ssx_arr = np.asarray(ssx, dtype=float)
    ssy_arr = np.asarray(ssy, dtype=float)
    out_shape = np.broadcast(ssx_arr, ssy_arr).shape
    u, v = normalize_ss(ssx_arr, ssy_arr, spline["bbox"])
    coords = torch.tensor(
        np.column_stack((np.ravel(u), np.ravel(v))), dtype=torch.float32
    )
    with torch.no_grad():
        residual = grid(coords).cpu().numpy().ravel()
    if out_shape == ():
        return float(residual[0])
    return residual.reshape(out_shape)


def evaluate_z(ssx, ssy, geometry):
    """Full height z(SSX, SSY) = plane + spline residual [um]."""
    z_plane = plane_z(ssx, ssy, geometry["pretilt"], geometry["rotation"])
    z_res = evaluate_spline_residual(ssx, ssy, geometry)
    return z_plane + z_res


def build_geometry(ssx, ssy, z, spline_resolution=(3, 3), spline_n_iter=800):
    """Fit plane then spline residual. z should be height relative to the origin."""
    ssx = np.asarray(ssx, dtype=float).ravel()
    ssy = np.asarray(ssy, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    plane = fit_plane(ssx, ssy, z)
    z_pace = plane_z(ssx, ssy, plane["pretilt"], plane["rotation"])
    residual = z - z_pace
    spline = fit_spline_residual(
        ssx, ssy, residual, resolution=spline_resolution, n_iter=spline_n_iter
    )
    z_fit = z_pace + np.asarray(evaluate_spline_residual(ssx, ssy, {"spline": spline}), dtype=float)
    return {
        "model": MODEL_NAME,
        "created": datetime.now().isoformat(timespec="seconds"),
        "pretilt": plane["pretilt"],
        "rotation": plane["rotation"],
        "center": plane["center"],
        "norm": plane["norm"],
        "plane_mean_error": plane["mean_error"],
        "plane_rms_um": plane["rms_um"],
        "spline_rms_um": spline["rms_um"],
        "fit_rms_um": float(np.sqrt(np.mean((z - z_fit) ** 2))),
        "n_points": int(ssx.size),
        "coord_order": ["ssx", "ssy"],
        "spline": spline,
        "points": {
            "ssx": [float(v) for v in ssx],
            "ssy": [float(v) for v in ssy],
            "z_um": [float(v) for v in z],
        },
    }


def save_geometry(path, geometry):
    with open(path, "w") as fh:
        json.dump(geometry, fh, indent=2)


def load_geometry(path):
    with open(path, "r") as fh:
        geometry = json.load(fh)
    if "spline" not in geometry or "grid_data" not in geometry["spline"]:
        raise ValueError(f"Geometry file is missing spline data: {path}")
    return geometry


def make_plot(path, geometry):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    ssx = np.array(geometry["points"]["ssx"], dtype=float)
    ssy = np.array(geometry["points"]["ssy"], dtype=float)
    z = np.array(geometry["points"]["z_um"], dtype=float)
    pretilt = geometry["pretilt"]
    rotation = geometry["rotation"]

    margin = 0.1
    x_min, x_max = float(np.min(ssx)), float(np.max(ssx))
    y_min, y_max = float(np.min(ssy)), float(np.max(ssy))
    x_pad = max((x_max - x_min) * margin, 0.5)
    y_pad = max((y_max - y_min) * margin, 0.5)
    gx = np.linspace(x_min - x_pad, x_max + x_pad, 25)
    gy = np.linspace(y_min - y_pad, y_max + y_pad, 25)
    GX, GY = np.meshgrid(gx, gy)
    Z_plane = plane_z(GX, GY, pretilt, rotation)
    Z_spline = np.asarray(evaluate_z(GX, GY, geometry), dtype=float)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(ssx, ssy, z, s=50, color="#ff0000", label="Geo points")
    for i in range(len(ssx)):
        ax.text(ssx[i], ssy[i], z[i], f" {i + 1}", size=12, color="#ff0000")
    ax.plot_wireframe(GX, GY, Z_plane, color="#888888", linewidth=0.6, label="Plane")
    ax.plot_surface(GX, GY, Z_spline, alpha=0.35, color="#4c78a8", label="Spline")
    ax.set_xlabel("SSX (um)")
    ax.set_ylabel("SSY (um)")
    ax.set_zlabel("z (um)")
    ax.set_title("Lamella geometry (plane + spline)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.show()
