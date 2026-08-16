#!Python
# ===================================================================
# ScriptName     decolace_LPP.test_stats
# Purpose:       Test-mode CTF statistics: center vs others, image
#                shift vs defocus error and astigmatism.
# ===================================================================

import csv
import os

import numpy as np


def _finite(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return np.nan
    return v if np.isfinite(v) else np.nan


def _json_num(value):
    v = _finite(value)
    return None if not np.isfinite(v) else float(v)


def mark_center_holes(holes):
    """Mark the hole nearest SS (0, 0) as is_center."""
    if not holes:
        return None
    dists = [float(np.hypot(float(h["ssx"]), float(h["ssy"]))) for h in holes]
    center_idx = int(np.argmin(dists))
    for i, hole in enumerate(holes):
        hole["is_center"] = i == center_idx
    return center_idx


def attach_hole_test(hole, result, science_target, probe_target):
    """Store per-hole test measurements (JSON-safe numbers)."""
    science = result.get("science_ctf") or {}
    probe = result.get("ctf") or {}
    sci_df = _finite(science.get("defocus_um"))
    prb_df = _finite(probe.get("defocus_um"))
    sci_err = sci_df - float(science_target) if np.isfinite(sci_df) else np.nan
    prb_err = prb_df - float(probe_target) if np.isfinite(prb_df) else np.nan
    hole["test"] = {
        "is_x": _json_num(result.get("is_x")),
        "is_y": _json_num(result.get("is_y")),
        "is_mag": _json_num(result.get("is_mag")),
        "science_target_um": float(science_target),
        "probe_target_um": float(probe_target),
        "science_defocus_um": _json_num(sci_df),
        "science_defocus_err_um": _json_num(sci_err),
        "science_astig_um": _json_num(science.get("astig_um")),
        "science_astig_x_um": _json_num(science.get("astig_x_um")),
        "science_astig_y_um": _json_num(science.get("astig_y_um")),
        "science_score": _json_num(science.get("fit_score")),
        "science_resolution_A": _json_num(science.get("resolution_A")),
        "probe_defocus_um": _json_num(prb_df),
        "probe_defocus_err_um": _json_num(prb_err),
        "probe_astig_um": _json_num(probe.get("astig_um")),
        "probe_astig_x_um": _json_num(probe.get("astig_x_um")),
        "probe_astig_y_um": _json_num(probe.get("astig_y_um")),
        "probe_score": _json_num(probe.get("fit_score")),
        "probe_resolution_A": _json_num(probe.get("resolution_A")),
    }
    return hole["test"]


def _group_holes(holes):
    acquired = [h for h in holes if h.get("acquired") and h.get("test")]
    center = [h for h in acquired if h.get("is_center")]
    others = [h for h in acquired if not h.get("is_center")]
    return {"all": acquired, "center": center, "others": others}


def _col(holes, key):
    return np.array([_finite((h.get("test") or {}).get(key)) for h in holes], dtype=float)


def _stats(values):
    v = np.asarray(values, dtype=float)
    finite = v[np.isfinite(v)]
    n = int(v.size)
    n_ok = int(finite.size)
    if n_ok == 0:
        return {
            "n": n,
            "n_ok": 0,
            "n_failed": n,
            "mean": np.nan,
            "std": np.nan,
            "rms": np.nan,
        }
    return {
        "n": n,
        "n_ok": n_ok,
        "n_failed": n - n_ok,
        "mean": float(np.mean(finite)),
        "std": float(np.std(finite, ddof=1)) if n_ok > 1 else 0.0,
        "rms": float(np.sqrt(np.mean(finite ** 2))),
    }


def _pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(mask)) < 2:
        return np.nan
    if np.std(x[mask]) == 0 or np.std(y[mask]) == 0:
        return np.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def summarize(holes):
    groups = _group_holes(holes)
    summary = {"groups": {}, "correlations": {}}
    for name, subset in groups.items():
        is_mag = _col(subset, "is_mag")
        block = {"n_holes": len(subset), "mean_is_mag": _stats(is_mag)["mean"]}
        for shot in ("science", "probe"):
            err = _col(subset, f"{shot}_defocus_err_um")
            astig = np.abs(_col(subset, f"{shot}_astig_um"))
            err_stats = _stats(err)
            astig_stats = _stats(astig)
            block[shot] = {
                "n": err_stats["n"],
                "n_failed": err_stats["n_failed"],
                "defocus_err_mean": err_stats["mean"],
                "defocus_err_std": err_stats["std"],
                "defocus_err_rms": err_stats["rms"],
                "astig_mean": astig_stats["mean"],
                "astig_std": astig_stats["std"],
            }
        summary["groups"][name] = block

    all_holes = groups["all"]
    is_mag = np.abs(_col(all_holes, "is_mag"))
    summary["correlations"] = {
        "is_vs_science_defocus_err": _pearson(is_mag, np.abs(_col(all_holes, "science_defocus_err_um"))),
        "is_vs_probe_defocus_err": _pearson(is_mag, np.abs(_col(all_holes, "probe_defocus_err_um"))),
        "is_vs_science_astig": _pearson(is_mag, np.abs(_col(all_holes, "science_astig_um"))),
        "is_vs_probe_astig": _pearson(is_mag, np.abs(_col(all_holes, "probe_astig_um"))),
    }
    return summary


def _fmt(value, digits=4):
    v = _finite(value)
    return "nan" if not np.isfinite(v) else f"{v:.{digits}f}"


def format_summary_table(summary):
    lines = [
        "group    shot     n  fail  df_err_mean  df_err_std  df_err_rms  astig_mean  astig_std  mean_|IS|",
    ]
    for group_name, block in summary["groups"].items():
        for shot in ("science", "probe"):
            s = block[shot]
            lines.append(
                f"{group_name:<8} {shot:<8} {s['n']:3d} {s['n_failed']:5d}  "
                f"{_fmt(s['defocus_err_mean']):>11} {_fmt(s['defocus_err_std']):>11} "
                f"{_fmt(s['defocus_err_rms']):>10} {_fmt(s['astig_mean']):>10} "
                f"{_fmt(s['astig_std']):>9} {_fmt(block['mean_is_mag']):>9}"
            )
    corr = summary["correlations"]
    lines.append(
        "corr |IS| vs |science df err|={0}  |probe df err|={1}  "
        "|science astig|={2}  |probe astig|={3}".format(
            _fmt(corr["is_vs_science_defocus_err"]),
            _fmt(corr["is_vs_probe_defocus_err"]),
            _fmt(corr["is_vs_science_astig"]),
            _fmt(corr["is_vs_probe_astig"]),
        )
    )
    return "\n".join(lines)


def _measurement_rows(holes):
    rows = []
    for hole in holes:
        if not hole.get("acquired") or not hole.get("test"):
            continue
        t = hole["test"]
        rows.append(
            {
                "index": hole.get("index"),
                "ssx": hole.get("ssx"),
                "ssy": hole.get("ssy"),
                "is_center": bool(hole.get("is_center")),
                "is_x": t.get("is_x"),
                "is_y": t.get("is_y"),
                "is_mag": t.get("is_mag"),
                "science_target_um": t.get("science_target_um"),
                "science_defocus_um": t.get("science_defocus_um"),
                "science_defocus_err_um": t.get("science_defocus_err_um"),
                "science_astig_um": t.get("science_astig_um"),
                "science_score": t.get("science_score"),
                "science_resolution_A": t.get("science_resolution_A"),
                "probe_target_um": t.get("probe_target_um"),
                "probe_defocus_um": t.get("probe_defocus_um"),
                "probe_defocus_err_um": t.get("probe_defocus_err_um"),
                "probe_astig_um": t.get("probe_astig_um"),
                "probe_score": t.get("probe_score"),
                "probe_resolution_A": t.get("probe_resolution_A"),
            }
        )
    return rows


def write_measurements_csv(path, holes):
    rows = _measurement_rows(holes)
    if not rows:
        return path
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_summary_csv(path, summary):
    rows = []
    for group_name, block in summary["groups"].items():
        for shot in ("science", "probe"):
            s = block[shot]
            rows.append(
                {
                    "group": group_name,
                    "shot": shot,
                    "n": s["n"],
                    "n_failed": s["n_failed"],
                    "defocus_err_mean": s["defocus_err_mean"],
                    "defocus_err_std": s["defocus_err_std"],
                    "defocus_err_rms": s["defocus_err_rms"],
                    "astig_mean": s["astig_mean"],
                    "astig_std": s["astig_std"],
                    "mean_is_mag": block["mean_is_mag"],
                }
            )
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    return path


def make_plot(path, holes, summary):
    import matplotlib.pyplot as plt

    acquired = [h for h in holes if h.get("acquired") and h.get("test")]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), tight_layout=True)

    ax = axes[0, 0]
    if acquired:
        ssx = np.array([h["ssx"] for h in acquired], dtype=float)
        ssy = np.array([h["ssy"] for h in acquired], dtype=float)
        err = _col(acquired, "science_defocus_err_um")
        sc = ax.scatter(ssx, ssy, c=err, cmap="coolwarm", s=70)
        fig.colorbar(sc, ax=ax, shrink=0.8, label="science df err (um)")
        for h in acquired:
            if h.get("is_center"):
                ax.scatter([h["ssx"]], [h["ssy"]], marker="*", s=180, color="k", label="center")
    ax.plot(0, 0, "k+", markersize=10)
    ax.set_xlabel("SSX (um)")
    ax.set_ylabel("SSY (um)")
    ax.set_title("Science defocus error")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    is_mag = np.abs(_col(acquired, "is_mag"))
    sci_err = np.abs(_col(acquired, "science_defocus_err_um"))
    prb_err = np.abs(_col(acquired, "probe_defocus_err_um"))
    ax.scatter(is_mag, sci_err, label="science", s=40)
    ax.scatter(is_mag, prb_err, label="probe", s=40, marker="s")
    for h in acquired:
        if h.get("is_center"):
            t = h["test"]
            ax.scatter(
                [abs(_finite(t.get("is_mag")))],
                [abs(_finite(t.get("science_defocus_err_um")))],
                marker="*",
                s=160,
                color="k",
                zorder=5,
                label="center",
            )
    ax.set_xlabel("|image shift|")
    ax.set_ylabel("|defocus error| (um)")
    ax.set_title("|IS| vs |defocus error|")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    sci_as = np.abs(_col(acquired, "science_astig_um"))
    prb_as = np.abs(_col(acquired, "probe_astig_um"))
    ax.scatter(is_mag, sci_as, label="science", s=40)
    ax.scatter(is_mag, prb_as, label="probe", s=40, marker="s")
    for h in acquired:
        if h.get("is_center"):
            t = h["test"]
            ax.scatter(
                [abs(_finite(t.get("is_mag")))],
                [abs(_finite(t.get("science_astig_um")))],
                marker="*",
                s=160,
                color="k",
                zorder=5,
                label="center",
            )
    ax.set_xlabel("|image shift|")
    ax.set_ylabel("|astig| (um)")
    ax.set_title("|IS| vs |astig|")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.axis("off")
    ax.set_title("Summary")
    ax.text(0.0, 1.0, format_summary_table(summary), va="top", ha="left", family="monospace", fontsize=7)

    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_outputs(holes, out_dir, stem, log=None, echo_table=True):
    mark_center_holes(holes)
    summary = summarize(holes)
    meas_path = os.path.join(out_dir, f"{stem}_test_measurements.csv")
    sum_path = os.path.join(out_dir, f"{stem}_test_summary.csv")
    plot_path = os.path.join(out_dir, f"{stem}_test_stats.png")
    write_measurements_csv(meas_path, holes)
    write_summary_csv(sum_path, summary)
    try:
        make_plot(plot_path, holes, summary)
    except Exception as exc:
        if log is not None:
            log(f"WARNING: Could not write test-mode plot: {exc}")
        plot_path = None
    table = format_summary_table(summary)
    if log is not None and echo_table:
        log("----- test_mode summary -----")
        log(table)
        log(f"Saved {meas_path}")
        log(f"Saved {sum_path}")
        if plot_path:
            log(f"Saved {plot_path}")
    return summary
