#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


# Keep the visual language aligned with after_routing_view.py.
plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

GRID_COLOR = "#E3E7EA"
AUXILIARY_COLOR = "#30B080"
BOUNDARY_COLOR = "#667078"
# Increase/decrease this value to adjust chiplet-region opacity.
CHIPLET_REGION_ALPHA = 0.20
# Bounding-box padding. Lower values make inferred chiplet boxes tighter.
CHIPLET_BOX_GRID_PADDING_RATIO = 0.004
CHIPLET_BOX_SPAN_PADDING_RATIO = 0.02
# Bottom space reserved for the external legend. Lower values move the main
# plot closer to the legend; increase them if a many-row legend overlaps it.
LEGEND_BASE_MARGIN = 0.016
LEGEND_ROW_MARGIN = 0.024
LEGEND_VERTICAL_ANCHOR = 0.010


def load_json(path):
    with open(path, "r", encoding="utf-8") as stream:
        return json.load(stream)


def iter_cases(benchmark_dir):
    for layout_path in sorted(Path(benchmark_dir).glob("*_grid_layout.json")):
        case_name = layout_path.name[: -len("_grid_layout.json")]
        yield case_name, layout_path


def owner_from_name(name):
    """Return the chiplet/IO identifier at the start of a bump name."""
    if not name:
        return "Unknown"

    name = str(name)
    numbered_owner = re.match(r"^(C\d+)(?:_|::)", name)
    if numbered_owner:
        return numbered_owner.group(1)

    hierarchical_owner = re.match(r"^(U::[^:_]+)", name)
    if hierarchical_owner:
        return hierarchical_owner.group(1)

    return name.split("_", 1)[0]


def make_owner_labels(owners, substrate_owners):
    """Keep original owner names, except owners identified on bottom_layer."""
    return {
        owner: "Package substrate" if owner in substrate_owners else owner
        for owner in owners
    }


def collect_layout_points(layout):
    point_by_name = {}
    point_kind = {}
    top_by_owner = defaultdict(list)

    for item in layout.get("top_layer", []):
        name = item.get("bump_name", "")
        point = (item["grid_coord_x"], item["grid_coord_y"])
        point_by_name[name] = point
        point_kind[name] = "top"
        top_by_owner[owner_from_name(name)].append(point)

    bottom_points = []
    for item in layout.get("bottom_layer", []):
        name = item.get("c4_name", "")
        point = (item["grid_coord_x"], item["grid_coord_y"])
        point_by_name[name] = point
        point_kind[name] = "bottom"
        bottom_points.append(point)

    return point_by_name, point_kind, top_by_owner, bottom_points


def collect_connections(netlist, point_by_name, point_kind):
    """Group straight, pre-routing connections by their endpoint owners."""
    grouped_segments = defaultdict(list)
    grouped_points = defaultdict(list)
    skipped = 0

    for net in netlist.get("nets", []):
        bump_names = [
            bump.get("bump_name", "")
            for bump in net.get("bumps", [])
            if bump.get("bump_name")
        ]
        if len(bump_names) != 2 or any(name not in point_by_name for name in bump_names):
            skipped += 1
            continue

        first_name, second_name = bump_names
        pair = tuple(sorted((owner_from_name(first_name), owner_from_name(second_name))))
        first = point_by_name[first_name]
        second = point_by_name[second_name]
        grouped_segments[pair].append((first, second))
        if point_kind[first_name] == "top":
            grouped_points[pair].append(first)
        if point_kind[second_name] == "top":
            grouped_points[pair].append(second)

    return grouped_segments, grouped_points, skipped


def padded_bounds(points, width, height):
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    span_x = max(xs) - min(xs)
    span_y = max(ys) - min(ys)
    pad_x = max(
        width * CHIPLET_BOX_GRID_PADDING_RATIO,
        span_x * CHIPLET_BOX_SPAN_PADDING_RATIO,
    )
    pad_y = max(
        height * CHIPLET_BOX_GRID_PADDING_RATIO,
        span_y * CHIPLET_BOX_SPAN_PADDING_RATIO,
    )
    left = max(0, min(xs) - pad_x)
    right = min(width, max(xs) + pad_x)
    bottom = max(0, min(ys) - pad_y)
    top = min(height, max(ys) + pad_y)
    return left, bottom, right - left, top - bottom


def build_case_figure(case_name, layout, netlist, label_limit=0):
    grid_info = layout["grid_info"]
    width = grid_info["grid_max_width"]
    height = grid_info["grid_max_height"]
    point_by_name, point_kind, top_by_owner, bottom_points = collect_layout_points(layout)
    grouped_segments, grouped_points, skipped = collect_connections(
        netlist, point_by_name, point_kind
    )
    top_by_name = {
        name: point
        for name, point in point_by_name.items()
        if point_kind[name] == "top"
    }

    pairs = sorted(grouped_segments)
    pair_cmap = plt.get_cmap("tab20", max(1, len(pairs)))
    pair_colors = {pair: pair_cmap(index) for index, pair in enumerate(pairs)}

    owners = sorted(top_by_owner)
    all_owners = sorted({owner for pair in pairs for owner in pair} | set(owners))
    substrate_owners = {
        owner_from_name(name)
        for name, kind in point_kind.items()
        if kind == "bottom"
    }
    owner_labels = make_owner_labels(all_owners, substrate_owners)
    owner_cmap = plt.get_cmap("tab20", max(1, len(owners)))
    owner_colors = {owner: owner_cmap(index) for index, owner in enumerate(owners)}

    aspect = width / max(1, height)
    fig_width = 12.0
    fig_height = max(7.0, min(11.0, fig_width / max(0.8, aspect)))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_facecolor("white")
    ax.set_aspect("equal")
    margin_x = width * 0.022
    margin_y = height * 0.022
    ax.set_xlim(-margin_x, width + margin_x)
    ax.set_ylim(-margin_y, height + margin_y)
    ax.grid(True, color=GRID_COLOR, linewidth=0.35, zorder=0)
    ax.tick_params(axis="both", labelsize=8, width=0.5)
    ax.set_xlabel("Grid X", fontsize=10)
    ax.set_ylabel("Grid Y", fontsize=10)
    ax.set_title(
        f"{case_name} (Before Routing)",
        fontsize=15,
        weight="bold",
        pad=24,
    )

    ax.add_patch(
        Rectangle(
            (0, 0),
            width,
            height,
            fill=False,
            edgecolor=BOUNDARY_COLOR,
            linewidth=1.0,
            linestyle=(0, (4, 3)),
            zorder=1,
        )
    )

    # Infer chiplet regions from each chiplet's complete top-bump extent.
    for owner in owners:
        left, bottom, box_width, box_height = padded_bounds(
            top_by_owner[owner], width, height
        )
        color = owner_colors[owner]
        ax.add_patch(
            Rectangle(
                (left, bottom),
                box_width,
                box_height,
                facecolor=color,
                edgecolor=color,
                linewidth=0.8,
                alpha=CHIPLET_REGION_ALPHA,
                zorder=0.5,
            )
        )
        ax.text(
            left + width * 0.006,
            bottom + box_height - height * 0.012,
            owner_labels[owner],
            color=color,
            fontsize=8,
            weight="bold",
            va="top",
            zorder=3,
        )

    # Draw fly-lines efficiently as one collection per chiplet-pair class.
    for pair in pairs:
        color = pair_colors[pair]
        ax.add_collection(
            LineCollection(
                grouped_segments[pair],
                colors=[color],
                linewidths=0.28,
                alpha=0.24,
                zorder=2,
            )
        )
        points = grouped_points[pair]
        ax.scatter(
            [point[0] for point in points],
            [point[1] for point in points],
            s=4.2,
            marker="o",
            facecolor=color,
            edgecolor="none",
            alpha=0.95,
            zorder=5,
        )

    if bottom_points:
        ax.scatter(
            [point[0] for point in bottom_points],
            [point[1] for point in bottom_points],
            s=3.2,
            marker="s",
            facecolor=AUXILIARY_COLOR,
            edgecolor="none",
            alpha=0.85,
            zorder=4,
        )

    if label_limit > 0:
        for index, (name, point) in enumerate(top_by_name.items()):
            if index >= label_limit:
                break
            ax.text(
                point[0] + width * 0.002,
                point[1] + height * 0.002,
                name,
                fontsize=3.5,
                color="#263238",
                zorder=7,
            )

    pair_counts = Counter({pair: len(grouped_segments[pair]) for pair in pairs})
    legend_elements = [
        Patch(
            facecolor=pair_colors[pair],
            edgecolor="none",
            label=(
                f"{owner_labels[pair[0]]} – {owner_labels[pair[1]]} "
                f"({pair_counts[pair]} nets)"
            ),
        )
        for pair in pairs
    ]
    legend_elements.extend(
        [
            Line2D(
                [0], [0], marker="o", color="none", markerfacecolor="#263238",
                markeredgecolor="none", markersize=5, label="microbump"
            ),
            Line2D(
                [0], [0], marker="s", color="none", markerfacecolor=AUXILIARY_COLOR,
                markeredgecolor="none", markersize=5, label="C4 bump"
            ),
        ]
    )
    legend_columns = min(5, len(legend_elements))
    legend_rows = math.ceil(len(legend_elements) / legend_columns)
    bottom_margin = min(
        0.30, LEGEND_BASE_MARGIN + LEGEND_ROW_MARGIN * legend_rows
    )
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, LEGEND_VERTICAL_ANCHOR),
        fontsize=7,
        ncol=legend_columns,
        frameon=False,
    )

    fig.tight_layout(rect=[0.01, bottom_margin, 0.99, 0.99])
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Visualize pre-routing netlist connections and bump placement."
    )
    parser.add_argument(
        "case", nargs="?", default="all",
        help="case name such as C2IO1, or all for every benchmark case"
    )
    parser.add_argument(
        "--benchmark-dir", default=None,
        help="directory containing *_grid_layout.json and *_netlist.json"
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="directory for generated PDF figures; default: visualizer/pic"
    )
    parser.add_argument(
        "--label-limit", type=int, default=0,
        help="label only the first N top bumps; 0 disables labels"
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    default_benchmark_dir = project_dir / "benchmarks"
    benchmark_dir = Path(args.benchmark_dir or default_benchmark_dir)
    output_dir = Path(args.output_dir or script_dir / "pic")
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.case.lower() == "all":
        cases = list(iter_cases(benchmark_dir))
        pdf_path = output_dir / "all_cases_before_routing.pdf"
    else:
        layout_path = benchmark_dir / f"{args.case}_grid_layout.json"
        cases = [(args.case, layout_path)]
        pdf_path = output_dir / f"{args.case}_before_routing.pdf"

    missing = []
    for case_name, layout_path in cases:
        netlist_path = benchmark_dir / f"{case_name}_netlist.json"
        if not layout_path.exists():
            missing.append(str(layout_path))
        if not netlist_path.exists():
            missing.append(str(netlist_path))
    if missing:
        raise FileNotFoundError("Missing benchmark file(s): " + ", ".join(missing))

    with PdfPages(pdf_path) as pdf:
        for case_name, layout_path in cases:
            netlist_path = benchmark_dir / f"{case_name}_netlist.json"
            fig = build_case_figure(
                case_name,
                load_json(layout_path),
                load_json(netlist_path),
                label_limit=args.label_limit,
            )
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)

    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
