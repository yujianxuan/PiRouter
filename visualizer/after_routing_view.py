#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import json
import re
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages


FONT_DIR = os.path.join(
    os.path.expanduser("~"), ".local", "share", "fonts", "msttcorefonts"
)
TIMES_NEW_ROMAN_FILES = ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf")
for font_name in TIMES_NEW_ROMAN_FILES:
    font_path = os.path.join(FONT_DIR, font_name)
    if os.path.isfile(font_path):
        font_manager.fontManager.addfont(font_path)

try:
    font_manager.findfont("Times New Roman", fallback_to_default=False)
    PLOT_FONT_FAMILY = "Times New Roman"
except ValueError:
    PLOT_FONT_FAMILY = "DejaVu Serif"

plt.rcParams.update(
    {
        "font.family": PLOT_FONT_FAMILY,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def parse_layer_index(layer_name, top_index=None):
    if isinstance(layer_name, int):
        return layer_name

    layer_str = str(layer_name).strip().lower()

    if layer_str == "bottom":
        return 0

    if layer_str == "top":
        if top_index is None:
            return None
        return top_index

    nums = re.findall(r"\d+", layer_str)
    if nums:
        return int(nums[0])

    return None


def get_layer_num(layer_name, layer_count):
    idx = parse_layer_index(layer_name, top_index=layer_count - 1)
    if idx is None:
        raise ValueError(f"Invalid layer name: {layer_name}")
    return idx


def collect_layer_indices_from_routing(routing_result):
    indices = []
    has_top = False
    max_middle = 0

    if not routing_result:
        return indices, has_top, max_middle

    for segments in routing_result.values():
        for seg in segments:
            for key in ("start_grid_coordinate", "end_grid_coordinate"):
                coord = seg.get(key)
                if not coord or len(coord) < 3:
                    continue

                layer = coord[2]
                if isinstance(layer, str) and layer.strip().lower() == "top":
                    has_top = True
                    continue

                idx = parse_layer_index(layer)
                if idx is not None:
                    indices.append(idx)
                    max_middle = max(max_middle, idx)

    return indices, has_top, max_middle


def collect_layer_indices_from_exit_data(exit_data):
    indices = []

    if not exit_data:
        return indices

    for item in exit_data:
        for key in ("bump_layer", "exit_layer"):
            if key not in item:
                continue
            idx = parse_layer_index(item.get(key))
            if idx is not None:
                indices.append(idx)

    return indices


def collect_layer_indices_from_team_bounds(team_bounds):
    indices = []

    if not team_bounds:
        return indices

    for bounds in team_bounds.values():
        if "layer" in bounds:
            indices.append(bounds["layer"])

    return indices


def infer_layer_count(routing_result=None, exit_data=None, team_bounds=None, user_layer_count=None):
    indices = []

    routing_indices, has_top, max_middle = collect_layer_indices_from_routing(routing_result)
    indices.extend(routing_indices)
    indices.extend(collect_layer_indices_from_exit_data(exit_data))
    indices.extend(collect_layer_indices_from_team_bounds(team_bounds))

    inferred = (max(indices) + 1) if indices else 1

    if has_top:
        inferred = max(inferred, max_middle + 2)

    if user_layer_count is not None:
        inferred = max(inferred, user_layer_count)

    return max(1, inferred)


def load_net_team_mapping(layer_assignment_path):
    net_to_team = {}
    team_bounds = {}
    current_team = None

    with open(layer_assignment_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("PairTeamID:"):
                current_team = int(line.split(":")[1].strip())
                team_bounds[current_team] = {"start": {}, "end": {}}
                continue

            if current_team is None:
                continue

            if line.startswith("Layer:"):
                team_bounds[current_team]["layer"] = int(line.split(":")[1].strip())
            elif line.startswith("NetName:"):
                net_name = line.split(":")[1].strip()
                net_to_team[net_name] = current_team
            elif line.startswith("StartLeft:"):
                team_bounds[current_team]["start"]["left"] = int(line.split(":")[1].strip())
            elif line.startswith("StartRight:"):
                team_bounds[current_team]["start"]["right"] = int(line.split(":")[1].strip())
            elif line.startswith("StartUp:"):
                team_bounds[current_team]["start"]["up"] = int(line.split(":")[1].strip())
            elif line.startswith("StartDown:"):
                team_bounds[current_team]["start"]["down"] = int(line.split(":")[1].strip())
            elif line.startswith("EndLeft:"):
                team_bounds[current_team]["end"]["left"] = int(line.split(":")[1].strip())
            elif line.startswith("EndRight:"):
                team_bounds[current_team]["end"]["right"] = int(line.split(":")[1].strip())
            elif line.startswith("EndUp:"):
                team_bounds[current_team]["end"]["up"] = int(line.split(":")[1].strip())
            elif line.startswith("EndDown:"):
                team_bounds[current_team]["end"]["down"] = int(line.split(":")[1].strip())
            elif line == "---":
                current_team = None

    return net_to_team, team_bounds


def load_exit_terminals(exit_terminals_path):
    if not exit_terminals_path or not os.path.exists(exit_terminals_path):
        return {}, {}, None

    data = load_json(exit_terminals_path)
    net_origin = {}
    net_exit = {}

    for item in data:
        net_name = item.get("net_name")
        if not net_name:
            continue

        origin = (
            item.get("bump_original_x"),
            item.get("bump_original_y"),
            item.get("bump_layer"),
        )
        exit_pt = (
            item.get("exit_x"),
            item.get("exit_y"),
            item.get("exit_layer"),
        )

        net_origin.setdefault(net_name, []).append(origin)
        net_exit.setdefault(net_name, []).append(exit_pt)

    return net_origin, net_exit, data


def layer_title(layer_idx, layer_count):
    if layer_idx == 0:
        return "Bottom Layer"
    if layer_idx == layer_count - 1:
        return "Top Layer"
    return f"M{layer_idx} Layer"


def visualize_case(
    case_name,
    layer_count,
    layout_path,
    netlist_path,
    routing_result_path=None,
    exit_terminals_path=None,
    label_bump=False,
    label_net=False,
    net_to_team=None,
    team_bounds=None,
    team_filter=None,
):
    layout = load_json(layout_path)
    routing_result = None
    if routing_result_path and os.path.exists(routing_result_path):
        routing_result = load_json(routing_result_path)

    net_origin, net_exit, exit_data = load_exit_terminals(exit_terminals_path)

    layer_count = infer_layer_count(
        routing_result=routing_result,
        exit_data=exit_data,
        team_bounds=team_bounds,
        user_layer_count=layer_count,
    )
    grid_info = layout["grid_info"]
    nx = grid_info["grid_max_width"]
    ny = grid_info["grid_max_height"]

    fig_width = max(10.5, 4.2 * layer_count)
    fig_height = 4.8
    fig, axes = plt.subplots(1, layer_count, figsize=(fig_width, fig_height))

    if layer_count == 1:
        axes = [axes]

    fig.suptitle(f"Routing Case: {case_name}", fontsize=9, weight="bold", y=0.98)

    for i in range(layer_count):
        ax = axes[i]

        ax.set_facecolor("white")
        ax.set_aspect("equal")
        ax.set_xlim(0, nx)
        ax.set_ylim(0, ny)
        ax.set_title(layer_title(i, layer_count), fontsize=8, weight="medium")
        ax.set_xlabel("Grid X", fontsize=7)
        ax.set_ylabel("Grid Y", fontsize=7)
        ax.tick_params(axis="both", labelsize=6, width=0.5)
        ax.grid(True, color="#E3E7EA", linewidth=0.35)

    team_colors_all = None
    if net_to_team:
        teams_all = sorted(set(net_to_team.values()))
        cmap = plt.get_cmap("tab20", max(1, len(teams_all)))
        team_colors_all = {team: cmap(i) for i, team in enumerate(teams_all)}

    valid_teams = set(team_colors_all.keys()) if team_colors_all else set()
    if team_filter is not None:
        valid_teams = {team_filter}

    if routing_result:
        route_linewidth = 0.02
        fallback_cmap = plt.get_cmap("tab20", max(1, len(routing_result)))

        for idx, (net_name, segments) in enumerate(routing_result.items()):
            if team_colors_all and net_to_team:
                if net_name not in net_to_team:
                    continue

                team = net_to_team[net_name]
                if team not in valid_teams:
                    continue

                color = team_colors_all[team]
            else:
                color = fallback_cmap(idx)

            for seg in segments:
                start = seg["start_grid_coordinate"]
                end = seg["end_grid_coordinate"]

                start_x, start_y, start_layer = start[0], start[1], start[2]
                end_x, end_y, end_layer = end[0], end[1], end[2]

                start_layer_num = get_layer_num(start_layer, layer_count)
                end_layer_num = get_layer_num(end_layer, layer_count)

                if start_layer_num != end_layer_num:
                    lo = min(start_layer_num, end_layer_num)
                    hi = max(start_layer_num, end_layer_num)

                    for layer_i in range(lo, hi + 1):
                        if layer_i < 0 or layer_i >= layer_count:
                            continue

                        axes[layer_i].scatter(
                            start_x,
                            start_y,
                            marker="s",
                            s=0.003,
                            color="#30B080",
                            edgecolor="none",
                            alpha=0.9,
                            zorder=5,
                        )
                    continue

                if start_layer_num < 0 or start_layer_num >= layer_count:
                    continue

                axes[start_layer_num].plot(
                    [start_x, end_x],
                    [start_y, end_y],
                    linestyle="-",
                    color=color,
                    linewidth=route_linewidth,
                    alpha=0.8,
                    zorder=2,
                )

    if net_to_team and team_colors_all:
        for net_name, team in net_to_team.items():
            if team not in valid_teams:
                continue

            color = team_colors_all[team]

            for x, y, layer in net_origin.get(net_name, []):
                layer_num = get_layer_num(layer, layer_count)
                if layer_num < 0 or layer_num >= layer_count:
                    continue

                axes[layer_num].scatter(
                    x,
                    y,
                    marker="o",
                    s=0.008,
                    facecolor=color,
                    edgecolor="none",
                    zorder=7,
                )

                if label_bump:
                    axes[layer_num].text(
                        x + 0.3,
                        y - 0.3,
                        net_name,
                        fontsize=0.05,
                        color=color,
                        zorder=10,
                    )

            for x, y, layer in net_exit.get(net_name, []):
                layer_num = get_layer_num(layer, layer_count)
                if layer_num < 0 or layer_num >= layer_count:
                    continue

                axes[layer_num].scatter(
                    x,
                    y,
                    marker="^",
                    s=0.01,
                    facecolor=color,
                    edgecolor="none",
                    zorder=7,
                )

                if label_net:
                    axes[layer_num].text(
                        x + 0.3,
                        y + 0.2,
                        net_name,
                        fontsize=0.05,
                        color=color,
                        weight="medium",
                        zorder=10,
                    )

    plt.tight_layout(rect=[0.02, 0.16, 1, 0.94])

    legend_elements = []

    if team_colors_all and valid_teams:
        if team_filter is not None:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    color=team_colors_all[team_filter],
                    linewidth=3,
                    label=f"Team {team_filter}",
                )
            )
        else:
            for team, color in team_colors_all.items():
                legend_elements.append(
                    Line2D([0], [0], color=color, linewidth=3, label=f"Team {team}")
                )

    legend_elements.extend(
        [
            Line2D(
                [0],
                [0],
                marker="s",
                color="w",
                markerfacecolor="#30B080",
                markersize=4,
                label="Via",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="black",
                markersize=5,
                label="Bump",
            ),
            Line2D(
                [0],
                [0],
                marker="^",
                color="w",
                markerfacecolor="black",
                markersize=5,
                label="Exit",
            ),
        ]
    )

    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.00),
        ncol=min(8, len(legend_elements)),
        fontsize=6,
        frameon=False,
    )

    return fig, layer_count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="2.5D routing visualization tool")
    parser.add_argument("case", type=str, help="case name, e.g. C8IO1")
    parser.add_argument(
        "layer",
        nargs="?",
        default=None,
        help="optional layer count; omit it to auto-detect",
    )
    parser.add_argument("--label-bump", action="store_true", help="show bump labels")
    parser.add_argument("--label-net", action="store_true", help="show net labels")
    parser.add_argument(
        "--result-dir",
        type=str,
        default="../output/test",
        help="directory to read routing result files from",
    )
    parser.add_argument(
        "--benchmark-dir",
        type=str,
        default=None,
        help="directory containing *_grid_layout.json and *_netlist.json",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="directory for the generated PDF (default: visualizer/pic)",
    )
    parser.add_argument(
        "--team-pages",
        action="store_true",
        help="also generate one PDF page per team",
    )

    args = parser.parse_args()

    user_layer_count = None
    if args.layer is not None and str(args.layer).lower() != "auto":
        user_layer_count = int(args.layer)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.abspath(os.path.join(script_dir, ".."))
    default_benchmark_dir = os.path.join(project_dir, "benchmarks")
    benchmark_dir = args.benchmark_dir or default_benchmark_dir
    benchmark_dir = os.path.abspath(benchmark_dir)

    if os.path.isabs(args.result_dir):
        result_dir = args.result_dir
    else:
        result_dir = os.path.abspath(os.path.join(script_dir, args.result_dir))

    layout_path = os.path.join(benchmark_dir, f"{args.case}_grid_layout.json")
    netlist_path = os.path.join(benchmark_dir, f"{args.case}_netlist.json")

    routing_result_path = os.path.join(result_dir, "routing_result.json")
    exit_terminals_path = os.path.join(result_dir, f"{args.case}_exit_terminals.json")
    team_file_path = os.path.join(result_dir, f"{args.case}_layer_assignment.txt")

    pic_path = os.path.abspath(args.output_dir or os.path.join(script_dir, "pic"))
    os.makedirs(pic_path, exist_ok=True)
    save_path = os.path.join(pic_path, f"{args.case}_after_routing.pdf")

    net_to_team = None
    team_bounds = None

    if os.path.exists(team_file_path):
        net_to_team, team_bounds = load_net_team_mapping(team_file_path)
        print(f"Loaded team mapping: {len(net_to_team)} nets")
    else:
        print(f"Warning: team mapping file not found: {team_file_path}")

    with PdfPages(save_path) as pdf:
        print("Generating overall visualization...")
        fig_overall, actual_layer_count = visualize_case(
            args.case,
            user_layer_count,
            layout_path,
            netlist_path,
            routing_result_path,
            exit_terminals_path=exit_terminals_path,
            label_bump=args.label_bump,
            label_net=args.label_net,
            net_to_team=net_to_team,
            team_bounds=team_bounds,
            team_filter=None,
        )

        print(f"Auto-detected layer count: {actual_layer_count}")
        pdf.savefig(fig_overall, dpi=300, bbox_inches="tight")
        plt.close(fig_overall)

        if args.team_pages and team_bounds:
            for team_id in sorted(team_bounds.keys()):
                print(f"Generating team {team_id} visualization...")
                fig_team, _ = visualize_case(
                    args.case,
                    actual_layer_count,
                    layout_path,
                    netlist_path,
                    routing_result_path,
                    exit_terminals_path=exit_terminals_path,
                    label_bump=args.label_bump,
                    label_net=args.label_net,
                    net_to_team=net_to_team,
                    team_bounds=team_bounds,
                    team_filter=team_id,
                )
                pdf.savefig(fig_team, dpi=300, bbox_inches="tight")
                plt.close(fig_team)
        elif team_bounds:
            print("Team pages disabled. Use --team-pages to generate per-team pages.")

    print(f"PDF saved to: {save_path}")
