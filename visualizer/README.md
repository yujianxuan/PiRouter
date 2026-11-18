# Visualization

The `visualizer/` directory contains two independent Matplotlib-based viewers:

- `before_routing_view.py` visualizes a case before routing. It shows chiplet regions, microbumps, package-substrate C4 bumps, and the netlist connections;
- `after_routing_view.py` visualizes a routed result layer by layer. It shows routing segments, vias, original bumps, and escape points.

Both viewers use the benchmark JSON files in `benchmarks/`. The after-routing viewer additionally reads a `routing_result.json` produced by the router. Generated figures are PDF files.

## Requirements

- Python 3;
- Matplotlib.

## Usage

Run the commands below from the project root, where CASE denotes the name of a specific test case.

Visualize a case before routing:

```bash
python visualizer/before_routing_view.py CASE
```

Visualize a routed result:

```bash
python visualizer/after_routing_view.py CASE \
  --result-dir results/CASE
```

PS: The generated PDF files are saved under `visualizer/pic/` by default.