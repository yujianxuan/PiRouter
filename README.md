# PiRouter: A Parallel Interposer Router Based on Resource Allocation Driven Escape for 2.5D ICs

Efficient interposer routing is critical since chiplet placement and routing are executed iteratively for layout optimization.
This work presents PiRouter, a deterministic parallel interposer router for 2.5D integrated circuits.
PiRouter forms routing teams using adaptive $k$-means and reduces contention among teams through conflict-aware resource allocation.
Within each team, it uses hybrid escape routing to generate pin-group candidates in parallel, and commits them in a fixed order to ensure
deterministic results.
The routing flow continues with pattern-based
area routing, rerouting, and refinement, using a virtual-node occupancy model to enforce legality on the X-routing grid throughout.
On six public chiplet designs from CIT-Bench, PiRouter achieves speedups of 5.4$\times$ to 36$\times$ over the first-place solution of the 2025 openDACS contest.
It reduces
wirelength and via count on all six designs, with maximum reductions of 17.8\% and 44.6\%, respectively.
With eight threads, PiRouter achieves up to 3.17$\times$ overall speedup and 4.19$\times$ speedup in escape routing over single-threaded execution, while producing identical routing results.

## Directory structure

```text
PiRouter/
├── src/                    C++ source code
├── include/                project header files
├── lib/                    third-party dependencies
├── build/router            compiled executable
├── benchmarks/             test datasets from CIT-Bench
├── visualizer/             visualization tools for layouts before and after routing
├── Makefile                automated build script
└── README.md               quick-start guide and usage
```

## Requirements

For compiling and running the router:

- Linux x86-64;
- GCC/G++ with C++17 and OpenMP support;
- GNU Make.

For the visualization programs:

- Python 3;
- Matplotlib.

## Compilation

From the project root, run:

```bash
make -j
```

The executable is generated at `build/router`. 

To remove the build output:

```bash
make clean
```

The release build uses C++17, OpenMP, optimization, and link-time optimization settings defined in the `Makefile`.

## Execution

The router is invoked directly with one case name, two input files, an output path, and an optional escape-routing thread count:

```bash
./build/router CASE \
  benchmarks/CASE_grid_layout.json \
  benchmarks/CASE_netlist.json \
  OUTPUT/routing_result.json \
  ESCAPE_THREADS
```

PS: `ESCAPE_THREADS` is a positive integer. The main output is `routing_result.json`. The output directory must already exist.

For example:

```bash
./build/router C2IO1 \
  benchmarks/C2IO1_grid_layout.json \
  benchmarks/C2IO1_netlist.json \
  results/C2IO1/routing_result.json \
  8
```

## Dataset

See [benchmarks/README.md](benchmarks/README.md) for the dataset description.

## Visualization

See [visualizer/README.md](visualizer/README.md) for the visualization description.