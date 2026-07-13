# One-way function via boson sampling — reproducible figures

Scripts and cached data reproducing the figures and speed tests of the paper
*Classical simulability of one-way functions based on boson sampling*.

The estimator itself lives in the package at
[`src/boson_samplers/one_way_function.jl`](../../../src/boson_samplers/one_way_function.jl)
(public entry points: `find_most_probable_bin`, `estimate_S`, `estimate_time`,
`estimate_time_to_precision`). The floating-point precision model referenced by
the `F64-<X>` tags in that source is documented in
[`numerical_precision.md`](numerical_precision.md).

All scripts run against the package environment. Invoke them from the package
root; each `cd`s to this directory (or reads/writes via `@__DIR__`) so outputs
land here regardless of launch directory:

```bash
# main paper figure (4 panels); reuses paper_figure_data.jls cache by default
julia --project=. docs/publication/one_way_function/paper_figure.jl
#   PAPER_FIGURE_RECOMPUTE=1 forces a fresh simulation.

# variance-scaling figure (Var(Z) vs n, power-law fit)
JULIA_NUM_THREADS=8 julia --project=. docs/publication/one_way_function/variance_scaling_clean.jl

# speed tests (print-only, no file output)
JULIA_NUM_THREADS=8 julia --project=. docs/publication/one_way_function/benchmark_cdf_value.jl
JULIA_NUM_THREADS=8 julia --project=. docs/publication/one_way_function/benchmark_fixed_error.jl
JULIA_NUM_THREADS=24 julia --project=. docs/publication/one_way_function/check_n100_1percent.jl
```

## Files

| File | Produces |
|---|---|
| `paper_figure.jl` | `paper_figure.pdf` / `.png` — the 4-panel paper figure |
| `paper_figure_data.jls` | serialized simulation cache for `paper_figure.jl` |
| `variance_scaling_clean.jl` | `variance_scaling_clean.pdf` / `.png` and `variance_vs_n.pdf` / `.png` |
| `benchmark_cdf_value.jl` | wall-clock cost of one `S(x0)` value vs `M` and `m·n` |
| `benchmark_fixed_error.jl` | cost of one `S(x0)` value at fixed statistical error `ε` |
| `check_n100_1percent.jl` | verified cost of one `S(x0)` value at `n=m=100`, `ε=1%` |
| `numerical_precision.md` | precision analysis of the `Float64` operations |
