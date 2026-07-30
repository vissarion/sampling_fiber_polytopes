# sampling_secondary

Uniform sampling from the **secondary polytope** Σ(A) of a point configuration
`A = {a₁, …, a_m} ⊂ ℤ^d`, using only a linear-optimization oracle.

Σ(A) is the convex hull of the GKZ vectors of all triangulations of `A`. Its vertices correspond
to the *regular* triangulations of `A` and its dimension is `n = m − d − 1`. The difficulty is
that it is given only **implicitly**: the number of triangulations is exponential in `m`, and
neither an H-representation nor a vertex list can be assumed.

This implements the three-step pipeline of the accompanying draft:

| step | what | how |
|---|---|---|
| (a) | optimization oracle | lifting + regular triangulation, via [`respol`](https://github.com/vissarion/respol) |
| (b) | separation oracle | away-step Frank–Wolfe on the projection problem — the linear subproblem *is* (a) |
| (c) | sampling | ball walk, via [`volesti`](https://github.com/GeomScale/volesti) |

Design rationale, the full Frank–Wolfe construction with its soundness and termination
arguments, and hand-worked numerical examples are in [`DESIGN.md`](DESIGN.md).

## Why sampling rather than enumeration

`respol` can already *enumerate* the vertices of Σ(A), but its beneath-and-beyond driver keeps a
triangulation of the n-dimensional hull in memory, whose full-cell count grows like `N^⌊n/2⌋`.
That caps enumeration at around `n ≈ 10`.

Sampling never builds a hull, so it reaches the regime where enumeration is impossible — in
particular **small d, large m**, because the cost of one oracle call depends on `d` and `m`, not
on `n`. At `d = 2, m = 50` you get `n = 47`: hopeless to enumerate, while each oracle call
triangulates only 52 points in dimension 5.

## Status

Working and validated end to end. `ctest` runs 10 tests in ~4 s.

- Sampling on the draft's Example 1 (`A = {1,2,4,6}`, Σ(A) a quadrilateral in ℝ⁴):
  20 000 samples, **20000/20000 strictly inside**, **χ² = 35.8 on ~26 dof** against an independent
  rejection sampler, mean error 0.11 % of span — from **423 oracle calls** plus 1.5 M lazy cache
  hits.
- Associahedron family (`A` = vertices of a convex n-gon, so Σ(A) is the Stasheff polytope):
  vertex counts match **Catalan(n−2)** exactly for n = 4…7 — 2, 5, 14, 42 — with the right
  dimension `n−3` each time.
- Optimization oracle reproduces the draft's four GKZ vectors for Example 1 exactly, `simple.txt`
  → 4 vertices and `cube3.txt` → 72/74 under random sweeps, matching `respol`'s own test suite.

Not yet implemented: log-concave targets, hit-and-run / billiard walks (they need a boundary
oracle), rounding, and the GLS ellipsoid reference oracle.

## Dependencies and directory layout

The CMake defaults assume this **sibling layout**. A clone with no siblings fails at configure
time, so either reproduce it or override the paths.

```
workspace/
├── volesti/                      ← VOLESTI_ROOT   (default ../../volesti)
├── cgal/                         ← pass -DCGAL_DIR explicitly
└── sampling_fiber_polytopes/
    ├── respol/                   ← RESPOL_ROOT    (default ../respol)
    └── sampling_secondary/       ← this repository
```

All three are overridable cache variables:

```bash
cmake -DRESPOL_ROOT=/path/to/respol \
      -DVOLESTI_ROOT=/path/to/volesti \
      -DCGAL_DIR=/path/to/cgal ..
```

| dependency | version tested | notes |
|---|---|---|
| respol | `develop` @ `de5004c` | must be built first — its `external/leda/libleda.a` is linked |
| volesti | `low_volatility` branch | header-only; Eigen is taken from its `external/_deps/eigen-src` |
| CGAL | 6.2 (header-only) | pass `-DCGAL_DIR` |
| GMP, MPFR, Boost headers | system | `libgmp-dev libmpfr-dev libboost-dev` |
| CMake | ≥ 3.20 | tested with 4.2 |
| C++ | 17 | tested with GCC 15 |

Eigen resolution order: `-DEIGEN_INCLUDE_DIR` → volesti's unpacked copy → system `Eigen3` →
fetch 3.4.0. `lp_solve` is **not** required — the Frank–Wolfe design removes the LP dependency.

## Build and test

```bash
mkdir -p build && cd build
cmake -DCGAL_DIR=/path/to/cgal ..
make -j"$(nproc)"
ctest --output-on-failure
```

Two constraints the build enforces mechanically, both checked by `ctest`:

- **`res_enum_functions.h` may be included by exactly one translation unit.** It has no include
  guards, defines non-inline globals, reads globals the includer must define, and `#include`s a
  `.cpp`. It is therefore isolated in `src/opt_oracle_respol.cpp` behind a PIMPL facade whose
  public header contains no CGAL. Enforced by the `source_purity` and `no_cgal_leak` tests.
- **Include order is load-bearing.** `respol/patches/include` must come *first*: both it and
  `respol/external/kernel_d/include` provide `CGAL/Kernel_d/Point_d.h`, and only the patched copy
  defines `set_index`/`set_hash`. Get it wrong and you see ~40 errors of the form
  *"`CPoint_d` has no member named `set_index`"* deep inside respol. See `cmake/RespolIsland.cmake`.

`respol` and `volesti` are consumed strictly **read-only** — this project never modifies either.
The `source_purity` test enforces that.

## Layout

```
include/secpoly/
  opt_oracle.h          PIMPL facade over respol — zero CGAL in this header
  affine_hull.h         intrinsic chart; keeps volesti in ℝⁿ, never ℝ^m
  rationalize.h         double ↔ exact-integer boundary
  bootstrap.h           draft Lemma 4 — certified relative-interior point
  separation_oracle.h   away-step Frank–Wolfe, with the lazy atom cache
  secondary_body.h      the volesti convex-body concept
src/opt_oracle_respol.cpp   the only TU that includes respol
test/                       10 tests; see DESIGN.md §7 for the worked examples
```

## Licensing

Released under the **GNU Lesser General Public License v3** (see `LICENSE`), matching `respol`
and `volesti`, which this is a derivative work of.

⚠️ **Note on LEDA.** `respol` links the LEDA *Free Edition*, whose terms ship only as
`external/leda/LICENSE.pdf` and normally restrict commercial use. That does not affect the
licensing of the source here, but it does constrain redistribution of built binaries — check
those terms before distributing anything linked against LEDA.

LEDA is pulled in only because `respol` `#include`s `tropli/tropli_disc.cpp` unconditionally and
LEDA installs a file-scope static initializer. The discriminant code path that actually uses it is
never executed here (`polytope_type` is always 1, secondary). See `cmake/RespolIsland.cmake` for
how the dependency could be dropped entirely.
