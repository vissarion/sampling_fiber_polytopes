# sampling_secondary

`sampling_secondary` is a `C++` library for sampling from the secondary polytope
`Sigma(A)` of a point configuration `A` in `R^d`. The vertices of `Sigma(A)` are the regular
triangulations of `A` and its dimension is `m-d-1`, where `m` is the number of points. Neither a
halfspace description nor a vertex list is available in general, since the number of
triangulations is exponential in `m`, so the polytope is accessed exclusively through an
optimization oracle.

The oracle is obtained by using a cost vector as a lifting function: the induced regular
subdivision is computed with [`respol`](https://github.com/vissarion/respol), and its GKZ vector
maximizes the cost over `Sigma(A)`. A separation oracle is derived from it by a conditional
gradient (Frank--Wolfe) method, and the random walks of
[`volesti`](https://github.com/GeomScale/volesti) are run on the resulting body.

Since the cost of one oracle call is governed by `d` and `m` and not by the dimension of
`Sigma(A)`, the method applies to configurations for which enumerating the vertices is out of
reach.

## Compile and use

### Compile dependencies

`respol` must be compiled first: its `external/leda/libleda.a` is linked here. Follow the
instructions in the `respol` repository. `volesti` is header-only and needs no compilation.
Besides those, you need `CGAL` (tested with 6.2, header-only), `GMP`, `MPFR` and the `Boost`
headers. On a Debian or Ubuntu system:

```
$ sudo apt install libgmp-dev libmpfr-dev libboost-dev
```

`Eigen` is taken from the copy `volesti` unpacks under `external/_deps`, from a system
installation, or fetched, in that order. `lp_solve` is not required.

### Directory layout

By default the sibling repositories are expected in the following arrangement:

```
volesti/
cgal/
sampling_fiber_polytopes/
    respol/
    sampling_secondary/
```

Any other arrangement works if the paths are given explicitly:

```
$ cmake -DRESPOL_ROOT=/path/to/respol \
        -DVOLESTI_ROOT=/path/to/volesti \
        -DCGAL_DIR=/path/to/cgal ..
```

### Compile the sources

In a build directory execute:

```
$ cmake -DCGAL_DIR=_YOUR_CGAL_PATH_ ..
$ make
```

The following command will execute the test-suite:

```
$ ctest --output-on-failure
```

Two constraints are enforced by the tests rather than by convention. The header
`res_enum_functions.h` of `respol` may be included by exactly one translation unit, so it is
isolated behind a facade whose public header contains no `CGAL`. The include directories of
`respol` must be given in the order used by its own build, with `patches/include` first;
otherwise the unpatched `CGAL/Kernel_d/Point_d.h` is found and compilation fails inside `respol`
with errors about a missing member `set_index`. Neither `respol` nor `volesti` is modified by
this project.

## Use

A configuration is given by its dimension and its points, with integer coordinates. The oracle,
the affine hull of `Sigma(A)`, a certified relative-interior point and the separation oracle are
built from it, and the body is then handed to a random walk of `volesti`:

```cpp
#include <secpoly/bootstrap.h>
#include <secpoly/opt_oracle.h>
#include <secpoly/secondary_body.h>
#include <secpoly/separation_oracle.h>

using namespace secpoly;

std::vector<std::vector<long long>> A{{1}, {2}, {4}, {6}};   // d = 1, m = 4

SecondaryOptOracle oracle(1, A);
AffineHull          hull(1, A);
BootstrapResult     boot = bootstrap(oracle, hull);

auto sep = std::make_shared<SeparationOracle>(&oracle, hull, boot);
SecondaryBody<Point> body(sep, boot.rho_simp);   // dimension() is m-d-1, not m

uniform_sampling(points, body, rng, BallWalk(L), walk_length, n, start, burns);
```

Samples are produced in coordinates on the affine hull of `Sigma(A)` and are mapped back to GKZ
vectors in `R^m` by `body.to_gkz()`. The optimization oracle may also be used on its own, to
compute the regular triangulation induced by a lifting.

The `proj` argument of `SecondaryOptOracle` selects which GKZ coordinates are kept. It defaults
to all of them; a proper subset gives a projection of `Sigma(A)`.

### Example

For `A = {1,2,4,6}` the secondary polytope is a quadrilateral in `R^4`, whose four vertices are
the GKZ vectors of the four triangulations of the interval `[1,6]`. Running the end-to-end test
prints them together with the corresponding triangulations, and writes the sampled points to
`samples.csv`:

```
(5,0,0,5)   {[1,6]}
(1,5,0,4)   {[1,2],[2,6]}
(3,0,5,2)   {[1,4],[4,6]}
(1,3,4,2)   {[1,2],[2,4],[4,6]}
```

For collinear `A` the secondary polytope is combinatorially an `(m-2)`-cube, so this family
gives instances of any dimension with a known vertex count, useful for validation and for
scaling experiments.

#### Credits

Copyright (c) 2026 Vissarion Fisikopoulos

Built on [`respol`](https://github.com/vissarion/respol) by Ioannis Z. Emiris, Vissarion
Fisikopoulos, Christos Konaxis and Luis Peñaranda, and on
[`volesti`](https://github.com/GeomScale/volesti), part of the
[GeomScale](https://geomscale.github.io) project.

You may redistribute or modify the software under the [GNU Lesser General Public
License](LICENSE) as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version. It is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY.

Note that `respol` links the LEDA Free Edition, whose terms are distributed with it and restrict
some uses. This affects the redistribution of binaries linked against it, not the licensing of
the sources here.

#### Publications

1. Sampling from secondary polytopes. *manuscript*

2. An oracle-based, output sensitive algorithm for projections of resultant polytopes.
   *I.Z. Emiris, V. Fisikopoulos, C. Konaxis, L. Peñaranda.*
   International Journal of Computational Geometry and Applications, vol. 23, pp. 397-423,
   World Scientific, 2013.
   https://vissarion.github.io/publications/EFKP_IJCGA_13.pdf

3. Discriminants, Resultants, and Multidimensional Determinants.
   *I.M. Gelfand, M.M. Kapranov, A.V. Zelevinsky.* Birkhäuser, 1994.
   Chapter 7 for secondary polytopes.
