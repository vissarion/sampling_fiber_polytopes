# Sampling from Secondary Polytopes — Design

This document records *why* the implementation looks the way it does. It is written before the
code so the reasoning can be reviewed independently of it. Companion documents:

- `../STATUS.md` — paper-level context, the draft's theorem numbering, open questions
- `../sampling_fiber_polytopes-2.pdf` — the draft this implements
- `../references/` — GKZ94, GLS88, Ziegler95 and the open-access references

---

## 1. The problem

Let `A = {a₁, …, a_m} ⊂ ℝ^d` be a finite point set affinely spanning `ℝ^d`. Its **secondary
polytope** `Σ(A)` is the convex hull of the GKZ vectors of all triangulations of `A`:

> `Φ_T(A) = (φ₁, …, φ_m)`, where `φᵢ = Σ_{σ ∈ T, aᵢ ∈ σ} vol(σ)`

`Σ(A)` has dimension **n = m − d − 1**; its vertices correspond exactly to the *regular*
triangulations of `A` and its edges to bistellar flips (GKZ94 Thm 7.1.7; draft Prop. 1).

The difficulty is that `Σ(A)` is available only **implicitly**. The number of triangulations is
exponential in `m`, and neither an H-representation nor a vertex list can be assumed. The single
primitive we have is a **linear optimization oracle**: given a cost vector, return the GKZ vector
maximizing it. Everything else must be built from that.

**Goal of v1:** sample uniformly from `Σ(A)`. Log-concave targets and hit-and-run come later.

### Why sampling rather than enumeration

`respol` already *enumerates* the vertices of `Σ(A)`, but its beneath-and-beyond driver keeps a
triangulation of the n-dimensional hull in memory — `Triangulation Res(PD)` at
`respol/respol/res_enum_d.cpp:158`, grown by `insert_new_Rvertex2`
(`res_enum_functions.h:1004`). The full-cell count of such a triangulation grows like `N^⌊n/2⌋`,
which caps enumeration around `n ≈ 10`.

Sampling never builds a hull, so it reaches the regime where enumeration is impossible — in
particular **small d, large m**, because the cost of one oracle call depends on `d` and `m`, not
on `n`. At `d = 2, m = 50` we get `n = 47`: hopeless to enumerate, while each oracle call
triangulates only 52 points in dimension 5.

---

## 2. Architecture

```
main_sample.cpp
      |
      v
SecondaryBody          volesti body concept: is_in / dimension / InnerBall   <- volesti BallWalk
      |
      v
SeparationOracle       away-step Frank-Wolfe (+ ellipsoid reference)
      |
      v
SecondaryOptOracle     PIMPL; public header has ZERO CGAL   ---- the only TU touching respol
```

**Hard rule: nothing under `respol/` or `volesti/` is modified.** Both are consumed read-only.
Several choices below exist solely to honour that, and are noted where they arise.

A second constraint is forced by respol itself: `respol/include/res_enum_functions.h` can be
included in **exactly one translation unit**. It has no include guards, defines non-inline globals
(`double conv_time = 0;` at `:117`), reads globals `D`/`CD`/`PD` that the *includer* must define,
and `#include`s `tropli/tropli_disc.cpp` — a `.cpp` — at `:27`. The PIMPL boundary is therefore
structural, not stylistic.

---

## 3. Part A — the optimization oracle

We wrap exactly one respol function (`res_enum_functions.h:839`):

```cpp
template <class NT_> std::vector<NT_> compute_res_vertex2(
    pointset, mi, RD, proj, dets, Pdets, Res, T, nli, conf);
```

`nli` is the cost vector (a `CGAL::Cartesian_d<Gmpq>::Direction_d`); the return value is the GKZ
vector maximizing it. `conf.polytope_type = 1` selects the secondary branch (`:529-536`).

**`Res` and `Pdets` are unused in the body** — in `:865-905` the only occurrences of the token
`Res` are inside the debug string `"new Res vertex (up)"`. Both take empty dummies. This is
precisely what keeps the n-dimensional hull out of memory.

### Why the Cayley trick is still applied

Mathematically `Σ(A)` needs a *single* point set; the Cayley trick belongs to resultants. respol
offers no way to say that. `read_pointset` hardcodes `CD = 2*D+1` (`parse_functions.h:116`) and
**requires exactly d+1 supports**, looping `for(i=0;i<d+1;++i)` at `:124` and calling `exit(-1)` at
`:130-135` otherwise; `res_enum_d.cpp:130` then calls `cayley_trick` unconditionally, even under
`-s` (the block at `:133-146` that would have skipped it is commented out). The shipped
`secondary_examples/cube3.txt`, whose second line is `8 1 1 1` for `d = 3`, confirms this is the
intended usage.

So we encode `Σ(A)` as a **degenerate resultant**:

| quantity | value |
|---|---|
| `pointset` | the `m` points of `A`, followed by `d` copies of the origin |
| `mi` | `[m, 1, 1, …, 1]` (d+1 entries) |
| `proj` | `[0 … m−1]` (kept *settable*, see §8) |
| `D`, `CD`, `PD` | `d`, `2d+1`, `m` |
| `RD` | `(m+d) − 2d − 1 = m − d − 1 = n` ✓ |

then Cayley coordinates are appended (support 0 → zeros, support j → `e_{j−1}`).

**Why this is correct.** With `A₁ = A` and `A₂ = … = A_{d+1} = {0}`, the Minkowski sum is `A`
itself, so regular *mixed* subdivisions of the sum are exactly regular subdivisions of `A`. Every
full-dimensional Cayley simplex must contain all `d` singleton pads — they are the only points with
nonzero Cayley coordinates — leaving exactly `d+1` points of `A`, which is a genuine cell of a
triangulation of `A`, and whose `det` is its normalized volume. Hence `project_upper_hull_r`'s
secondary branch, which sums `det` into each `A`-vertex of each upper simplex, yields precisely the
GKZ vector.

**Known inefficiency, accepted.** This triangulates `m+d` points in dimension `2d+1` where the
mathematics needs only `m` points in dimension `d+1` — mild at `d=2` (5 vs 3), meaningful at `d=5`
(11 vs 6). A non-Cayley entry point would be a worthwhile upstream contribution to respol, but
would require modifying it, so it is out of scope.

We re-implement `cayley_trick`'s ~10 lines in our wrapper (mirroring `:159-185`) — **not** to avoid
the Cayley trick, which we need, but because the original writes `topcom_cayley.txt` into the
current directory at `:188-190`.

### Two numerical facts that shape everything downstream

1. **GKZ vectors are integers.** The `/d!` normalization is commented out (`:570-575`), so entries
   are un-normalized volumes. Integer input gives integer output — exactly representable in
   `double` while `V_norm ≲ 2⁴⁰`, and exactly comparable in `int64` arithmetic for certificates.
2. **Degenerate liftings return an empty vector** (`:518-527`), when the induced regular
   subdivision fails to be a triangulation. Callers *must* detect this; see §5.

---

## 4. Part B — separation from optimization, via Frank–Wolfe

### 4.1 The one idea: projecting onto the body *is* separating from it

Let `x*` be the point of `Σ(A)` nearest to a query `y`, and `δ = ‖y − x*‖`. If `δ = 0` then
`y ∈ Σ(A)`. If `δ > 0`, then

> **γ = y − x\***

is a separating hyperplane normal. Geometrically: stand at the point of the body closest to `y`;
the arrow to `y` points straight out through the surface. A wall perpendicular to that arrow,
slid until it touches at `x*`, has the whole body behind it and `y` in front at distance `δ`.

Algebraically, optimality of `x*` for `min ½‖x − y‖²` over a convex set gives
`(x* − y)ᵀ(x − x*) ≥ 0` for all `x ∈ Σ(A)`, i.e.

> `γᵀx ≤ γᵀx*` for all `x ∈ Σ(A)`, while `γᵀy − γᵀx* = ‖y − x*‖² = δ² > 0`

So `γ`'s maximum over the body is attained at `x*`, and `y` beats it by exactly `δ²`.
**Separation reduces to projection.** What remains is computing the projection with only an
optimization oracle.

### 4.2 Frank–Wolfe does exactly that

Minimize `f(x) = ½‖x − y‖²` over `Σ(A)`, where `∇f(x) = x − y`. Frank–Wolfe linearizes `f` at the
current iterate and minimizes the linear model over the body:

> `s = argmin_{x ∈ Σ(A)} ∇f(x_k)ᵀx = argmax (y − x_k)ᵀx = OPT(γ_k)`, with `γ_k = y − x_k`

— which is **one call to our optimization oracle**, in the direction pointing from where we stand
toward `y`. Then step partway toward it: `x_{k+1} = x_k + α(s − x_k)`, `α ∈ [0,1]`.

Two structural properties make this work. Every iterate stays inside the body automatically — we
start at a vertex and only ever take convex combinations of vertices, so feasibility is free and no
membership test is needed. And the body is touched *only* through `OPT`: no facets, no hull, no
triangulation.

Picture it as an invisible gem. The only question you may ask is "which corner sticks out furthest
in direction γ?" You try to walk to `y` from inside, and each answer either proves `y` is outside
or pulls you closer to it.

### 4.3 The exit test, and why it is exact

**You never need to run Frank–Wolfe to convergence.** After the oracle returns `s`, test

> **`cᵀs < cᵀy` ?**

where `c` is the (rounded, integral) vector actually sent to the oracle — see §5. If yes: for every
`x ∈ Σ(A)`, `cᵀx ≤ cᵀs < cᵀy`, because `s` is the *true* maximizer of `c`. So `c` separates `y`
from the entire body. This is a complete proof standing on its own; it does not care whether `x_k`
was anywhere near the projection, nor how many iterations have run. The oracle's exactness does all
the work.

### 4.4 One identity that makes it transparent

The **Frank–Wolfe gap** is `g = ∇f(x)ᵀ(x − s) = γᵀs − γᵀx ≥ 0`. Substituting,

> `γᵀy − γᵀs = ‖y − x‖² − g`

so

> **the separation test fires ⟺ `g < ‖y − x‖²`**

In words: `‖y − x‖²` is how far you still are from `y`; `g` bounds how much progress remains
available. *Once the progress still possible is less than the distance still remaining, you can
never reach `y`, so `y` is outside.* That is the whole algorithm in one sentence.

### 4.5 The line search

With `d = s − x` and `γ = y − x`, minimizing `φ(α) = ½‖x + αd − y‖²` gives

> **`α* = γᵀd / ‖d‖² = g / ‖d‖²`**, clamped to `[0,1]`

`α*‖d‖` is the length of the shadow `γ` casts on `d`: walk along the segment exactly as far as
`y`'s perpendicular shadow reaches, i.e. stop when `y` is directly beside you, at right angles to
your direction of travel (`φ′(α) = (x(α) − y)ᵀd`, so `φ′ = 0` says exactly that).

The **upper clamp** genuinely binds: the perpendicular foot can lie beyond `s`, meaning `y` is
further along this direction than the body extends; we walk to the vertex and stop. The **lower
clamp never binds** in exact arithmetic, because `s` maximizes `γ` over a body containing `x`, so
`g = γᵀ(s−x) ≥ 0`; it exists only to absorb floating-point noise when `g ≈ 0`.

One step reduces the squared distance by exactly `g²/‖d‖²`. Since `‖d‖ ≤ D` (the diameter), each
step buys at least `g²/D²` — which is where the convergence rate comes from.

### 4.6 Completeness

**Soundness — OUTSIDE is never wrong.** If `y ∈ Σ(A)` then for any direction `c`,
`max_{Σ(A)} cᵀx ≥ cᵀy`, because `y` is itself a competitor. So `cᵀs ≥ cᵀy` always and the test can
never fire. No tolerance, no assumption.

**Termination — if `y` is outside, the test must fire.** Let `δ = dist(y, Σ(A)) > 0`. FW converges,
so `‖y − x_k‖² → δ² > 0` while `g_k → 0`. One limit is strictly positive and the other is zero, so
`g < ‖y − x‖²` holds after finitely many iterations. Geometrically: as `x_k` approaches `x*`, the
direction `γ_k` approaches the outward normal of a supporting hyperplane, and `y` lies strictly on
the far side of it.

**Iteration count.** Vanilla FW has `g_k = O(D²/k)`, firing within `O(D²/δ²)` oracle calls; away
steps give linear convergence, `O(log(D/δ))`. Either way the cost grows as `y` approaches the
boundary — intrinsic, not a defect, since a point *on* `∂Σ(A)` cannot be separated at all.

**If `y` is inside**, `δ = 0` and `‖y − x_k‖ → 0`; we stop at `‖y − x_k‖ ≤ tol` and return `x_k` as
witness — an explicit convex combination of GKZ vectors, hence a genuine point of `Σ(A)` within
`tol` of `y`. This is GLS *weak* membership, exactly the model the draft's Theorem 8 assumes. The
only ambiguous region is the `tol`-thin shell around `∂Σ(A)`, which is unavoidable: deciding exact
membership for a boundary point of an implicitly-defined body requires exact arithmetic throughout.

### 4.7 The loop

```
x <- x0            # bootstrap relative-interior point; DETERMINISTIC (see 4.8)
loop:
    gamma <- y - x                               # true steepest direction, double
    c     <- rationalize(gamma)                  # integer vector actually sent to respol
    s     <- OPT(c)                              # exact argmax of c (perturbed on degeneracy)
    if c^T s < c^T y:                            # UNCONDITIONAL PROOF, zero tolerance
        return OUTSIDE, hyperplane (c, c^T s)
    if ||y - x|| <= tol:
        return INSIDE, witness x
    alpha <- clamp(gamma^T (s-x) / ||s-x||^2, 0, 1)
    if alpha == 0: goto STALL
    x <- x + alpha (s - x)                       # stays inside; no membership test needed
```

### 4.8 Four details that are easy to get wrong

**Test with `c`, never with `γ`.** `s` is the exact maximizer of the *rounded integer* `c` that
respol received, not of the real-valued `γ`. Comparing `γᵀs` against `γᵀy` would quietly reintroduce
a rounding-sized tolerance into what is billed as an unconditional proof. `cᵀs` versus `cᵀy` is
computable in exact integer arithmetic against the stored ambient GKZ vectors, so the OUTSIDE
verdict stays a genuine proof. `γ` remains the right direction for the *line search*, which is pure
optimization and carries no correctness content.

**Determinism.** Always start at `x₀`. Different starting vertices give different — equally valid —
separating hyperplanes and call counts, so the traced examples in §7 are reproducible only with the
start pinned.

**Ties need no rule.** A tie in `argmax cᵀ` means `c` is non-generic — precisely when
`project_upper_hull_r` returns empty. Ties are therefore already routed into the perturbation path
and never reach the FW loop as ambiguity.

**Stall.** `α = 0` means `γᵀ(s−x) ≤ 0`: the rounded direction disagrees with the true steepest
direction by more than the remaining gap. That can only happen once the true gap is below rounding
noise, i.e. at numerical convergence. Policy: re-`rationalize` at +8 bits and retry once; if it
stalls again, decide by `‖y − x‖` against `tol` and record the query as boundary-ambiguous.

### 4.9 Memory, and an overstatement to avoid

Away steps move *away* from the worst-weighted active atom, dropping it at weight zero, which is
what upgrades `O(1/k)` to linear convergence (Lacoste-Julien & Jaggi 2015).

It is tempting to say Carathéodory caps the active set at `n+1`. **It does not.** Carathéodory says
such a representation *exists*; it does not mean away-step FW *maintains* one. Each FW step adds an
atom and only drop steps remove one, so the active set is bounded by the iteration count of the
current query. A hard `n+1` cap needs a fully-corrective variant (re-solving a small QP over the
active set) or an explicit periodic Carathéodory reduction.

Neither is required here: each separation query runs a fresh, short FW loop that is early-exited by
the separation test, so the active set stays small, and it can be capped explicitly if profiling
ever shows growth. State is `O(|S|·n)` plus a bounded atom cache.

The load-bearing claim is unaffected and worth stating precisely: **no hull, no face lattice, no
triangulation at any point.** That — not any particular `n+1` bound — is what removes respol's
`n ≈ 10` ceiling.

### 4.10 Why not the ellipsoid method?

The draft's Theorem 5 obtains separation from optimization via the GLS central-cut ellipsoid
reduction (GLS88 Thm 4.2.2 + Remark 4.2.5 → Cor. 4.2.7). That is what makes Theorem 8 a *theorem*:
it is polynomial, `O(n² log(R/ρ))` oracle calls. But it cannot amortize — every query restarts from
scratch — and the ball walk issues on the order of `10⁶` queries. FW keeps ellipsoid's `O(n²)`
memory profile while reusing work across queries.

The honest trade: **FW is finite but worst-case exponential; the polynomiality claim rests on the
ellipsoid variant.** We therefore ship the ellipsoid oracle too, as a validation-only reference that
must agree with FW on every non-ambiguous query on small instances. GLS make the same
theory-versus-practice point in their own preface: the ellipsoid method "has not proved to be
competitive with the simplex method in practice."

---

## 5. Exact/inexact boundary

`respol` computes in exact rationals; the sampler runs in `double`. The interface is managed as
follows.

**Outgoing (`γ` → integer `c`).** Reduce `γ` modulo `span{1, rows of A}` first — adding an affine
function to a lifting provably does not change the induced regular subdivision, so this is *exactly*
free and only shrinks the dynamic range. Then scale so `max|cᵢ| ≈ 2⁵⁰` and round to `int64`.
Integer lifts keep every determinant integral, which is the fastest path under
`USE_HASHED_DETERMINANTS`; exact dyadic conversion would produce enormous mixed denominators for no
benefit.

**Degeneracy.** On an empty return, perturb `c` by a small *random* integer vector, retrying with
geometrically shrinking magnitude (cap ~30 attempts). Randomize rather than using a fixed direction:
the degenerate set is a finite union of hyperplanes — the walls of the secondary fan — which a
random direction misses with probability 1.

**What perturbation can and cannot affect.** Three consequences, only the third needing care:

1. *Validity is unconditional.* For **any** `c''` whatsoever, `(c'', c''ᵀs)` is an exactly valid
   halfspace for `Σ(A)`, because `s` is the true maximizer of `c''`. Cached halfspaces and OUTSIDE
   proofs never depend on how far the perturbation drifted.
2. *The separation test stays exact*, because we test with the same `c''` we called with.
3. *The convergence rate does depend on it.* The `O(log(D/δ))` away-step rate assumes an **exact**
   linear-minimization oracle for the current direction; a perturbed answer makes it an approximate
   LMO (Jaggi 2013 covers this, with an additive term). The saving grace is that for a sufficiently
   small perturbation the maximizer of `c''` lies in the *optimal face* of `c`, so `cᵀs = max_Σ cᵀ`
   holds exactly and the step is a true FW step. The geometric shrink reaches that regime, and the
   stall guard catches the case where it has not. **The quoted rate is contingent on this** — it is
   not the exact-oracle theorem applied verbatim.

**Incoming (integral `Gmpq` → `double`).** GKZ coordinates are bounded by `V_norm`, and the
intrinsic coordinates are linear combinations, so usable precision degrades as
`ε_mach · V_norm · √m`. The operational threshold is `V_norm < 2⁴⁰`, not `2⁵³`. Compute `V_norm`
exactly as an `mpz` at setup; if it exceeds the threshold, rescale the intrinsic chart globally
(uniform scaling is irrelevant for uniform sampling) while keeping **exact ambient copies** for all
certificates.

---

## 6. Part C — why volesti never sees ℝ^m

`Σ(A)` has dimension `n = m−d−1` but lives in `ℝ^m` — for `A = {1,2,4,6}` it is a *quadrilateral in
ℝ⁴*, of Lebesgue measure zero. Handing volesti an m-dimensional body would not merely be slow:
`BallWalk` proposes `y` in an m-ball and accepts iff `is_in(y)`
(`uniform_ball_walk.hpp:78,82`), which happens with probability **0**. The chain would never move
and would silently return `rnum` copies of the starting point.

The lever is `dimension()`. Every vector volesti creates is sized by it —
`GetPointInDsphere<Point>::apply(P.dimension(), _delta, rng)` at `:78` forwards to `GetDirection`,
which builds `Point p(dim)` (`include/sampling/sphere.hpp:28`). So the body reports
`dimension() = n`, its `Point`s are n-dimensional, and the walk lives entirely in `ℝⁿ`, where
`Q = {y ∈ ℝⁿ : x₀ + Uy ∈ Σ(A)}` is genuinely full-dimensional. Ambient `ℝ^m` appears only inside
`is_in`, when `y` is lifted to `x = x₀ + Uy` for the oracle, and once at the end to report samples
as GKZ vectors.

**Computing `U`.** With `M = [1ᵀ; A_mat]` the `(d+1)×m` constraint matrix, take the *full*
Householder QR of `Mᵀ`; the **last `m−(d+1) = n` columns** of the `m×m` orthogonal factor are an
orthonormal basis of `ker(M)`. Assert `rank(M) = d+1` exactly — if `A` fails to affinely span
`ℝ^d`, then `dim Σ(A) ≠ m−d−1` and the encoding is wrong, so fail loudly.

**Why orthonormal, and why the obvious reason is not the reason.** Uniformity alone does *not*
require it: any affine bijection has constant Jacobian, so uniform pushes forward to uniform under
any linear chart. Orthonormality makes `y ↦ x₀ + Uy` an **isometry** onto the affine hull, which
buys two things a skewed chart loses — the body keeps its true shape (a distorted chart makes it
artificially elongated and the ball walk crawls), and the radius bounds `ρ_simp` and `V√m` transfer
unchanged, both being Euclidean quantities in ambient space.

**The body interface.** `BallWalk` (`uniform_ball_walk.hpp:17-94`) touches exactly three methods,
all of which must be `const` because `apply` takes `BallPolytope const& P` at `:71`:

| method | line | our implementation |
|---|---|---|
| `InnerBall()` → `pair<Point,NT>` | `:67` | `(0ₙ, ρ)` — the chart's origin *is* `x₀` |
| `dimension()` | `:67`, `:78` | **`n = m−d−1`**, never `m` |
| `is_in(y)` → `int`, inside iff `== -1` | `:82` | lift `x = x₀ + Uy`, ask the separation oracle |

The starting point falls out for free: volesti's `Point p(n)` zero-initializes
(`cartesian_geom/point.h:30`), and `y = 0` is exactly `x₀`, which the bootstrap certifies to lie in
`relint Σ(A)`.

**Step size.** `ρ_simp` is a simplex inradius and badly under-estimates the body's, so the derived
`δ = 4ρ/√n` would be far too small. Tune by calling `Walk::update_delta` (`:86`) to drive the
acceptance rate into `[0.2, 0.5]`, then use the explicit-`L` constructor `BallWalk(L)` (`:19`) so
`InnerBall()` is never consulted in production. Acceptance rate is free instrumentation, since
`is_in` is called exactly once per step.

volesti has this same idiom for explicitly-given bodies
(`include/preprocess/get_full_dimensional_polytope.hpp`); we apply it to a body with no
H-representation to reduce.

---

## 7. Worked examples

Using the draft's own Example 1: `A = {1, 2, 4, 6} ⊂ ℝ¹`, so `m = 4`, `d = 1`, `n = 2`.

The four triangulations of `[1,6]` and their GKZ vectors (`φᵢ` = total length of intervals having
`aᵢ` as an endpoint):

| triangulation | GKZ vector |
|---|---|
| `{[1,6]}` | `v₁ = (5, 0, 0, 5)` |
| `{[1,2],[2,6]}` | `v₂ = (1, 5, 0, 4)` |
| `{[1,4],[4,6]}` | `v₃ = (3, 0, 5, 2)` |
| `{[1,2],[2,4],[4,6]}` | `v₄ = (1, 3, 4, 2)` |

All four satisfy `1ᵀx = 10` and `(1,2,4,6)ᵀx = 35`, so `Σ(A)` is a quadrilateral in `ℝ⁴`. Using
`(x₂, x₃)` as chart coordinates the vertices become

> `v₁ = (0,0)`, `v₂ = (5,0)`, `v₃ = (0,5)`, `v₄ = (3,4)`

For *checking only* — the algorithm never sees this — that polygon is
`{u ≥ 0, 2u₁+u₂ ≤ 10, u₁+3u₂ ≤ 15}`.

> **Note.** This `(x₂,x₃)` chart is a coordinate projection, not an isometry. That is fine for
> illustrating separation logic, which is affine-invariant, but the implementation uses the
> orthonormal `U` of §6 for sampling.

### 7.1 A point that is OUTSIDE

Query `y = (1,5)`, i.e. `(2.2, 1, 5, 1.8)` in `ℝ⁴`. Start at `x₀ = v₂ = (5,0)`.

| iter | `γ = y − x` | `OPT(γ)` | `γᵀs` vs `γᵀy` | fires? | `α` | new `x` |
|---|---|---|---|---|---|---|
| 1 | `(−4, 5)` | `v₃ = (0,5)`, val 25 | 25 vs 21 | no | 0.9 | `(0.5, 4.5)` |
| 2 | `(0.5, 0.5)` | `v₄ = (3,4)`, val 3.5 | 3.5 vs 3 | no | 2/13 | `(23/26, 115/26)` |
| 3 | `(3/26)·(1,5)` | `v₃ = (0,5)`, val 75/26 | 75/26 vs 3 | **YES** | — | — |

**Verdict: OUTSIDE**, separating hyperplane `γ ∝ (1,5)`, i.e. in `ℝ⁴`:

> **`x₂ + 5x₃ ≤ 25`**

Verify against all four GKZ vertices — 0, 5, **25**, 23, so the max is 25 — while the query gives
`1 + 25 = 26 > 25`. Genuine separation, margin 1. ✓

Three things worth noticing:

- **It stopped before converging.** The true nearest point is `(0.9, 4.7)` at distance 0.316; we
  were still 0.588 away when the certificate appeared. You never need the projection, only a
  direction good enough to clear the body.
- **The separator is not a facet.** The tight facet is `x₂ + 3x₃ ≤ 15`; FW returned
  `x₂ + 5x₃ ≤ 25`. Any valid separating hyperplane is a complete answer.
- **The trace is start-dependent.** From other starting vertices a different valid hyperplane
  appears — commonly `x₂ + 3x₃ ≤ 15`, the true nearest-edge direction. Only with the start pinned
  is the call count reproducible; tests must assert invariants, not just this trace.

### 7.2 A point that is INSIDE

Query `y = (2,2)`, i.e. `(2.6, 2, 2, 3.4)` in `ℝ⁴`. Start at `x₀ = v₁ = (0,0)`, `tol = 0.01`.

| iter | `γ = y − x` | `OPT(γ)` | `γᵀs` vs `γᵀy` | fires? | `α` | new `x` | `‖y−x‖` |
|---|---|---|---|---|---|---|---|
| 1 | `(2, 2)` | `v₄ = (3,4)`, val 14 | 14 vs 8 | no | 0.560 | `(1.680, 2.240)` | 0.400 |
| 2 | `(0.32, −0.24)` | `v₂ = (5,0)`, val 1.6 | 1.6 vs 0.16 | no | 0.0998 | `(2.011, 2.017)` | 0.0200 |
| 3 | `(−0.011, −0.017)` | `v₁ = (0,0)`, val 0 | 0 vs −0.055 | no | 0.00689 | `(1.997, 2.003)` | 0.00378 |

At the top of iteration 4, `‖y − x₃‖ = 0.00378 ≤ tol` → **INSIDE**.

Notice the test result in every row: `γᵀs ≥ γᵀy` *always*. That is not luck — it is the soundness
argument of §4.6. `y` is in the body, so it competes in the maximization and the maximum can never
fall below `γᵀy`. The OUTSIDE branch is structurally unreachable for an interior point.

The witness accumulates from the step weights:

> `y ≈ 0.400·v₁ + 0.099·v₂ + 0.501·v₄`

Check in `ℝ⁴`: `0.400(5,0,0,5) + 0.099(1,5,0,4) + 0.501(1,3,4,2) = (2.600, 1.998, 2.004, 3.398)`
≈ `(2.6, 2, 2, 3.4)`. ✓

Across both traces the algorithm made **6 oracle calls** and held at most three vertices at a time.
It never enumerated the four triangulations, never built a hull, and never looked at a facet
description — which is the point, since in a real instance those four vertices would instead be an
exponentially large set.

---

## 7bis. The associahedron — a family with closed-form ground truth

GKZ94 Ch. 7 §3B: when `A` is the vertex set of a **convex n-gon**, `Σ(A)` is the
**associahedron** (Stasheff polytope). Its vertices biject with the triangulations of the n-gon,
equivalently with the parenthesizations of `x₁⋯x_{n−1}`, so:

> `dim Σ(A) = m − d − 1 = n − 3`   and   `#vertices = Catalan(n−2)`

This is the strongest validation in the suite, for three reasons: it is **exact and wholly
independent** (Catalan numbers owe nothing to respol, to the chart, or to Frank–Wolfe — Euler
counted these long before secondary polytopes existed); it is a **family**, so an off-by-one or a
boundary case shows up as a pattern rather than as one lucky pass; and it exercises **d = 2**,
whereas Example 1 is d = 1 and cube3 is d = 3.

**Choice of point configuration.** We take the n points on the parabola `y = x²`, i.e.
`aᵢ = (i, i²)` for `i = 0..n−1`. The parabola is strictly convex, so all n points are genuine
vertices of their hull — no three collinear, nothing in the interior — the coordinates stay
integral as respol requires, and index order is already the order around the perimeter, which is
what GKZ's statement assumes.

**Measured** (`test/test_associahedron.cpp`, 3.1 s):

| n-gon | dim Σ(A) | vertices found | Catalan(n−2) | V_norm | degenerate liftings |
|---|---|---|---|---|---|
| 4 | 1 (an interval) | 2 | 2 ✓ | 8 | 0 |
| 5 | 2 (a plane pentagon) | 5 | 5 ✓ | 20 | 0 |
| 6 | 3 | 14 | 14 ✓ | 40 | 0 |
| 7 | 4 | 42 | 42 ✓ | 70 | 0 |

GKZ94 works these small cases explicitly and they agree exactly: n=4 an interval, n=5 a plane
pentagon, n=6 the 3-dimensional polytope with 14 vertices.

The pentagon case is additionally driven through the whole separation stack: 200/200 interior
points are never declared outside (soundness), and 200/200 exterior cuts dominate every vertex of
Σ(A) (each certificate a genuine proof). Notably the atom cache stayed at **3 atoms** throughout —
the bootstrap simplex alone sufficed to separate every far exterior point, so the lazy step made
those queries free.

*Note on trial counts:* the sweep needs more random directions as n grows, because a random
direction hits a given vertex with probability proportional to the solid angle of its normal cone,
and those cones shrink as the vertex count climbs (2 → 42 over this range).

## 8. Scope

**In v1:** uniform sampling by ball walk, over `Σ(A)` given by the respol optimization oracle,
with FW separation and an ellipsoid reference oracle for validation.

**Deliberately deferred:** log-concave targets (draft Remark 3); hit-and-run and billiard, which
need a boundary oracle obtainable by bisection on the separation oracle; rounding and
preprocessing; and general fiber polytopes. `proj` is kept *settable* rather than hardcoded to
`[0…m−1]`, because respol's machinery is already a projection oracle — so projections of `Σ(A)`,
and eventually the fiber-polytope generalization (draft future-work item 1, Ziegler95 Thm 9.6),
need no redesign.

## 9. References

- **[GKZ94]** Gelfand, Kapranov, Zelevinsky, *Discriminants, Resultants, and Multidimensional
  Determinants*, Birkhäuser 1994 — Ch. 7 for secondary polytopes, coherent triangulations, the
  secondary fan, and the associahedron/skew-cube examples.
- **[GLS88]** Grötschel, Lovász, Schrijver, *Geometric Algorithms and Combinatorial Optimization*,
  Springer 1988 — Ch. 2 for the weak/strong oracle hierarchy, Ch. 4 for optimization ↔ separation.
- **[Zie95]** Ziegler, *Lectures on Polytopes*, Springer 1995 — §2.3 polarity, Ch. 9 fiber
  polytopes (Thm 9.6: secondary polytopes are the simplex case).
- **[EFKP13]** Emiris, Fisikopoulos, Konaxis, Peñaranda, IJCGA 23(04n05):397–423, 2013 — respol.
- **[CF+20]** Chalkis, Fisikopoulos et al., *volesti*, arXiv:2007.01578.
- Frank & Wolfe (1956); Lacoste-Julien & Jaggi, *On the global linear convergence of Frank–Wolfe
  optimization variants*, NeurIPS 2015; Jaggi, *Revisiting Frank–Wolfe*, ICML 2013 (approximate
  LMO); Braun, Pokutta & Zink, *Lazifying conditional gradient algorithms*, ICML 2017 (the cached
  weak-separation oracle pattern).
