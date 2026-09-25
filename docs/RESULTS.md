# Milestone 1 — Results

**Status: validated.** The reworked pipeline computes physically correct cavity
resonant frequencies end-to-end:

```
Gmsh/OpenCASCADE CSG + mesh  →  MeshAdapter (MicroGrid::loadTetraMesh)
  →  reused edge-element FEM assembly (MicroEngine)  →  LAPACK dsygv eigen-solve
```

## Validation vs. analytic modes

Run `resonator_cli [cyl|box] <dims...> <meshSize>`.

| Cavity                | Mode              | Analytic  | Computed  | Error   |
|-----------------------|-------------------|-----------|-----------|---------|
| Box 1 × 0.8 × 1.5     | TE₁₀₁ (dominant)  | 180.2 MHz | 179.0 MHz | −0.62 % |
| Box                   | TE₀₁₁             | 213 MHz   | 210.8 MHz | −1.0 %  |
| Box                   | TM₁₁₀             | 240 MHz   | 236.8 MHz | −1.3 %  |
| Cylinder r=1, h=2     | TM₀₁₀ (dominant)  | 114.7 MHz | 113.8 MHz | −1.0 %  |
| Cylinder              | TE₁₁₁ (degenerate)| 115.4 MHz | 115.8 MHz | +0.3 %  |
| Cylinder              | TM₀₁₁             | 135.6 MHz | 136.0 MHz | +0.3 %  |

Errors are consistent with first-order edge elements on coarse meshes and shrink
with refinement. This confirms the mesh adapter, PEC boundary detection
("face with no neighbour"), FEM assembly, and eigen-solve are all correct.

## Key engineering decisions made during bring-up

1. **Dropped ARPACK++/SuperLU, use LAPACK `dsygv`.** The vendored ARPACK++ 1.2
   calls an obsolete SuperLU ABI (`dgstrf` lost `drop_tol`, gained `GlobalLU_t*`
   in SuperLU ≥5), which crashed at runtime. A dense generalized symmetric solve
   is robust and exact for the moderate free-edge system, and returns the full
   spectrum including the edge-element null space. See `MathEighValues.cpp`.

2. **Headless core.** Added a `MicroEngine()` constructor (null grapher) and
   `generateFromMesh()` / `calculateSync()` so the solver runs without the GUI.

3. **Adapter face convention.** `MicroGrid::loadTetraMesh` links tetrahedron
   neighbours using surface *i* = nodes `{i,(i+1)%4,(i+2)%4}` — the convention
   `calculateBoundariesSurface`/`calculateMetallEdges` rely on — so PEC detection
   is self-consistent.

## Milestone 2 — CSG core (coaxial cavity)

`resonator_cli coax <rInner> <rOuter> <height> <meshSize>` meshes `outer − inner`
cylinders. The inner cylinder is removed from the domain, so its surface plus the
outer wall and both end caps all become PEC automatically.

| Cavity                        | Mode         | Analytic  | Computed  | Error   |
|-------------------------------|--------------|-----------|-----------|---------|
| Coax rIn=0.3 rOut=1 h=2       | shorted TEM₁ | 75.0 MHz  | 75.4 MHz  | +0.47 % |

The shorted-TEM fundamental `f = c/(2h)` is independent of the radii; the inner
conductor correctly pulls the fundamental down from the empty-cavity TM₀₁₀
(~114 MHz) to 75 MHz. This validates the boolean `cavity − core` cut and
PEC-on-core detection.

## Milestone 3 (partial) — spiral core with elliptical cross-section

`resonator_cli spiral <helixR> <pitch> <turns> <ellA> <ellB> <meshSize>` builds a
cylindrical cavity minus a helical core: an ellipse (semi-axes ellA radial, ellB
axial) swept along an OCC B-spline helix via `occ::addPipe`, then boolean-cut from
the cavity. **The geometry meshes and solves** (e.g. helixR=0.5, pitch=0.5,
turns=2, ellipse 0.12×0.12: 1127 nodes / 4023 tets, solved at n=3546 free edges).

This is the headline shape and the hardest to build — it works. **However**, it
also exposes the eigensolver limitation sharply: on this fully-curved geometry at
a mesh coarse relative to the thin core, the near-zero spurious band rises to
k²≈2, overlapping the expected helical fundamental (~12 MHz quarter-wave
estimate). Reliable spiral frequency extraction therefore needs finer meshes,
which the dense solver cannot reach — making the sparse shift-invert solver
(below) the critical path for this shape. The empty/coaxial cavities validate
cleanly because their spurious band stays well below the physical modes.

Still to add for full shape coverage: rectangular core, and the cavity×core
matrix wired behind one parametric entry point for the GUI.

## Known issues / follow-ups

- **Near-zero spurious band on curved geometry.** The cylinder/coax show a few
  gradient (null-space) modes that don't collapse exactly to zero (k² up to ~1 on
  coarse meshes, always well below the physical modes). The box (flat mesh) has a
  clean null space. Crucially this band **converges toward zero under refinement**
  (cylinder band top: k²≈1.03 at meshSize 0.5 → ≈0.25 at 0.3), confirming it is a
  benign first-order edge-element discretization artifact, not a pipeline bug —
  the discrete gradient kernel is only approximately preserved on distorted
  meshes. A tree-cotree / divergence-cleaning gauge would push it to machine zero.

## Sparse shift-invert eigensolver (Spectra) — done

`MathEighValVectorShiftInvert` (Spectra `SymGEigsShiftSolver`, shift-invert mode)
finds the `nev` eigenvalues nearest a shift σ = k² of a target frequency
(`MicroEngine::setResonatorSolveTarget`). It only factorizes `(A − σB)` (never
`B`), so it tolerates a singular mass matrix, targets the physical range directly
(skipping the null space), and scales far past the dense solver. Robustness:
retries with small multiplicative shift perturbations if a factorization hits an
eigenvalue. The dense `dsygv` path remains as a fallback (σ < 0).

Cylindrical cavity TM₀₁₀ convergence with mesh refinement (now reachable):

| meshSize | tets | TM₀₁₀ error |
|----------|------|-------------|
| 0.35     | 952  | −1.0 %      |
| 0.20     | 3840 | −0.31 %     |
| 0.15     | 8853 | −0.11 %     |

Remaining hard cases (both are the classic "interior eigenvalues near a large
null space" problem, whose clean fix is a **tree-cotree / divergence gauge** that
removes the gradient null space so the spectrum near zero is clean):
- Coax at very fine mesh (0.25): σ near the fundamental lands in a dense
  eigenvalue cluster and `(A − σB)` won't factorize for any nearby shift. Coax
  validates fine at meshSize 0.3.
- Spiral: its helical fundamental sits at very low k² (long wire), overlapping the
  spurious band until the mesh is fine enough to push the null space below it.

## Tree-cotree gauge — attempted, and why it is the wrong tool here

Implemented (`MicroEngine::buildTreeCotreeGauge`, union-find spanning tree with
PEC edges unioned first, non-PEC tree edges gauged to zero) and **tested** —
which revealed it corrupts the eigenproblem, so it is **disabled by default**
(`setGauge(false)`), kept behind a flag for study.

Why it fails: zeroing the tree DOFs is the correct gauge for the *driven/source*
problem `Ax = b` (it fixes the non-unique gradient part of the solution). But the
resonator is a *generalized eigenproblem* `S e = k² T e`. Zeroing tree DOFs
Galerkin-restricts the trial space to `{e : e_tree = 0}`. The physical eigenmodes
are `T`-orthogonal to the gradient null space (`Gᵀ T e = 0`) — they are *not* zero
on the tree edges — so they do not live in that restricted space, and the
restricted eigenvalues are wrong. Measured directly: the cylinder's TM₀₁₀ moved
from the correct 114.6 MHz (ungauged shift-invert) to a spurious ~106–110 MHz
cluster when gauged. The DOF bookkeeping was correct (965 nodes → 360 gauged tree
edges, 3238 cotree DOFs, no disconnection), so this is the method, not a bug.

Correct null-space removal for the eigenproblem is a **mass-metric projection**
(deflate `Range(G)` in the `T` inner product) or a **grad-div penalty**
`(S + s·Cᵀ M_n⁻¹ C) e = k² T e` (lifts the divergent/gradient modes out of the
way while leaving solenoidal physical modes). Both are larger, tunable additions —
the recommended next numerical step for the spiral's very-low fundamental.

Meanwhile the **shift-invert solver (default) already handles standard cavities
well** (cylinder −0.11 %, box, coaxial), which covers most of the shape matrix.

## Grad-div penalty (null-space removal) — working

`MathEighValVectorShiftInvertGauged` implements a **mass-metric grad-div penalty**:
solve `(S + s·(TG) D⁻¹ (GᵀT)) e = k² T e`, `G` = discrete gradient, `D = diag(GᵀTG)`.
Because physical modes satisfy `GᵀT e = 0`, the penalty is **exactly zero on them**
(preserved at any `s`) while the gradient null-space modes are lifted to ~`s`. So the
**physical modes become the smallest eigenvalues** — no target frequency needed, just
take the smallest. (Correct, unlike the Euclidean `GGᵀ` penalty, which perturbs
physical modes, or tree-cotree zeroing, which corrupts the eigenproblem — both above.)

Two details made `G` the true null space (`‖S·G‖/‖S‖ ≈ 0.01–0.07` for cavities):
- **Length scaling** — basis is `len·Whitney`, so a nodal gradient's DOF coefficient
  is `(φ_head − φ_tail)/len` → G entries `∓1/len`.
- **Interior-only columns** — only interior nodes are columns of `G`; a boundary-node
  gradient leaks onto the excluded PEC edges and is not a null vector. Including
  boundary columns was what spoiled the earlier attempts.

ON by default in `resonator_cli` (factor 50); `PENALTY=<factor>` (0 disables) and
`SIGMA=<k²>` are study hooks.

**Results (smallest computed mode = dominant physical mode):**

| Cavity   | Mode  | Analytic  | Computed | Error   |
|----------|-------|-----------|----------|---------|
| Cylinder | TM₀₁₀ | 114.7 MHz | 114.38   | −0.32 % |
| Box      | TE₁₀₁ | 180.2 MHz | 179.28   | −0.49 % |
| Coaxial  | TEM   | 74.9 MHz  | 75.51    | +0.75 % |

**Spiral core:** the fundamental now emerges as the smallest eigenvalue and converges
toward the helical quarter-wave estimate (~11.8 MHz) under refinement — 10.07 MHz at
meshSize 0.18 → 11.87 MHz at 0.13. The gauge residual is larger on the fully-curved
helix (0.39 → 0.27 as the mesh refines: some boundary nodes on the thin helical wall
are imperfectly resolved at coarse meshes) but shrinks with refinement and the
fundamental is physically consistent.

## Dense eigen-solve: robustness and scaling — background

- **Dense eigen-solve: robustness and scaling.** `dsygv` is O(n³) / O(n²) memory
  (n = free edges) and, being Cholesky-based, requires a well-conditioned mass
  matrix. At finer meshes Gmsh can emit sliver tetrahedra whose edges carry ≈0
  mass, making the assembled mass matrix numerically singular; `dsygv` then
  reports "B not positive definite" and the solve is refused (e.g. coax at
  meshSize 0.25, n=1827). Validation therefore runs at 0.3–0.35. The proper fix
  is a sparse **shift-invert** solver (Spectra/SLEPc) targeting the smallest
  nonzero eigenvalues — it tolerates a singular mass matrix via the spectral
  shift and scales to fine meshes. This is the top numerical follow-up.

## Feature completion — done

- **Parametric mesher** `CsgGmshMesher::buildResonator(ResonatorSpec)` covers the
  full matrix: cavity {cylinder, box} × core {none, cylinder, box, spiral}. All 8
  combinations mesh successfully. The core is centred on the cavity axis and
  boolean-cut from the cavity.
- **Qt + VTK GUI** (`MainDialog` + `ResonatorView`): cavity/core dropdowns,
  per-shape parameter fields, mesh size / #modes / penalty controls, a "Compute
  modes" button, and a resonant-frequency list. The hand-rolled `MathGrapher`
  (legacy fixed-function OpenGL) was removed; the engine is now graphics-free
  (`sampleFieldAtCentroids` feeds the view). The **VTK** view (VTK 9,
  `QVTKOpenGLNativeWidget`) has three toggleable layers — **solids** (exterior
  PEC surface, semi-transparent; default on), **fields** (mode field as arrow
  glyphs coloured by magnitude + scalar bar; default on; click a frequency to
  switch modes), **mesh** (tet wireframe; default off). Verified headless under
  Xvfb + software GL: the default cylinder renders TM₀₁₀ with its correct
  on-axis field peak. Env-gated self-test (`RESONATOR_AUTORUN`,
  `RESONATOR_SHOT`) drives compute + screenshot for CI.

## Next milestones

Numerical core and features are in place. Remaining polish:

- Reduce the spiral gauge residual (better boundary-node handling on curved
  walls) and run the solve off the UI thread for large meshes.
- Contribution hygiene: license, example gallery, CI.
