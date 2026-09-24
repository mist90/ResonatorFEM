# ResonatorFEM — Design & Implementation Plan

An edge-element (Whitney/vector) finite-element solver for electromagnetic cavity
resonators, with a Qt/OpenGL GUI. Reworked from the original **MicroLBWave** so that
geometry and meshing are handled by a mature open-source stack (**Gmsh** built on the
**OpenCASCADE** kernel) instead of a hand-written Delaunay triangulator.

## 1. Background

Two earlier attempts live alongside this project:

- **MicroLBWave (v1)** — working FEM solver + GUI, but its mesh is an *unconstrained*
  Bowyer–Watson Delaunay tetrahedralization of a point cloud. Element faces do not
  conform to object surfaces, so the mesh does not faithfully represent the geometry.
- **MicroLBWave2 (v2)** — an attempt at a *constrained conforming* Delaunay mesher
  (robust exact predicates, R-tree index, triangulated surfaces). Correct in ambition
  but too large a problem to finish; it stalled.

**This project keeps the FEM solver and GUI from v1 and replaces the entire geometry +
mesh layer with Gmsh/OpenCASCADE.**

## 2. What we reuse vs. replace

| Keep (from v1)                                             | Replace / add                         |
|-----------------------------------------------------------|---------------------------------------|
| `MicroEngine` — FEM assembly + ARPACK eigen-solve         | `CsgGmshMesher` — CSG geometry + mesh |
| `MicroTetraedr`, `MicroNode`, `MicroEdge` — mesh structs  | `MeshAdapter` — Gmsh → mesh structs   |
| `MicroGrid` — mesh **container** (Delaunay code retired)  | GUI selectors for shape/core          |
| `MathEighValues` (ARPACK++), `MathMatrix*`, `MathComplex` |                                       |
| `MathPoint3D`, `MathVector3D`, `MicroPort`                |                                       |
| `MainDialog`, `MathGrapher`, `GLFunctions` — GUI          |                                       |

## 3. Key architectural insight

`MicroEngine`'s resonator eigen-path never asks *how* the mesh was built. It reads, per
tetrahedron: 4 node coordinates, the 6-edge enumeration, a per-tet material (ε), and
which edges are PEC (perfect electric conductor). So the whole job is a **thin adapter**
that fills the existing `MicroNode`/`MicroTetraedr` structures from Gmsh output.

Better still: PEC walls are detected as **"a tetrahedron face with no neighbour."** If we
mesh the CSG solid `cavity − core`, then *every* exterior face — the outer cavity wall
**and** the inner core surface — is automatically a PEC boundary, with zero manual
tagging. This is exactly the conforming behaviour v1 and v2 lacked.

## 4. The adapter contract (resonator eigen-path)

Verified against `MicroEngine.cpp`, `MicroGrid.cpp`, `MicroTetraedr.cpp`, `MicroNode.cpp`.

**Per node (`MicroNode`)** — required: `_point` (x,y,z) via `setPoint`; flag
`NODE_IS_BOUNDARY` on nodes lying on any boundary face; `NODE_IS_METALL` if on a PEC
surface. `_serialNumber` is assigned later by `setNumberNodesTetraedrs()`. The
`_neighbourNodes` / `_areaFacesDirihle` arrays are **only** used by the Dirichlet path —
**not** the resonator path — so the adapter may leave them empty.

**Per tetrahedron (`MicroTetraedr`)** — required: the 4 node iterators `_nodes[0..3]`
(any orientation; volume uses `fabs(det)`), the 4 neighbour links `_tetraedrs[0..3]`
(neighbour tet or `EmptyIterator()`), and `_serialNumber` (assigned by
`setNumberNodesTetraedrs()`). Circumsphere is **not** needed (Delaunay-only).

**Face-numbering convention (load-bearing).** The engine has two inconsistent surface
conventions. The resonator path (`calculateBoundariesSurface` + `calculateMetallEdges`)
uses **surface `i` = nodes `{i, (i+1)%4, (i+2)%4}`**. The adapter MUST set `_tetraedrs[i]`
to the neighbour across *that* face, so boundary→PEC detection stays self-consistent.
(`getNumPointFromSurface` uses a different convention but is not used in this path.)

**6-edge enumeration** (fixed, `getNumBeginEdge`/`getNumEndEdge`):
`e0:0→1  e1:0→2  e2:0→3  e3:1→2  e4:1→3  e5:2→3`.

**Global assembly** (`calculateResonatorMode`): iterate tets; for each of 6 local edges,
look it up in `globalTable` by undirected node-pair identity (`findEdge`,
`MicroEdge::operator==`). New edge → fresh DOF (`EDGE_FREE`); shared edge → reuse DOF,
flag `EDGE_REVERSE` if opposite orientation (sign flip on assembly); PEC edge →
`EDGE_NULL` (excluded). Result is the generalized eigenproblem **T·x = k²·R·x** solved by
`MathEighValVector` (ARPACK++). T = curl-curl (stiffness), R = ε-weighted mass.

**Materials.** ε is per-tet, indexed by `getSerialNumber()`. In v1 it's assigned by a
point-in-object test; here it comes directly from **Gmsh region (physical volume) tags**.

**Build sequence the adapter must reproduce:** create nodes → create tets w/ node
iterators → build face adjacency to set `_tetraedrs[i]` (convention above) and flag
boundary nodes `NODE_IS_BOUNDARY` → `setNumberNodesTetraedrs()` →
`calculateBoundariesSurface()` per tet → set `SURFACE_IS_METALL` on boundary faces →
`calculateMetallEdges()`.

## 5. Geometry: shapes via OpenCASCADE (through Gmsh's `occ` API)

- **Cavity:** cylinder (`occ::addCylinder`) or rectangular box (`occ::addBox`).
- **Core (kernel):** cylinder, box, or **spiral with elliptical cross-section** — an
  ellipse (`occ::addDisk` scaled, or an ellipse wire) swept along a helix path
  (`occ::addWire` over a helix) via `occ::addPipe`.
- **Domain to mesh:** `occ::cut(cavity, core)` → the vacuum/dielectric region. Mesh with
  `mesh::generate(3)`. Region tags → ε; all exterior faces → PEC.

## 6. Build system

CMake (replaces qmake). Targets:
- `resonator_core` — static lib: FEM + math + mesh adapter (no Qt GUI).
- `resonator_cli` — headless: mesh a cavity, solve, print eigen-frequencies. **The
  milestone-1 validation harness.**
- `ResonatorFEM` — the Qt/OpenGL GUI app.

Dependencies (Ubuntu): `libgmsh-dev`, `liblapack-dev`, `libblas-dev`,
`libeigen3-dev`, `libspectra-dev`, `qtbase5-dev`, `libqt5opengl5-dev`. The
eigensolver uses LAPACK (dense `dsygv`) and Spectra (sparse shift-invert); the
original ARPACK++/SuperLU path was dropped (see `docs/RESULTS.md`).

## 7. Milestones

1. **Empty cylindrical cavity, end-to-end, validated.** Gmsh meshes a cylinder →
   adapter → `MicroEngine` resonator solve → compare lowest eigen-frequencies to the
   analytic TM/TE modes of a cylindrical cavity. Headless `resonator_cli`. *Proves the
   whole numerical pipeline.*
2. **CSG core.** Add `cavity − core` (start: coaxial cylinder core → coaxial resonator,
   also analytic). Confirm PEC-on-core works.
3. **All shapes + GUI.** Cavity {cyl, rect} × core {cyl, rect, spiral/elliptical};
   wire the selectors into `MainDialog`.
4. **Polish for contribution.** Docs, license, examples, CI.

## 8. Validation references (analytic)

- **Cylindrical cavity** radius `a`, height `d`. TM₀₁₀ (dominant, no z-variation):
  `f = c·χ₀₁ / (2π a)`, `χ₀₁ = 2.404826`. General TMₙₘₗ / TEₙₘₗ from Bessel-root
  formulas — used to check the first several computed eigen-frequencies.
