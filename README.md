# ResonatorFEM

An edge-element (Whitney / vector) finite-element solver for the resonant modes of
electromagnetic cavities, with a Qt/OpenGL GUI. Geometry and meshing use the
open-source **Gmsh** library on its **OpenCASCADE** kernel; the FEM assembly and
GUI are reworked from the original *MicroLBWave* education project, whose
hand-written Delaunay mesher is replaced here by conforming CSG-based meshes.

## Capabilities

- **Cavity shapes:** cylindrical, rectangular.
- **Cores (kernels):** cylindrical (coaxial), and a **helical spiral with an
  elliptical cross-section** (ellipse swept along a helix, boolean-cut from the
  cavity). Every exterior face of the meshed `cavity − core` solid becomes a PEC
  wall automatically.
- **Solvers:** dense LAPACK `dsygv` (small problems, full spectrum) and a sparse
  **Spectra shift-invert** solver (scales to fine meshes, targets a frequency).

## Build

Requires a C++17 compiler, CMake, and: `libgmsh-dev`, `liblapack-dev`,
`libblas-dev`, `libeigen3-dev`, `libspectra-dev`, `qtbase5-dev`,
`libqt5opengl5-dev`.

```sh
cmake -S . -B build -G Ninja
cmake --build build
```

Targets: `resonator_core` (library), `resonator_cli` (headless validation
harness), `ResonatorFEM` (GUI).

## Validation

`resonator_cli` meshes a cavity, runs the eigen-solve, and compares to analytic
modes:

```sh
./build/resonator_cli cyl  1 2 0.2         # cylindrical cavity  (TM010 -0.3%)
./build/resonator_cli box  1 0.8 1.5 0.25  # rectangular cavity  (TE101 -0.6%)
./build/resonator_cli coax 0.3 1 2 0.3     # coaxial cavity      (TEM   +0.5%)
./build/resonator_cli spiral 0.5 0.5 2 0.12 0.12 0.18   # spiral core
```

See `docs/PLAN.md` (design + FEM↔mesh contract) and `docs/RESULTS.md` (validation
results, numerical findings, and known limitations).

## GUI

`ResonatorFEM` is a Qt/OpenGL app: pick the cavity and core shapes from dropdowns,
set the parameters, click **Compute modes**, and read the resonant frequencies
while the mesh renders in the 3D view. The parametric backend
(`CsgGmshMesher::buildResonator`) covers cavity {cylinder, box} × core {none,
cylinder, box, spiral}.

## Status

Working end to end: geometry → Gmsh/OCC mesh → edge-element FEM → grad-div-gauged
sparse eigensolve → GUI. Cavities validated against analytic modes (cyl −0.3%,
box −0.5%, coax +0.8%); the spiral fundamental converges to the helical
quarter-wave estimate. Remaining polish (spiral gauge residual on curved walls,
off-thread solve, license/CI) is tracked in `docs/RESULTS.md`.
