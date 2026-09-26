# ResonatorFEM

An edge-element (Whitney / vector) finite-element solver for the resonant modes of
electromagnetic cavities, with a Qt/OpenGL GUI. Geometry and meshing use the
open-source **Gmsh** library on its **OpenCASCADE** kernel; the FEM assembly and
GUI are reworked from the original *MicroLBWave* education project, whose
hand-written Delaunay mesher is replaced here by conforming CSG-based meshes.

<img width="1386" height="941" alt="image" src="https://github.com/user-attachments/assets/b3215309-989a-486a-9236-c842846da18a" />


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
`libblas-dev`, `libeigen3-dev`, `libspectra-dev`, `qtbase5-dev`, `libvtk9-dev`,
`libvtk9-qt-dev` (VTK also pulls in `libopenmpi-dev`).

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

## GUI

`ResonatorFEM` is a Qt + **VTK** app: pick the cavity and core shapes from
dropdowns, set the parameters, click **Compute modes**, and read the resonant
frequencies. The parametric backend (`CsgGmshMesher::buildResonator`) covers
cavity {cylinder, box} × core {none, cylinder, box, spiral}.

The VTK view has three independently toggleable layers (checkboxes):
- **Solids** — the geometry (semi-transparent PEC surface), *on* by default;
- **Fields** — the selected mode's field as arrow glyphs coloured by magnitude
  with a scalar bar, *on* by default (click a frequency to switch modes);
- **Mesh** — the tetrahedral wireframe, *off* by default.


## License

MIT — see [LICENSE](LICENSE).
