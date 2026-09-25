cmake_minimum_required(VERSION 3.16)
project(ResonatorFEM LANGUAGES C CXX)   # C needed so find_package(MPI) defines MPI::MPI_C (VTK)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Qt5 REQUIRED COMPONENTS Core Widgets)

# --- external libraries -----------------------------------------------------
# Geometry+meshing via Gmsh (OpenCASCADE kernel); eigen-solve via LAPACK (dense)
# and Spectra/Eigen (sparse shift-invert). 3D view via VTK.
find_library(GMSH_LIB   gmsh   REQUIRED)
find_library(LAPACK_LIB lapack REQUIRED)
find_library(BLAS_LIB   blas   REQUIRED)

# --- reusable core: FEM assembly + math + eigen-solve + mesh backend --------
# Depends only on Qt Core (MicroEngine is a QObject); no GUI/GL.
add_library(resonator_core STATIC
  src/MicroNode.cpp
  src/MicroTetraedr.cpp
  src/MicroGrid.cpp
  src/MicroEngine.cpp
  src/MathPoint3D.cpp
  src/MathVector3D.cpp
  src/MathComplex.cpp
  src/MathEighValues.cpp
  src/CsgGmshMesher.cpp
)
target_include_directories(resonator_core PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  /usr/include/eigen3
)
target_link_libraries(resonator_core PUBLIC
  Qt5::Core
  ${GMSH_LIB} ${LAPACK_LIB} ${BLAS_LIB}
)

# --- headless validation harness -------------------------------------------
add_executable(resonator_cli src/cli_resonator.cpp)
target_link_libraries(resonator_cli PRIVATE resonator_core)

# --- Qt + VTK GUI application ----------------------------------------------
# VTK 9.1's packaged config references the MPI::MPI_C imported target, so MPI
# must be found first.
find_package(MPI REQUIRED)
find_package(VTK 9 REQUIRED COMPONENTS
  CommonCore CommonDataModel
  FiltersCore FiltersGeometry FiltersSources
  RenderingCore RenderingOpenGL2 RenderingAnnotation
  InteractionStyle GUISupportQt
)
add_executable(ResonatorFEM
  src/main.cpp
  src/MainDialog.cpp
  src/ResonatorView.cpp
)
target_link_libraries(ResonatorFEM PRIVATE resonator_core Qt5::Widgets ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS ResonatorFEM MODULES ${VTK_LIBRARIES})
