cmake_minimum_required(VERSION 3.16)
project(ResonatorFEM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets OpenGL)
find_package(OpenGL REQUIRED)

# --- external libraries -----------------------------------------------------
# Geometry+meshing via Gmsh (OpenCASCADE kernel); dense generalized symmetric
# eigen-solve via LAPACK. (The old ARPACK++/SuperLU path was dropped: the
# vendored ARPACK++ 1.2 is ABI-incompatible with modern SuperLU.)
find_library(GMSH_LIB   gmsh   REQUIRED)
find_library(LAPACK_LIB lapack REQUIRED)
find_library(BLAS_LIB   blas   REQUIRED)

# --- reusable core: FEM assembly + math + eigen-solve + mesh backend -------
# (still depends on Qt because MicroEngine is a QThread/QObject; a fully
#  Qt-free numeric core is planned — see docs/PLAN.md milestone 4.)
add_library(resonator_core STATIC
  src/MicroNode.cpp
  src/MicroTetraedr.cpp
  src/MicroGrid.cpp
  src/MicroEngine.cpp
  src/MicroPort.cpp
  src/MathPoint3D.cpp
  src/MathVector3D.cpp
  src/MathComplex.cpp
  src/MathMatrix.cpp
  src/MathMatrixSparse.cpp
  src/MathEighValues.cpp
  src/MathObject.cpp
  src/MathGeomObjects.cpp
  src/MathGrapher.cpp
  src/GLFunctions.cpp
  src/CsgGmshMesher.cpp
)
target_include_directories(resonator_core PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  /usr/include/eigen3
)
target_link_libraries(resonator_core PUBLIC
  Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL
  OpenGL::GL OpenGL::GLU
  ${GMSH_LIB} ${LAPACK_LIB} ${BLAS_LIB}
)

# --- headless validation harness (milestone 1) ----------------------------
add_executable(resonator_cli src/cli_resonator.cpp)
target_link_libraries(resonator_cli PRIVATE resonator_core)

# --- Qt/OpenGL GUI application ---------------------------------------------
add_executable(ResonatorFEM
  src/main.cpp
  src/MainDialog.cpp
)
target_link_libraries(ResonatorFEM PRIVATE resonator_core)
