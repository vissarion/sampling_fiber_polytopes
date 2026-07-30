# Locate Eigen and expose it as the imported target secpoly::eigen.
#
# Deliberately does NOT include volesti's own external/cmake-files/Eigen.cmake:
# that file (and LPSolve.cmake) uses the single-argument FetchContent_Populate(<name>)
# form, which is an ERROR on CMake 4.x. volesti is consumed read-only, so we cannot
# fix it there -- we just avoid it.
#
# Resolution order, first hit wins.

if(TARGET secpoly::eigen)
  return()
endif()

set(_eigen_dir "")

# 1. explicit override
if(EIGEN_INCLUDE_DIR AND EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/Core")
  set(_eigen_dir "${EIGEN_INCLUDE_DIR}")
endif()

# 2. the copy volesti already unpacked (the usual hit on a machine with volesti built)
if(NOT _eigen_dir)
  set(_v "${VOLESTI_ROOT}/external/_deps/eigen-src")
  if(EXISTS "${_v}/Eigen/Core")
    set(_eigen_dir "${_v}")
  endif()
endif()

# 3. a system / config-mode Eigen3
if(NOT _eigen_dir)
  find_package(Eigen3 3.3 QUIET NO_MODULE)
  if(Eigen3_FOUND)
    add_library(secpoly::eigen INTERFACE IMPORTED)
    target_link_libraries(secpoly::eigen INTERFACE Eigen3::Eigen)
    message(STATUS "Eigen: system package (Eigen3::Eigen)")
    return()
  endif()
endif()

# 4. fetch, using the MODERN api (never FetchContent_Populate(<name>))
if(NOT _eigen_dir)
  include(FetchContent)
  FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(eigen)
  set(_eigen_dir "${eigen_SOURCE_DIR}")
endif()

if(NOT _eigen_dir)
  message(FATAL_ERROR "Eigen not found. Pass -DEIGEN_INCLUDE_DIR=/path/to/eigen")
endif()

add_library(secpoly::eigen INTERFACE IMPORTED)
target_include_directories(secpoly::eigen SYSTEM INTERFACE "${_eigen_dir}")
message(STATUS "Eigen: ${_eigen_dir}")
