cmake_minimum_required(VERSION 3.20)
project(sampling_secondary CXX)

# --- global, never per-target: libstdc++ ABI must be uniform across every TU ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- dependency roots, all overridable, defaulted to the local checkouts -------
get_filename_component(_here "${CMAKE_CURRENT_SOURCE_DIR}" ABSOLUTE)
set(RESPOL_ROOT  "${_here}/../respol"                CACHE PATH "respol checkout")
set(VOLESTI_ROOT "${_here}/../../volesti"            CACHE PATH "volesti checkout")

foreach(_r RESPOL_ROOT VOLESTI_ROOT)
  if(NOT EXISTS "${${_r}}")
    message(FATAL_ERROR "${_r} not found at ${${_r}} -- pass -D${_r}=/path")
  endif()
endforeach()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindEigenLocal)
include(RespolIsland)

# --- the respol island: one TU, everything about CGAL/LEDA kept PRIVATE -------
add_library(secpoly_optoracle STATIC src/opt_oracle_respol.cpp)
target_include_directories(secpoly_optoracle PUBLIC include)   # <- contains NO CGAL
secpoly_respol_island(secpoly_optoracle)

# --- the CGAL-free core: headers + volesti + Eigen ---------------------------
add_library(secpoly_core INTERFACE)
target_include_directories(secpoly_core INTERFACE
                           include
                           ${VOLESTI_ROOT}/include
                           ${VOLESTI_ROOT}/external)
target_link_libraries(secpoly_core INTERFACE secpoly::eigen)

enable_testing()
add_subdirectory(test)
