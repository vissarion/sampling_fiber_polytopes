# Everything CGAL / LEDA / patched-Kernel_d lives here, and every setting is PRIVATE.
#
# This module encapsulates the one translation unit that is allowed to touch respol.
# `respol/include/res_enum_functions.h` can be included exactly ONCE in the whole
# program: it has no include guards, defines non-inline globals (`double conv_time = 0;`
# at :117), reads globals D/CD/PD that the *includer* must define, and #includes
# `tropli/tropli_disc.cpp` -- a .cpp -- at :27.
#
# respol itself is NEVER modified. We replicate its build settings here rather than
# building through its own CMakeLists.

find_package(CGAL REQUIRED COMPONENTS Core)

find_library(LEDA_LIB NAMES leda
             PATHS ${RESPOL_ROOT}/external/leda
             NO_DEFAULT_PATH REQUIRED)
find_library(GMP_LIB  NAMES gmp  REQUIRED)
find_library(MPFR_LIB NAMES mpfr REQUIRED)

# ---------------------------------------------------------------------------
# INCLUDE ORDER IS LOAD-BEARING -- do not "tidy" this list.
#
# respol/respol/CMakeLists.txt:62-69 issues eight successive
# include_directories(BEFORE ...) calls, and each one PREPENDS. The effective
# search order is therefore the REVERSE of its source order, which is what is
# reproduced below.
#
# Why patches/include must come first: both patches/include and
# external/kernel_d/include provide CGAL/Kernel_d/Point_d.h, and ONLY the patched
# copy defines set_index / set_hash. List them the other way round and the
# unpatched header wins, producing ~40 errors of the form
#     "CPoint_d has no member named 'set_index'"
# starting at res_enum_functions.h:370 -- deep inside respol, far from the setting
# that actually caused it.
# ---------------------------------------------------------------------------
set(RESPOL_INCLUDE_DIRS
    ${RESPOL_ROOT}/patches/include                      # MUST be first
    ${RESPOL_ROOT}/external/leda/incl
    ${RESPOL_ROOT}/include
    ${RESPOL_ROOT}/external
    ${RESPOL_ROOT}/external/kernel_d/include
    ${RESPOL_ROOT}/external/extreme_points_d/include
    ${RESPOL_ROOT}/external/spatial_sorting/include
    ${RESPOL_ROOT}/external/triangulation/include)

set(RESPOL_COMPILE_DEFS
    USE_HACKED_KERNEL_ORIENTATION
    USE_HACKED_GAUSSIAN_ELIMINATION
    USE_HASHED_DETERMINANTS
    USE_SORTED_INDICES
    USE_ONLY_CAYLEY_DET_HASH
    COMPUTE_VOL)

# Apply the respol island settings to a target. All PRIVATE: nothing about CGAL,
# LEDA or the patched headers may leak to anything that links this target.
#
# LEDA note: it is only *used* for discriminant polytopes (polytope_type==2, which
# we never set), but it is still required at LINK time. tropli is #included
# unconditionally and LEDA/system/memory_std.h:38 declares a file-scope static
# `memory_manager_init std_memory_mgr_init;` whose constructor lands in .init_array,
# a GC root. Measured: -Wl,--gc-sections cuts 32 undefined LEDA symbols to 4, and no
# linker flag removes those. See DESIGN.md.
function(secpoly_respol_island target)
  target_include_directories(${target} BEFORE PRIVATE ${RESPOL_INCLUDE_DIRS})
  target_compile_definitions(${target} PRIVATE ${RESPOL_COMPILE_DEFS})

  # respol MUST be compiled with assertions off, which is also what its own
  # README mandates (-DCMAKE_BUILD_TYPE=Release).
  #
  # Reason: HashedDeterminantBase::compute_determinant
  # (respol/include/hashed_determinant_base_impl.h:463-476) writes one element
  # past the end of `idx2`:
  #
  #     Index idx2;                                  // ends up with n-1 entries
  #     for (size_t i=1; i<n; ++i) idx2.push_back(idx[i]);
  #     for (size_t i=0; i<n; ++i) { ...; idx2[i]=idx[i]; }   // i==n-1 overflows
  #
  # It is benign in practice -- the last iteration's write is never read back --
  # but it is UB, and _GLIBCXX_ASSERTIONS (on in a Debug build) aborts on it.
  # We cannot fix it: respol is consumed read-only. So we pin the island to
  # assertions-off regardless of the project's build type, which keeps a Debug
  # build of OUR code usable.
  target_compile_options(${target} PRIVATE
                         -U_GLIBCXX_ASSERTIONS -U_GLIBCXX_DEBUG -DNDEBUG)
  target_link_libraries(${target} PRIVATE
                        CGAL::CGAL ${LEDA_LIB} ${GMP_LIB} ${MPFR_LIB})
  target_compile_options(${target} PRIVATE
                         -fvisibility=hidden -fvisibility-inlines-hidden)
endfunction()
