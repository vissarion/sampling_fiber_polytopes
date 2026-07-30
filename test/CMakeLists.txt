add_executable(test_respol_island test_respol_island.cpp)
target_link_libraries(test_respol_island PRIVATE secpoly_optoracle)
add_test(NAME respol_island COMMAND test_respol_island)

add_executable(test_opt_oracle test_opt_oracle.cpp)
target_link_libraries(test_opt_oracle PRIVATE secpoly_optoracle)
add_test(NAME opt_oracle COMMAND test_opt_oracle)
set_tests_properties(opt_oracle PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Compiled WITHOUT the respol include dirs, on purpose: it must succeed using
# only the public headers.
add_executable(test_no_cgal_leak test_no_cgal_leak.cpp)
target_link_libraries(test_no_cgal_leak PRIVATE secpoly_core)
add_test(NAME no_cgal_leak COMMAND test_no_cgal_leak)

# --- structural guards, run as tests -----------------------------------------

# D/CD/PD must have internal linkage in the respol island. respol's own
# res_enum_d.cpp:28-30 declares them non-static, which would export symbols
# literally named D, CD and PD and collide at link time.
add_test(NAME globals_not_exported
         COMMAND ${CMAKE_COMMAND}
                 -DLIB=$<TARGET_FILE:secpoly_optoracle>
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_globals.cmake)

# res_enum_functions.h may be included by exactly one file, and no public header
# may mention CGAL or LEDA.
add_test(NAME source_purity
         COMMAND ${CMAKE_COMMAND}
                 -DROOT=${CMAKE_SOURCE_DIR}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_purity.cmake)

add_executable(test_affine_hull test_affine_hull.cpp)
target_link_libraries(test_affine_hull PRIVATE secpoly_optoracle secpoly_core)
add_test(NAME affine_hull COMMAND test_affine_hull)

add_executable(test_bootstrap test_bootstrap.cpp)
target_link_libraries(test_bootstrap PRIVATE secpoly_optoracle secpoly_core)
add_test(NAME bootstrap COMMAND test_bootstrap)

add_executable(test_separation test_separation.cpp)
target_link_libraries(test_separation PRIVATE secpoly_optoracle secpoly_core)
add_test(NAME separation COMMAND test_separation)
