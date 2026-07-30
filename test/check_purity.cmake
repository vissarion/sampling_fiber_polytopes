# Structural invariants that the compiler cannot express.
#
# 1. res_enum_functions.h may be included by EXACTLY ONE translation unit. It has
#    no include guards, defines non-inline globals, and #includes a .cpp, so a
#    second includer is an ODR violation / duplicate-symbol link error.
# 2. No public header may mention CGAL or LEDA -- that is the PIMPL boundary.

file(GLOB_RECURSE sources "${ROOT}/src/*.cpp" "${ROOT}/test/*.cpp"
                          "${ROOT}/include/*.h")
set(includers "")
foreach(f ${sources})
  file(READ "${f}" text)
  # Match an actual #include directive, not a mention in a comment.
  if(text MATCHES "#include[ \t]*[<\"]res_enum_functions\\.h[>\"]")
    list(APPEND includers "${f}")
  endif()
endforeach()

list(LENGTH includers n)
if(NOT n EQUAL 1)
  message(FATAL_ERROR
          "res_enum_functions.h must be included by exactly 1 file, found ${n}: ${includers}")
endif()
if(NOT includers MATCHES "opt_oracle_respol\\.cpp$")
  message(FATAL_ERROR "the sole includer must be src/opt_oracle_respol.cpp, got ${includers}")
endif()

file(GLOB_RECURSE headers "${ROOT}/include/*.h" "${ROOT}/include/*.hpp")
foreach(h ${headers})
  file(READ "${h}" text)
  if(text MATCHES "#include[ \t]*<CGAL/" OR text MATCHES "#include[ \t]*<LEDA/")
    message(FATAL_ERROR "public header leaks CGAL/LEDA: ${h}")
  endif()
endforeach()

message(STATUS "source_purity: PASS (1 respol includer, no CGAL/LEDA in public headers)")
