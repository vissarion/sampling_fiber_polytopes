# Assert that D, CD and PD are NOT exported from the respol island.
#
# res_enum_functions.h reads these as globals and respol's res_enum_d.cpp:28-30
# defines them non-static. We declare them `static` in src/opt_oracle_respol.cpp
# so they get internal linkage and vanish from the symbol table entirely.

find_program(NM_EXE nm)
if(NOT NM_EXE)
  message(STATUS "nm not available; skipping globals check")
  return()
endif()

execute_process(COMMAND ${NM_EXE} ${LIB}
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "nm failed on ${LIB}: ${err}")
endif()

# Match a defined (B/D/G = bss/data/global-data) symbol named exactly D, CD or PD.
string(REGEX MATCHALL "\n[0-9a-fA-F]+ [BDGbdgSs] (D|CD|PD)\r?\n" bad "\n${out}\n")
list(LENGTH bad n_bad)
if(n_bad GREATER 0)
  message(FATAL_ERROR
          "respol globals leaked into the object file (${n_bad}): ${bad}\n"
          "Declare them `static` before #include <res_enum_functions.h>.")
endif()

message(STATUS "globals_not_exported: PASS (no D/CD/PD exported)")
