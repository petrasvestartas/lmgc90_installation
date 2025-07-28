# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_2D/ConfinedCompression
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_2D/ConfinedCompression
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(poroMAILx_2D_ConfinedCompression "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "poroMAILx_2D_ConfinedCompression" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_2D/ConfinedCompression/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_2D/ConfinedCompression/command.py" "/examples/poroMAILx_2D/ConfinedCompression")
  set_tests_properties(poroMAILx_2D_ConfinedCompression PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_2D/ConfinedCompression/CMakeLists.txt;11;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_2D/ConfinedCompression/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_poroMAILx_2D_ConfinedCompression "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "poroMAILx_2D_ConfinedCompression" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_2D/ConfinedCompression/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_2D/ConfinedCompression/command.py" "/examples/poroMAILx_2D/ConfinedCompression")
  set_tests_properties(save_poroMAILx_2D_ConfinedCompression PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_2D/ConfinedCompression/CMakeLists.txt;11;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_2D/ConfinedCompression/CMakeLists.txt;0;")
endif()
