# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/marble
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_RIGID_2D/marble
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(mecaMAILx_RIGID_2D_marble "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_RIGID_2D_marble" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_RIGID_2D/marble/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_RIGID_2D/marble/command.py" "/examples/mecaMAILx_RIGID_2D/marble")
  set_tests_properties(mecaMAILx_RIGID_2D_marble PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/marble/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/marble/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_mecaMAILx_RIGID_2D_marble "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_RIGID_2D_marble" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_RIGID_2D/marble/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_RIGID_2D/marble/command.py" "/examples/mecaMAILx_RIGID_2D/marble")
  set_tests_properties(save_mecaMAILx_RIGID_2D_marble PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/marble/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/marble/CMakeLists.txt;0;")
endif()
