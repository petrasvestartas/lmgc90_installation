# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/2Q4_PressurePostCracking
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/2Q4_PressurePostCracking
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(mecaMAILx_2D_2Q4_PressurePostCracking "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_2D_2Q4_PressurePostCracking" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/2Q4_PressurePostCracking/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/2Q4_PressurePostCracking/command.py" "/examples/mecaMAILx_2D/2Q4_PressurePostCracking")
  set_tests_properties(mecaMAILx_2D_2Q4_PressurePostCracking PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/2Q4_PressurePostCracking/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/2Q4_PressurePostCracking/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_mecaMAILx_2D_2Q4_PressurePostCracking "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_2D_2Q4_PressurePostCracking" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/2Q4_PressurePostCracking/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/2Q4_PressurePostCracking/command.py" "/examples/mecaMAILx_2D/2Q4_PressurePostCracking")
  set_tests_properties(save_mecaMAILx_2D_2Q4_PressurePostCracking PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/2Q4_PressurePostCracking/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/2Q4_PressurePostCracking/CMakeLists.txt;0;")
endif()
