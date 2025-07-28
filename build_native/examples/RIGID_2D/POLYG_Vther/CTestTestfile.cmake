# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/POLYG_Vther
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/POLYG_Vther
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(RIGID_2D_POLYG_Vther "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "RIGID_2D_POLYG_Vther" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/POLYG_Vther/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/POLYG_Vther/command.py" "/examples/RIGID_2D/POLYG_Vther")
  set_tests_properties(RIGID_2D_POLYG_Vther PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/POLYG_Vther/CMakeLists.txt;4;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/POLYG_Vther/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_RIGID_2D_POLYG_Vther "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "RIGID_2D_POLYG_Vther" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/POLYG_Vther/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/POLYG_Vther/command.py" "/examples/RIGID_2D/POLYG_Vther")
  set_tests_properties(save_RIGID_2D_POLYG_Vther PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/POLYG_Vther/CMakeLists.txt;4;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/POLYG_Vther/CMakeLists.txt;0;")
endif()
