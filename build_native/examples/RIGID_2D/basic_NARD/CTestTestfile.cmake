# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/basic_NARD
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/basic_NARD
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(RIGID_2D_basic_NARD "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "RIGID_2D_basic_NARD" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/basic_NARD/command.py")
  set_tests_properties(RIGID_2D_basic_NARD PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/basic_NARD/CMakeLists.txt;10;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/basic_NARD/CMakeLists.txt;0;")
endif()
