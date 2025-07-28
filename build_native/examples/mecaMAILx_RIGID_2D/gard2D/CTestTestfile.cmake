# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/gard2D
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_RIGID_2D/gard2D
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(mecaMAILx_RIGID_2D_gard2D "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_RIGID_2D_gard2D" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_RIGID_2D/gard2D/command.py")
  set_tests_properties(mecaMAILx_RIGID_2D_gard2D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/gard2D/CMakeLists.txt;11;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_RIGID_2D/gard2D/CMakeLists.txt;0;")
endif()
