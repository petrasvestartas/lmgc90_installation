# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/2H8_traction
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/2H8_traction
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(mecaMAILx_3D_2H8_traction "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_3D_2H8_traction" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/2H8_traction/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/2H8_traction/command.py")
  set_tests_properties(mecaMAILx_3D_2H8_traction PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;correct" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/2H8_traction/CMakeLists.txt;6;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/2H8_traction/CMakeLists.txt;0;")
endif()
