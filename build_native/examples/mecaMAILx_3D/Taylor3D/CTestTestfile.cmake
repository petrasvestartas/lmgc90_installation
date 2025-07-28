# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/Taylor3D
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Taylor3D
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(mecaMAILx_3D_Taylor3D "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_3D_Taylor3D" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Taylor3D/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Taylor3D/command.py")
  set_tests_properties(mecaMAILx_3D_Taylor3D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;long" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/Taylor3D/CMakeLists.txt;10;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/Taylor3D/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(mecaMAILx_3D_Taylor3DExplicit "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_3D_Taylor3DExplicit" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Taylor3D/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Taylor3D/command_explicit.py")
  set_tests_properties(mecaMAILx_3D_Taylor3DExplicit PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;long" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/Taylor3D/CMakeLists.txt;16;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/Taylor3D/CMakeLists.txt;0;")
endif()
