# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/cubes
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubes
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(mecaMAILx_3D_cubes "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_3D_cubes" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubes/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubes/command.py" "/examples/mecaMAILx_3D/cubes")
  set_tests_properties(mecaMAILx_3D_cubes PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/cubes/CMakeLists.txt;7;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/cubes/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_mecaMAILx_3D_cubes "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_3D_cubes" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubes/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubes/command.py" "/examples/mecaMAILx_3D/cubes")
  set_tests_properties(save_mecaMAILx_3D_cubes PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/cubes/CMakeLists.txt;7;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/cubes/CMakeLists.txt;0;")
endif()
