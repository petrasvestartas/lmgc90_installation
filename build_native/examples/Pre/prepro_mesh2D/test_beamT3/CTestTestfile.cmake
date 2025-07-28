# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_mesh2D/test_beamT3
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_mesh2D/test_beamT3
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(Pre_mesh2D_test_beamT3 "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "Pre_mesh2D_test_beamT3" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_mesh2D/test_beamT3/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_mesh2D/test_beamT3/command.py" "/examples/Pre/prepro_mesh2D/test_beamT3")
  set_tests_properties(Pre_mesh2D_test_beamT3 PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_mesh2D/test_beamT3/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_mesh2D/test_beamT3/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_Pre_mesh2D_test_beamT3 "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "Pre_mesh2D_test_beamT3" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_mesh2D/test_beamT3/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_mesh2D/test_beamT3/command.py" "/examples/Pre/prepro_mesh2D/test_beamT3")
  set_tests_properties(save_Pre_mesh2D_test_beamT3 PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_mesh2D/test_beamT3/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_mesh2D/test_beamT3/CMakeLists.txt;0;")
endif()
