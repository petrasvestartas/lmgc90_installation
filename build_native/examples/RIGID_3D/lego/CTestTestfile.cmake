# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_3D/lego
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_3D/lego
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(RIGID_3D_lego "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "RIGID_3D_lego" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_3D/lego/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_3D/lego/command.py" "/examples/RIGID_3D/lego")
  set_tests_properties(RIGID_3D_lego PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_3D/lego/CMakeLists.txt;7;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_3D/lego/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_RIGID_3D_lego "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "RIGID_3D_lego" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_3D/lego/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_3D/lego/command.py" "/examples/RIGID_3D/lego")
  set_tests_properties(save_RIGID_3D_lego PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_3D/lego/CMakeLists.txt;7;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_3D/lego/CMakeLists.txt;0;")
endif()
