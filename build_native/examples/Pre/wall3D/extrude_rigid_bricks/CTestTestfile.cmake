# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall3D/extrude_rigid_bricks
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall3D/extrude_rigid_bricks
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(Pre_wall3D_extrude_rigid_bricks "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "Pre_wall3D_extrude_rigid_bricks" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall3D/extrude_rigid_bricks/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall3D/extrude_rigid_bricks/command.py" "/examples/Pre/wall3D/extrude_rigid_bricks")
  set_tests_properties(Pre_wall3D_extrude_rigid_bricks PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall3D/extrude_rigid_bricks/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall3D/extrude_rigid_bricks/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_Pre_wall3D_extrude_rigid_bricks "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "Pre_wall3D_extrude_rigid_bricks" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall3D/extrude_rigid_bricks/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall3D/extrude_rigid_bricks/command.py" "/examples/Pre/wall3D/extrude_rigid_bricks")
  set_tests_properties(save_Pre_wall3D_extrude_rigid_bricks PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall3D/extrude_rigid_bricks/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall3D/extrude_rigid_bricks/CMakeLists.txt;0;")
endif()
