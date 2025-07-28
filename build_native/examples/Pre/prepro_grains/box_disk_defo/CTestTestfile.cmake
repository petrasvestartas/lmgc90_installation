# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/box_disk_defo
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/box_disk_defo
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(Pre_grains_box_disk_defo "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "Pre_grains_box_disk_defo" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/box_disk_defo/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/box_disk_defo/command.py" "/examples/Pre/prepro_grains/box_disk_defo")
  set_tests_properties(Pre_grains_box_disk_defo PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/box_disk_defo/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/box_disk_defo/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_Pre_grains_box_disk_defo "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "Pre_grains_box_disk_defo" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/box_disk_defo/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/box_disk_defo/command.py" "/examples/Pre/prepro_grains/box_disk_defo")
  set_tests_properties(save_Pre_grains_box_disk_defo PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/box_disk_defo/CMakeLists.txt;3;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/box_disk_defo/CMakeLists.txt;0;")
endif()
