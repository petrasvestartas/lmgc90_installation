# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/extruded_box_cylinders
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/extruded_box_cylinders
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(Pre_grains_extruded_box_cylinders "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "Pre_grains_extruded_box_cylinders" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/extruded_box_cylinders/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/prepro_grains/extruded_box_cylinders/command.py")
  set_tests_properties(Pre_grains_extruded_box_cylinders PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;slow" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/extruded_box_cylinders/CMakeLists.txt;6;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/prepro_grains/extruded_box_cylinders/CMakeLists.txt;0;")
endif()
