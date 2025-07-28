# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall2D/test_exploded_deformable_bricks
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall2D/test_exploded_deformable_bricks
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(Pre_wall2D_exploded_deformable_bricks "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "Pre_wall2D_exploded_deformable_bricks" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall2D/test_exploded_deformable_bricks/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/Pre/wall2D/test_exploded_deformable_bricks/command.py")
  set_tests_properties(Pre_wall2D_exploded_deformable_bricks PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall2D/test_exploded_deformable_bricks/CMakeLists.txt;6;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/Pre/wall2D/test_exploded_deformable_bricks/CMakeLists.txt;0;")
endif()
