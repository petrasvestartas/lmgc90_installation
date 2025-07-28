# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/1024DK_BoxJC
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/1024DK_BoxJC
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(RIGID_2D_1024DK_BoxJC "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "RIGID_2D_1024DK_BoxJC" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/1024DK_BoxJC/command.py")
  set_tests_properties(RIGID_2D_1024DK_BoxJC PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;slow" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/1024DK_BoxJC/CMakeLists.txt;10;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/1024DK_BoxJC/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(RIGID_2D_1024DK_BoxJC_post "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "RIGID_2D_1024DK_BoxJC_post" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/1024DK_BoxJC/command_post.py")
  set_tests_properties(RIGID_2D_1024DK_BoxJC_post PROPERTIES  DEPENDS "RIGID_2D_1024DK_BoxJC" ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;slow" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/1024DK_BoxJC/CMakeLists.txt;11;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/1024DK_BoxJC/CMakeLists.txt;0;")
endif()
