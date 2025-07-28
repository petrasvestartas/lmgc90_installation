# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/Drum_THERMORIGID
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/Drum_THERMORIGID
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(RIGID_2D_Drum_THERMORIGID "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "RIGID_2D_Drum_THERMORIGID" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/Drum_THERMORIGID/gen_drum.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/RIGID_2D/Drum_THERMORIGID/command.py")
  set_tests_properties(RIGID_2D_Drum_THERMORIGID PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;slow" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/Drum_THERMORIGID/CMakeLists.txt;6;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/RIGID_2D/Drum_THERMORIGID/CMakeLists.txt;0;")
endif()
