# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/therMAILx_1D/Reseau_1D_with_Source
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/therMAILx_1D/Reseau_1D_with_Source
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(therMAILx_1D_Reseau1D_with_source "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "therMAILx_1D_Reseau1D_with_source" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/therMAILx_1D/Reseau_1D_with_Source/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/therMAILx_1D/Reseau_1D_with_Source/command.py")
  set_tests_properties(therMAILx_1D_Reseau1D_with_source PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/therMAILx_1D/Reseau_1D_with_Source/CMakeLists.txt;10;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/therMAILx_1D/Reseau_1D_with_Source/CMakeLists.txt;0;")
endif()
