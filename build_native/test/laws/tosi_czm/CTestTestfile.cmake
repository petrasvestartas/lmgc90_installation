# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/laws/tosi_czm
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/laws/tosi_czm
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Ll][Oo][Ww])$")
  add_test(laws_tosi_czm "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/laws/tosi_czm/run_all.py" "--novisu")
  set_tests_properties(laws_tosi_czm PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/laws/tosi_czm/CMakeLists.txt;23;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/laws/tosi_czm/CMakeLists.txt;0;")
endif()
