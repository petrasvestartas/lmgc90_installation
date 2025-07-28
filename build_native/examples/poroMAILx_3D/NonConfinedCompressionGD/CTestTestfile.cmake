# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_3D/NonConfinedCompressionGD
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_3D/NonConfinedCompressionGD
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(poroMAILx_3D_NonConfinedCompressionGD "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "poroMAILx_3D_NonConfinedCompressionGD" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_3D/NonConfinedCompressionGD/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/poroMAILx_3D/NonConfinedCompressionGD/command.py")
  set_tests_properties(poroMAILx_3D_NonConfinedCompressionGD PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;long" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_3D/NonConfinedCompressionGD/CMakeLists.txt;13;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/poroMAILx_3D/NonConfinedCompressionGD/CMakeLists.txt;0;")
endif()
