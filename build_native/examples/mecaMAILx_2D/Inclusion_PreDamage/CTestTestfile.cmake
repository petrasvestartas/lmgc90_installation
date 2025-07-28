# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/Inclusion_PreDamage
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/Inclusion_PreDamage
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Jj][Uu][Ss][Tt]_[Rr][Uu][Nn])$")
  add_test(mecaMAILx_2D_Inclusion_PreDamage "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/just_run_it.sh" "/opt/homebrew/bin/python3.11" "mecaMAILx_2D_Inclusion_PreDamage" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/Inclusion_PreDamage/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_2D/Inclusion_PreDamage/command.py")
  set_tests_properties(mecaMAILx_2D_Inclusion_PreDamage PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "long" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/Inclusion_PreDamage/CMakeLists.txt;8;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_2D/Inclusion_PreDamage/CMakeLists.txt;0;")
endif()
