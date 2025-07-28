# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/triaxial_compression
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/triaxial_compression
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(demmefi_triaxial_compression "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "demmefi_triaxial_compression" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/triaxial_compression/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/triaxial_compression/command.py" "/examples/mecaMAILx_3D/demmefi/demmefi_triaxial_compression")
  set_tests_properties(demmefi_triaxial_compression PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/triaxial_compression/CMakeLists.txt;6;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/triaxial_compression/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_demmefi_triaxial_compression "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "demmefi_triaxial_compression" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/triaxial_compression/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/triaxial_compression/command.py" "/examples/mecaMAILx_3D/demmefi/demmefi_triaxial_compression")
  set_tests_properties(save_demmefi_triaxial_compression PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/triaxial_compression/CMakeLists.txt;6;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/triaxial_compression/CMakeLists.txt;0;")
endif()
