# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/compression_temperature
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/compression_temperature
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Nn][Oo][Nn]_[Rr][Ee][Gg])$")
  add_test(demmefi_compression_temperature "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/run_test.sh" "/opt/homebrew/bin/python3.11" "demmefi_compression_temperature" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/compression_temperature/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/compression_temperature/command.py" "/examples/mecaMAILx_3D/demmefi/demmefi_compression_temperature")
  set_tests_properties(demmefi_compression_temperature PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;40;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/compression_temperature/CMakeLists.txt;5;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/compression_temperature/CMakeLists.txt;0;")
endif()
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ss][Aa][Vv][Ee]_[Rr][Ee][Gg])$")
  add_test(save_demmefi_compression_temperature "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/save_test.sh" "/opt/homebrew/bin/python3.11" "demmefi_compression_temperature" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/compression_temperature/gen_sample.py" "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/compression_temperature/command.py" "/examples/mecaMAILx_3D/demmefi/demmefi_compression_temperature")
  set_tests_properties(save_demmefi_compression_temperature PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native:" LABELS "all;quick" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/ci_scripts/make_test.cmake;59;add_test;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/compression_temperature/CMakeLists.txt;5;createTest;/Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D/demmefi/compression_temperature/CMakeLists.txt;0;")
endif()
