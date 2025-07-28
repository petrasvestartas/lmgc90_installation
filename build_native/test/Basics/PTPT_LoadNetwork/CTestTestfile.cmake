# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/Basics/PTPT_LoadNetwork
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/PTPT_LoadNetwork
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_basic_PTPT_LoadNetwork_3D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/PTPT_LoadNetwork/test.py" "3")
set_tests_properties(python_basic_PTPT_LoadNetwork_3D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/Basics/PTPT_LoadNetwork/CMakeLists.txt;15;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/Basics/PTPT_LoadNetwork/CMakeLists.txt;0;")
add_test(python_basic_PTPT_LoadNetwork_2D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/PTPT_LoadNetwork/test.py" "2")
set_tests_properties(python_basic_PTPT_LoadNetwork_2D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/Basics/PTPT_LoadNetwork/CMakeLists.txt;27;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/Basics/PTPT_LoadNetwork/CMakeLists.txt;0;")
