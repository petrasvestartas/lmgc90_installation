# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/Basics/recup_rloc
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/recup_rloc
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_basic_recuprloc_2D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/recup_rloc/test_2d.py" "&&")
set_tests_properties(python_basic_recuprloc_2D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/Basics/recup_rloc/CMakeLists.txt;14;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/Basics/recup_rloc/CMakeLists.txt;0;")
add_test(python_basic_recuprloc_3D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/recup_rloc/test_3d.py" "&&")
set_tests_properties(python_basic_recuprloc_3D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/Basics/recup_rloc/CMakeLists.txt;24;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/Basics/recup_rloc/CMakeLists.txt;0;")
