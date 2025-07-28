# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/src/contribs/ann-1.1.2/test
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/ann-1.1.2/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ann_python "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/ann-1.1.2/test/test_pair.py")
set_tests_properties(ann_python PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/ann-1.1.2/test/CMakeLists.txt;9;add_test;/Users/petras/brg/2_code/lmgc90_installation/src/contribs/ann-1.1.2/test/CMakeLists.txt;0;")
