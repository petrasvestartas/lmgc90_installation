# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/laws/elastic_repell
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/laws/elastic_repell
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(laws_elastic_repell "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/laws/elastic_repell/run_all.py" "--novisu")
set_tests_properties(laws_elastic_repell PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/laws/elastic_repell/CMakeLists.txt;15;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/laws/elastic_repell/CMakeLists.txt;0;")
