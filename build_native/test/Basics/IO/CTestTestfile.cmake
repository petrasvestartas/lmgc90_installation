# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/Basics/IO
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/IO
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_basic_io "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/IO/run_all.py" "--with-hdf5")
set_tests_properties(python_basic_io PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/Basics/IO/CMakeLists.txt;25;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/Basics/IO/CMakeLists.txt;0;")
subdirs("hdf5")
