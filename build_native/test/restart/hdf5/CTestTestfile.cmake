# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/restart/hdf5
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_restart_hdf5_2D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/restart/hdf5/run_test.py" "2")
set_tests_properties(python_restart_hdf5_2D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/CMakeLists.txt;12;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/CMakeLists.txt;0;")
add_test(python_restart_hdf5_3D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/restart/hdf5/run_test.py" "3")
set_tests_properties(python_restart_hdf5_3D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/CMakeLists.txt;15;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/CMakeLists.txt;0;")
subdirs("version")
