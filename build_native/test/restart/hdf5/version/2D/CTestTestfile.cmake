# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/version/2D
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/restart/hdf5/version/2D
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_restart_hfile_2D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/restart/hdf5/version/2D/read_file.py" "&&")
set_tests_properties(python_restart_hfile_2D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/version/2D/CMakeLists.txt;16;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/version/2D/CMakeLists.txt;0;")
add_test(python_read_pre_hfile_2D "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/restart/hdf5/version/2D/read_pre.py" "&&")
set_tests_properties(python_read_pre_hfile_2D PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/version/2D/CMakeLists.txt;27;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/restart/hdf5/version/2D/CMakeLists.txt;0;")
