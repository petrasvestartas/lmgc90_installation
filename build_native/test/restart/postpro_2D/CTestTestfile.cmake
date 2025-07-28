# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/restart/postpro_2D
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/restart/postpro_2D
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_restart_postpro_2d "/opt/homebrew/bin/python3.11" "test.py" "--novisu")
set_tests_properties(python_restart_postpro_2d PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/restart/postpro_2D/CMakeLists.txt;13;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/restart/postpro_2D/CMakeLists.txt;0;")
