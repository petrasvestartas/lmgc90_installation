# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/RIGID_3D/CpF2F
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/RIGID_3D/CpF2F
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_cpf2f "/opt/homebrew/bin/python3.11" "/Users/petras/brg/2_code/lmgc90_installation/build_native/test/RIGID_3D/CpF2F/test.py")
set_tests_properties(python_cpf2f PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/RIGID_3D/CpF2F/CMakeLists.txt;11;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/RIGID_3D/CpF2F/CMakeLists.txt;0;")
