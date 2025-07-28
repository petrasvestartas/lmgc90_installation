# CMake generated Testfile for 
# Source directory: /Users/petras/brg/2_code/lmgc90_installation/test/Basics/contact_3D/stono/2polys
# Build directory: /Users/petras/brg/2_code/lmgc90_installation/build_native/test/Basics/contact_3D/stono/2polys
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(python_basic_stono_2polys "/opt/homebrew/bin/python3.11" "run_all.py")
set_tests_properties(python_basic_stono_2polys PROPERTIES  ENVIRONMENT "PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build_native::;OMP_SCHEDULE=STATIC;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1" _BACKTRACE_TRIPLES "/Users/petras/brg/2_code/lmgc90_installation/test/Basics/contact_3D/stono/2polys/CMakeLists.txt;18;add_test;/Users/petras/brg/2_code/lmgc90_installation/test/Basics/contact_3D/stono/2polys/CMakeLists.txt;0;")
