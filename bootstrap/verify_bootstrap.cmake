# LMGC90 Bootstrap Verification Script
# This script verifies that all dependencies were installed correctly

set(BOOTSTRAP_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/deps")

message(STATUS "🔍 Verifying LMGC90 Bootstrap Installation...")
message(STATUS "Bootstrap prefix: ${BOOTSTRAP_PREFIX}")

set(VERIFICATION_FAILED FALSE)

# Function to check if a file exists and report
function(check_file_exists filepath description)
    if(EXISTS "${filepath}")
        message(STATUS "✅ ${description}: ${filepath}")
    else()
        message(STATUS "❌ ${description}: MISSING - ${filepath}")
        set(VERIFICATION_FAILED TRUE PARENT_SCOPE)
    endif()
endfunction()

# Function to check if a program works
function(check_program_works program description)
    execute_process(
        COMMAND "${program}" --version
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_QUIET
    )
    if(result EQUAL 0)
        string(REGEX REPLACE "\n.*" "" first_line "${output}")
        message(STATUS "✅ ${description}: ${first_line}")
    else()
        message(STATUS "❌ ${description}: FAILED - ${program}")
        set(VERIFICATION_FAILED TRUE PARENT_SCOPE)
    endif()
endfunction()

message(STATUS "")
message(STATUS "📋 Checking GCC Toolchain...")
check_file_exists("${BOOTSTRAP_PREFIX}/gcc/bin/gcc" "GCC Compiler")
check_file_exists("${BOOTSTRAP_PREFIX}/gcc/bin/g++" "G++ Compiler")
check_file_exists("${BOOTSTRAP_PREFIX}/gcc/bin/gfortran" "GFortran Compiler")

if(EXISTS "${BOOTSTRAP_PREFIX}/gcc/bin/gcc")
    check_program_works("${BOOTSTRAP_PREFIX}/gcc/bin/gcc" "GCC Version")
endif()
if(EXISTS "${BOOTSTRAP_PREFIX}/gcc/bin/g++")
    check_program_works("${BOOTSTRAP_PREFIX}/gcc/bin/g++" "G++ Version")
endif()
if(EXISTS "${BOOTSTRAP_PREFIX}/gcc/bin/gfortran")
    check_program_works("${BOOTSTRAP_PREFIX}/gcc/bin/gfortran" "GFortran Version")
endif()

message(STATUS "")
message(STATUS "📋 Checking HDF5 Installation...")
check_file_exists("${BOOTSTRAP_PREFIX}/lib/libhdf5.dylib" "HDF5 Library (macOS)")
check_file_exists("${BOOTSTRAP_PREFIX}/lib/libhdf5.so" "HDF5 Library (Linux)")
check_file_exists("${BOOTSTRAP_PREFIX}/lib/libhdf5_fortran.dylib" "HDF5 Fortran Library (macOS)")
check_file_exists("${BOOTSTRAP_PREFIX}/lib/libhdf5_fortran.so" "HDF5 Fortran Library (Linux)")
check_file_exists("${BOOTSTRAP_PREFIX}/include/hdf5.h" "HDF5 Headers")
check_file_exists("${BOOTSTRAP_PREFIX}/bin/h5dump" "HDF5 Tools")

if(EXISTS "${BOOTSTRAP_PREFIX}/bin/h5dump")
    check_program_works("${BOOTSTRAP_PREFIX}/bin/h5dump" "HDF5 Version")
endif()

message(STATUS "")
message(STATUS "📋 Checking OpenBLAS Installation...")
check_file_exists("${BOOTSTRAP_PREFIX}/lib/libopenblas.dylib" "OpenBLAS Library (macOS)")
check_file_exists("${BOOTSTRAP_PREFIX}/lib/libopenblas.so" "OpenBLAS Library (Linux)")
check_file_exists("${BOOTSTRAP_PREFIX}/include/openblas_config.h" "OpenBLAS Headers")

message(STATUS "")
message(STATUS "📋 Checking SWIG Installation...")
check_file_exists("${BOOTSTRAP_PREFIX}/bin/swig" "SWIG Executable")

if(EXISTS "${BOOTSTRAP_PREFIX}/bin/swig")
    check_program_works("${BOOTSTRAP_PREFIX}/bin/swig" "SWIG Version")
    
    # Check SWIG version specifically (must be 3.0.12)
    execute_process(
        COMMAND "${BOOTSTRAP_PREFIX}/bin/swig" -version
        OUTPUT_VARIABLE swig_output
        ERROR_QUIET
    )
    if(swig_output MATCHES "SWIG Version 3\\.0\\.12")
        message(STATUS "✅ SWIG Version: Correct (3.0.12)")
    else()
        message(STATUS "❌ SWIG Version: Wrong version (expected 3.0.12)")
        set(VERIFICATION_FAILED TRUE)
    endif()
endif()

message(STATUS "")
message(STATUS "📋 Checking Python Dependencies...")
find_program(PYTHON_EXECUTABLE python3)
if(PYTHON_EXECUTABLE)
    message(STATUS "✅ Python3: ${PYTHON_EXECUTABLE}")
    
    # Check Python packages
    set(PYTHON_PACKAGES "numpy" "scipy" "matplotlib" "vtk" "compas")
    foreach(package ${PYTHON_PACKAGES})
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c "import ${package}; print('${package}:', ${package}.__version__)"
            RESULT_VARIABLE result
            OUTPUT_VARIABLE output
            ERROR_QUIET
        )
        if(result EQUAL 0)
            string(STRIP "${output}" output)
            message(STATUS "✅ Python Package: ${output}")
        else()
            message(STATUS "❌ Python Package: ${package} - MISSING")
            set(VERIFICATION_FAILED TRUE)
        endif()
    endforeach()
else()
    message(STATUS "❌ Python3: NOT FOUND")
    set(VERIFICATION_FAILED TRUE)
endif()

message(STATUS "")
message(STATUS "📋 Checking Environment Scripts...")
check_file_exists("${CMAKE_CURRENT_BINARY_DIR}/activate_bootstrap.sh" "Activation Script")
check_file_exists("${CMAKE_CURRENT_BINARY_DIR}/lmgc90_toolchain.cmake" "CMake Toolchain")

message(STATUS "")
if(VERIFICATION_FAILED)
    message(STATUS "❌ Bootstrap Verification FAILED!")
    message(STATUS "Some dependencies are missing or not working correctly.")
    message(STATUS "Please check the build logs and try running the bootstrap again.")
    message(FATAL_ERROR "Bootstrap verification failed")
else()
    message(STATUS "🎉 Bootstrap Verification SUCCESSFUL!")
    message(STATUS "All dependencies are properly installed and functional.")
    message(STATUS "")
    message(STATUS "Next steps:")
    message(STATUS "1. Source environment: source ${CMAKE_CURRENT_BINARY_DIR}/activate_bootstrap.sh")
    message(STATUS "2. Configure LMGC90: cmake -DCMAKE_TOOLCHAIN_FILE=${CMAKE_CURRENT_BINARY_DIR}/lmgc90_toolchain.cmake ..")
    message(STATUS "3. Build LMGC90: make -j8")
endif()
