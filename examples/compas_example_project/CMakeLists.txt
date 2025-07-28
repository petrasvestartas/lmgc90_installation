cmake_minimum_required(VERSION 3.10)
project(OCCT_Example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the installed OpenCASCADE library paths (from fork's install directory)
set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../install")
set(OpenCASCADE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../install/lib/cmake/opencascade")

# Find OpenCASCADE package
find_package(OpenCASCADE REQUIRED)

# Include directories
include_directories(
    ${OpenCASCADE_INCLUDE_DIR}
)

# Add executable
add_executable(write_step write_step.cpp)

# Link against OpenCASCADE libraries
target_link_libraries(write_step
    # Base modules
    TKernel TKMath 
    # Geometry modules
    TKG2d TKG3d TKGeomBase TKGeomAlgo TKBRep 
    # Topology modules
    TKTopAlgo TKShHealing 
    # Boolean operations
    TKBool TKBO 
    # Modeling modules
    TKFillet TKOffset TKFeat TKPrim TKMesh 
    # STEP modules
    TKCDF TKLCAF TKCAF TKV3d
    TKService 
    TKXSBase TKDE 
    # Data exchange modules
    TKXCAF TKDESTEP
    # Binary and XML persistence
    TKBin TKXml TKStd TKStdL TKBinL TKXmlL
    # XCAF Binary persistence
    TKBinXCAF TKXmlXCAF
    # Other modules
    TKTObj TKBinTObj TKXmlTObj
)

# Print some information about the OpenCASCADE libraries found
message(STATUS "OpenCASCADE_INCLUDE_DIR: ${OpenCASCADE_INCLUDE_DIR}")
message(STATUS "OpenCASCADE_LIBRARIES: ${OpenCASCADE_LIBRARIES}")
