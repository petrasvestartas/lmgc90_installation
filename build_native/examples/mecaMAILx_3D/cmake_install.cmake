# Install script for directory: /Users/petras/brg/2_code/lmgc90_installation/examples/mecaMAILx_3D

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/opt/homebrew/lib")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/One_BAR/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Poutre_BAR/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/1J_compression+shear/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/1J_traction/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/jointFE_Arche/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/jointFE_traction/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/jointFE_traction+shear/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/jointFE_compression+shear/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/demmefi/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/SHB/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubes/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubes_H8/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/mtl_poutre_H8/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/mtl_poutre_H8_gd/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/tet4_tet4/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/discreteFE_suspendedblock/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/polyd_arch_dila/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/2H8_traction/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Autocontact/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/Taylor3D/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/cubetet4_cubetet4/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/umat_and_fields/cmake_install.cmake")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/petras/brg/2_code/lmgc90_installation/build_native/examples/mecaMAILx_3D/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
