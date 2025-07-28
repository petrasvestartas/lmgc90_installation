# Install script for directory: /Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/petras/brg/2_code/lmgc90_installation/build_native")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
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

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MatLib4-2021-build/lib/libmatlib.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmatlib.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmatlib.dylib")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -id "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libmatlib.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmatlib.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmatlib.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmatlib.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/include/MatLib.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/include/c_matlib.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/include/matlib_macros.h"
    "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MatLib4-2021-build/include/matlib_config.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Chronometer.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Cloneable.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/ConstantFunction.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Copiable.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Exceptions.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/FileException.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Function.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Function2.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Property.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Rotation.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Rotation2D.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/Rotation3D.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/ShortArray.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/ShortMatrix.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/ShortSqrMatrix.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/StringMap.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/SyntaxError.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/TabulatedFunction.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/data/TabulatedFunction2.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/math/MathUtils.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/math/Vector1D.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/math/Vector2D.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/math/Vector3D.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/opti/OptiFunction.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/opti/OptiMethod.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/opti/OptiProblem.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/matl/ConstitutiveModel.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/matl/CriterionDictionary.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/matl/MaterialCriterion.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/matl/MaterialModel.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/matl/MaterialProperties.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/matl/MaterialState.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MatLib4-LMGC-2021-1/src/matl/ModelDictionary.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MatLib4-2021-build/appl/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MatLib4-2021-build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MatLib4-2021-build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
