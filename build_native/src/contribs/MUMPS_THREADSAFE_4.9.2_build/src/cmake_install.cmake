# Install script for directory: /Users/petras/brg/2_code/lmgc90_installation/src/contribs/MUMPS_THREADSAFE_4.9.2/src

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/lib/libmumps_common.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmumps_common.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmumps_common.dylib")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -id "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libmumps_common.dylib"
      -change "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/lib/libmpiseq.dylib" "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libmpiseq.dylib"
      -change "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/lib/libpord.dylib" "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libpord.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmumps_common.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmumps_common.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmumps_common.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/src/CMakeFiles/mumps_common.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/lib/libdmumps.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdmumps.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdmumps.dylib")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -id "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libdmumps.dylib"
      -change "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/lib/libmpiseq.dylib" "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libmpiseq.dylib"
      -change "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/lib/libmumps_common.dylib" "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libmumps_common.dylib"
      -change "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/lib/libpord.dylib" "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libpord.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdmumps.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdmumps.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdmumps.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/src/CMakeFiles/dmumps.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build/src/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
