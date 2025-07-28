# Install script for directory: /Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP

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
    set(CMAKE_INSTALL_CONFIG_NAME "")
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libClipper2.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper2" TYPE FILE FILES
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.version.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.core.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.engine.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.export.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.minkowski.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.offset.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.rectclip.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libClipper2Z.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2Z.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2Z.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2Z.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper2" TYPE FILE FILES
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.version.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.core.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.engine.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.export.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.minkowski.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.offset.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.rectclip.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES
    "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/Clipper2.pc"
    "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/Clipper2Z.pc"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libClipper2.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper2" TYPE FILE FILES
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.version.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.core.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.engine.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.export.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.minkowski.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.offset.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.rectclip.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/lib/libClipper2Z.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2Z.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2Z.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libClipper2Z.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper2" TYPE FILE FILES
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.version.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.core.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.engine.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.export.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.minkowski.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.offset.h"
    "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP/Clipper2Lib/include/clipper2/clipper.rectclip.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES
    "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/Clipper2Config.cmake"
    "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/Clipper2ConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets.cmake"
         "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets-noconfig.cmake")
  endif()
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/install_local_manifest.txt"
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
  file(WRITE "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
