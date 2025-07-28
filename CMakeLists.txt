CMAKE_MINIMUM_REQUIRED(VERSION 3.16.0)
# policy management
if(POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif(POLICY CMP0053)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif(POLICY CMP0077)
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif(POLICY CMP0078)
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif(POLICY CMP0086)
if(POLICY CMP0148)
  cmake_policy(SET CMP0148 NEW)
endif(POLICY CMP0148)

# path to homemade/taken macros/finds
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# Force out-of-source build
include(OutOfSourcesBuild)

project(LMGC90_dev Fortran C CXX)

SET(LMGC90_MAJOR_VERSION 3)
SET(LMGC90_MINOR_VERSION 0)
SET(LMGC90_PATCH_VERSION 0)

MESSAGE(STATUS "You try to build LMGC90_dev")

# Set they current year to help set version
set(YEAR 2025)

# Set release candidate number (used in sphinx doc)
set(RC 1)
if( RC )
  set(RCWITHDOT ".rc${RC}")
endif( RC )

include(CTest)

# maybe a little ugly...
if( EXISTS ${VENV_PATH} )
  # convert to PATH
  file(TO_CMAKE_PATH ${VENV_PATH} VENV_PATH)
  # check if not first
  if( DEFINED VENV_LAST )
  # check if value changed from last time...
    if( NOT (${VENV_PATH} STREQUAL ${VENV_LAST}) )
      unset(_Python_DEVELOPMENT_EMBED_SIGNATURE CACHE)
      unset(_Python_DEVELOPMENT_MODULE_SIGNATURE CACHE)
      unset(_Python_INTERPRETER_SIGNATURE CACHE)
      unset(_Python_NUMPY_SIGNATURE CACHE)
      set(CMAKE_INSTALL_PREFIX ${VENV_PATH})
      set(VENV_LAST ${VENV_PATH} CACHE PATH "Last venv path used")
    endif()
  else()
    set(CMAKE_INSTALL_PREFIX ${VENV_PATH})
    set(VENV_LAST ${VENV_PATH} CACHE PATH "Last venv path used")
  endif()
endif()

# Not that good... all other options are in src :s
option(NO_TEST "Deactivate any kind of test generation by CMake" OFF)

add_subdirectory(src)
if( NOT ${NO_TEST} )
  add_subdirectory(test)
  add_subdirectory(examples)
endif( NOT ${NO_TEST} )

