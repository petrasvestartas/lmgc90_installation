# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MUMPS_THREADSAFE_4.9.2")
  file(MAKE_DIRECTORY "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/MUMPS_THREADSAFE_4.9.2")
endif()
file(MAKE_DIRECTORY
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/MUMPS_THREADSAFE_4.9.2_build"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/mumps_project-prefix"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/mumps_project-prefix/tmp"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/mumps_project-prefix/src/mumps_project-stamp"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/mumps_project-prefix/src"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/mumps_project-prefix/src/mumps_project-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/mumps_project-prefix/src/mumps_project-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/mumps_project-prefix/src/mumps_project-stamp${cfgdir}") # cfgdir has leading slash
endif()
