# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP")
  file(MAKE_DIRECTORY "/Users/petras/brg/2_code/lmgc90_installation/src/contribs/clipper2-1.4.0/CPP")
endif()
file(MAKE_DIRECTORY
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_build"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_project-prefix"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_project-prefix/tmp"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_project-prefix/src/clipper_project-stamp"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_project-prefix/src"
  "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_project-prefix/src/clipper_project-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_project-prefix/src/clipper_project-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/petras/brg/2_code/lmgc90_installation/build_native/src/contribs/clipper_project-prefix/src/clipper_project-stamp${cfgdir}") # cfgdir has leading slash
endif()
