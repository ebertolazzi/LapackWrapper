######################################################################
#
# CMakeLists.txt for SUPERLU
#
######################################################################

# Required version
cmake_minimum_required(VERSION 2.8.12)

# Project Version
project(SuperLU NONE)
set(VERSION_MAJOR "5")
set(VERSION_MINOR "2")
set(VERSION_BugFix "1")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

######################################################################
#
# Usual initialization stuff
#
######################################################################
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)    ## ????
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

IF ( APPLE OR UNIX )
  set(PROJECT_NAME_LIB_EXPORT libsuperlu.a)
ELSE()
  set(PROJECT_NAME_LIB_EXPORT libsuperlu.lib)
ENDIF()

enable_language(C)
#enable_language (Fortran)
set(NOFORTRAN TRUE)
set(SUPERLU_VERSION "${PROJECT_VERSION}")
set(SUPERLU_REV "${PROJECT_REV}")

# setup options
option(enable_matlabmex "Build the Matlab mex library" OFF)
option(enable_tests     "Build tests" OFF)
option(enable_doc       "Build doxygen documentation" OFF)
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)

# setup required compiler defines and options.
## add_definitions(-DPRNTlevel=0 -DAdd_)
## get_directory_property( DirDefs COMPILE_DEFINITIONS )
set(CMAKE_C_FLAGS         "-fPIC -DPRNTlevel=0 -DAdd_ ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "-O3" CACHE STRING "")

######################################################################
#
# Include directories
#
######################################################################

include_directories( superlu/SRC )

######################################################################
#
# Add subdirectories
#
######################################################################

add_subdirectory( superlu/SRC )


# file(WRITE "make.defs" "# can be exposed to users"
#  ${CMAKE_C_COMPILER}  )
configure_file( superlu/make.inc.in superlu/make.inc )
