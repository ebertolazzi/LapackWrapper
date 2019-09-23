######################################################################
#
# CMakeLists.txt for SUPERLU
#
######################################################################

# Required version
cmake_minimum_required(VERSION 2.8.12)
SET( CMAKE_VERBOSE_MAKEFILE OFF )
SET( CMAKE_INSTALL_MESSAGE NEVER )

# Project Version
PROJECT( SuperLU C CXX )
set(VERSION_MAJOR "5")
set(VERSION_MINOR "2")
set(VERSION_BugFix "1")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

######################################################################
#
# Usual initialization stuff
#
######################################################################
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)    ## ????
SET(CMAKE_INSTALL_PREFIX     "${CMAKE_CURRENT_SOURCE_DIR}/lib")
SET(CMAKE_INSTALL_LIBDIR     "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_BINDIR     "${CMAKE_INSTALL_PREFIX}/bin")
SET(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_PREFIX}/include")

IF ( APPLE OR UNIX )
  SET(PROJECT_NAME_LIB_EXPORT libsuperlu.a)
ELSE()
  ADD_DEFINITIONS( -D_CRT_SECURE_NO_WARNINGS )
  SET(PROJECT_NAME_LIB_EXPORT libsuperlu.lib)
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

INCLUDE( ../../CMakeLists-cflags.txt )
## add_definitions(-DPRNTlevel=0 -DAdd_)
## get_directory_property( DirDefs COMPILE_DEFINITIONS )
ADD_DEFINITIONS( -DPRNTlevel=0 -DAdd_ )
# suppress all warnings
SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-everything" )


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
