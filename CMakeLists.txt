############################################################################
#                                                                          #
#  LapackWrapper -- modern, self-contained CMake build.                    #
#                                                                          #
#  No dependency on external cmake_utils fragments.                         #
#  Structure follows the modern Splines CMake style: explicit targets,      #
#  FetchContent/local dependency handling, target-scoped compile options.   #
#                                                                          #
############################################################################

cmake_minimum_required( VERSION 3.25 )

# --------------------------------------------------------------------------
# Apple FAT/universal build setup
# --------------------------------------------------------------------------
# On Apple platforms a single library can contain slices for several CPU
# architectures, for example arm64 and x86_64.  This is controlled by
# CMAKE_OSX_ARCHITECTURES.
#
# The variable must be set before project(), otherwise CMake may initialize
# compilers and targets before the requested architectures are known.
#
# Default behavior:
#   - ON  on Apple hosts or Apple cross-compilation targets.
#   - OFF on Linux/Windows and other platforms, where there is no equivalent
#     single-file universal library format.
#
# User overrides:
#   -DUTILS_BUILD_FAT=OFF
#   -DUTILS_FAT_ARCHITECTURES=arm64
#   -DCMAKE_OSX_ARCHITECTURES=arm64\;x86_64
set(UTILS_BUILD_FAT_DEFAULT OFF)
if(CMAKE_HOST_APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin|iOS|tvOS|watchOS|visionOS")
  set(UTILS_BUILD_FAT_DEFAULT ON)
endif()
option(UTILS_BUILD_FAT "Build Apple universal/fat binaries when supported" ${UTILS_BUILD_FAT_DEFAULT})
set(UTILS_FAT_ARCHITECTURES "arm64;x86_64" CACHE STRING "Apple architectures to build when UTILS_BUILD_FAT is ON")

if(UTILS_BUILD_FAT AND NOT CMAKE_OSX_ARCHITECTURES)
  set(
    CMAKE_OSX_ARCHITECTURES
    "${UTILS_FAT_ARCHITECTURES}"
    CACHE STRING "Build architectures for Apple platforms"
    FORCE
  )
endif()

# ============================================================================
# Project version discovery
# ============================================================================
# The public project version is taken from the most recent Git tag when the
# source tree is a Git checkout.  This keeps archives, install metadata, CPack
# packages, and generated CMake package files aligned with the repository tag.
#
# If Git is not available, or if the checkout has no usable tag, the build falls
# back to 0.0.0.  The fallback keeps configuration deterministic for source
# snapshots, vendored copies, and CI jobs that build from an exported tree.
#
# Tags that start with a leading 'v' are normalized from forms such as v1.2.3
# to the CMake-compatible version string 1.2.3 before project() is called.
# ============================================================================
find_package( Git QUIET )
if ( Git_FOUND )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if ( GIT_VERSION MATCHES "^v" )
    string( SUBSTRING ${GIT_VERSION} 1 -1 GIT_VERSION )
  endif()
endif()
if ( NOT GIT_VERSION )
  set( GIT_VERSION "0.0.0" )
endif()

# project() enables the requested languages, initializes compiler detection,
# and makes PROJECT_VERSION/PROJECT_DESCRIPTION available to later install,
# target, and packaging rules.  The Apple FAT configuration above must precede
# this call so compiler checks are performed for the requested architecture set.
project(
  LapackWrapper
  VERSION ${GIT_VERSION}
  DESCRIPTION "A C++ wrapper for BLAS and LAPACK"
  HOMEPAGE_URL "https://ebertolazzi.github.io/LapackWrapper/"
  LANGUAGES CXX C
)

# ============================================================================
# Platform suffix used in real library artifact names
# ============================================================================
# The build emits platform-qualified library names, for example:
#   liblapack_wrapper_osx_static.a
#   liblapack_wrapper_linux_static.a
#   liblapack_wrapper_win64.dll / .lib, depending on the generator
#
# Unsuffixed names are created later as compatibility aliases.  The suffix makes
# side-by-side local builds easier to inspect and avoids overwriting artifacts
# produced for another platform or ABI.
#
# The logic distinguishes MinGW from native Windows generators because their
# file naming and import-library conventions differ in practice.  The pointer
# size test is used only to select win32/win64 or mingw32/mingw64 naming.
# ============================================================================
if ( APPLE )
  set( LAPACK_WRAPPER_PLATFORM_SUFFIX osx )
elseif ( WIN32 )
  if ( MINGW )
    if ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
      set( LAPACK_WRAPPER_PLATFORM_SUFFIX mingw64 )
    else()
      set( LAPACK_WRAPPER_PLATFORM_SUFFIX mingw32 )
    endif()
  elseif ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set( LAPACK_WRAPPER_PLATFORM_SUFFIX win64 )
  else()
    set( LAPACK_WRAPPER_PLATFORM_SUFFIX win32 )
  endif()
elseif ( UNIX )
  set( LAPACK_WRAPPER_PLATFORM_SUFFIX linux )
else()
  string( TOLOWER "${CMAKE_SYSTEM_NAME}" LAPACK_WRAPPER_PLATFORM_SUFFIX )
endif()

# ============================================================================
# User-facing build options and cache variables
# ============================================================================
# This section defines all options that are expected to be controlled from the
# CMake command line, presets, CI scripts, or parent projects using add_subdirectory.
#
# Single-configuration generators, such as Unix Makefiles and Ninja, need an
# explicit CMAKE_BUILD_TYPE.  Multi-configuration generators, such as Visual
# Studio and Ninja Multi-Config, ignore CMAKE_BUILD_TYPE and select the build
# configuration at build time instead.
#
# BUILD_SHARED_LIBS follows the CMake convention: add_library() calls without an
# explicit STATIC/SHARED keyword are shared when this option is ON and static
# when it is OFF.  Both LapackWrapper and HSL use this convention below.
#
# LAPACK_WRAPPER_BACKEND chooses the BLAS/LAPACK provider.  It intentionally is a
# STRING cache variable, not several booleans, so exactly one backend is active.
# ============================================================================
if ( NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES )
  set( CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE )
  set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel" )
endif()

option( BUILD_SHARED_LIBS "Build LapackWrapper as a shared library" OFF )
option( BUILD_TESTING "Build LapackWrapper test executables" ${PROJECT_IS_TOP_LEVEL} )
option( LAPACK_WRAPPER_INSTALL "Generate install target" ${PROJECT_IS_TOP_LEVEL} )
option( LAPACK_WRAPPER_STRICT_WARNINGS "Enable strict warnings" OFF )
option( LAPACK_WRAPPER_FMT_HEADER_ONLY "Use fmt in header-only mode" ON )

set(
  LAPACK_WRAPPER_BACKEND "OpenBLAS"
  CACHE STRING "BLAS/LAPACK backend: Accelerate, OpenBLAS, MKL, Netlib"
)
set_property(
  CACHE LAPACK_WRAPPER_BACKEND
  PROPERTY STRINGS "Accelerate" "OpenBLAS" "MKL" "Netlib"
)

set(
  LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite"
  CACHE PATH "Local UtilsLite checkout"
)

set(
  LAPACK_WRAPPER_EIGEN_LOCAL_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/../Eigen"
  CACHE PATH "Local Eigen checkout"
)

set(
  LAPACK_WRAPPER_EIGEN_VERSION
  "5.0.1"
  CACHE STRING "Eigen version to fetch if no local checkout is found"
)

# Default local install layout:
#   <source>/lib/include
#   <source>/lib/lib
#
# The project keeps the historical local layout used by the surrounding build
# scripts: headers and libraries are installed inside the source tree under
# lib/.  A user can still override CMAKE_INSTALL_PREFIX from the command line.
# GNUInstallDirs is included later so the standard variables remain available.
if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE )
endif()
set( CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE )
set( CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE )

# Language mode.
# The wrapper requires C++20.  The standard is stored in the cache with FORCE so
# stale cache entries inherited from older cmake_utils-based builds cannot keep
# the project on C++17.  CXX_EXTENSIONS is disabled to request a strict standard
# mode such as -std=c++20 rather than compiler extensions such as -std=gnu++20.
set( CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard" FORCE )
set( CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require selected C++ standard" FORCE )
set( CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Disable compiler-specific C++ extensions" FORCE )

# GNUInstallDirs provides standard install-directory variables such as
# CMAKE_INSTALL_INCLUDEDIR and CMAKE_INSTALL_LIBDIR.  FetchContent is used to
# bring source dependencies into the build when a local checkout is not found.
include( GNUInstallDirs )
include( FetchContent )

# Do not set CMAKE_OSX_ARCHITECTURES here.  Apple universal/fat builds are
# configured before project() by the UTILS_BUILD_FAT block above.  Reassigning
# CMAKE_OSX_ARCHITECTURES after include()/project() would either be too late for
# already-created compiler state or would narrow a requested universal build
# back to the host architecture only.

# ============================================================================
# Dependency helper: prefer a local checkout, otherwise use FetchContent
# ============================================================================
# lapack_wrapper_use_local_or_fetch(NAME LOCAL_DIR) centralizes the dependency
# policy used for source dependencies that provide their own CMake project.
#
# Resolution order:
#   1. If LOCAL_DIR contains a CMakeLists.txt, expose it through the standard
#      FETCHCONTENT_SOURCE_DIR_<NAME> cache variable and call
#      FetchContent_MakeAvailable(NAME).  This makes local sibling checkouts
#      reproducible and avoids unnecessary downloads.
#   2. Otherwise call FetchContent_MakeAvailable(NAME), using the declaration
#      created by the caller.
#
# The helper temporarily sets CMAKE_SKIP_INSTALL_RULES while the dependency is
# added.  This avoids exporting or installing dependency rules as part of the
# LapackWrapper install target.  A small cmake_install.cmake stub is also written
# for the dependency build directory because some workflows expect it to exist.
#
# The function reports two variables to the parent scope:
#   <NAME>_SOURCE     - human-readable source summary for the final report.
#   <NAME>_SOURCE_DIR - actual source directory, when known.
# ============================================================================
function( lapack_wrapper_use_local_or_fetch NAME LOCAL_DIR )
  # Convert the user-provided path to an absolute path relative to this project.
  # This makes the behavior independent of the caller's current working
  # directory and keeps cache entries stable across configure invocations.
  get_filename_component(
    _local_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  string( TOUPPER "${NAME}" _name_upper )
  string( TOLOWER "${NAME}" _name_lower )

  if ( FETCHCONTENT_BASE_DIR )
    set( _fetchcontent_base_dir "${FETCHCONTENT_BASE_DIR}" )
  else()
    set( _fetchcontent_base_dir "${CMAKE_BINARY_DIR}/_deps" )
  endif()

  set( _dependency_install_stub "${_fetchcontent_base_dir}/${_name_lower}-build/cmake_install.cmake" )
  set( _saved_skip_install_rules "${CMAKE_SKIP_INSTALL_RULES}" )
  set( CMAKE_SKIP_INSTALL_RULES ON )

  if ( EXISTS "${_local_abs}/CMakeLists.txt" )
    message( STATUS "Using local ${NAME}: ${_local_abs}" )
    set( "FETCHCONTENT_SOURCE_DIR_${_name_upper}" "${_local_abs}" CACHE PATH "Local ${NAME} checkout" FORCE )
    FetchContent_MakeAvailable( ${NAME} )
    set( "${NAME}_SOURCE" "local: ${_local_abs}" PARENT_SCOPE )
    set( "${NAME}_SOURCE_DIR" "${_local_abs}" PARENT_SCOPE )
  else()
    message( STATUS "Fetching ${NAME}" )
    FetchContent_MakeAvailable( ${NAME} )
    set( "${NAME}_SOURCE" "FetchContent" PARENT_SCOPE )
    if ( DEFINED "${_name_lower}_SOURCE_DIR" )
      set( "${NAME}_SOURCE_DIR" "${${_name_lower}_SOURCE_DIR}" PARENT_SCOPE )
    endif()
  endif()

  file( MAKE_DIRECTORY "${_fetchcontent_base_dir}/${_name_lower}-build" )
  if ( NOT EXISTS "${_dependency_install_stub}" )
    file(
      WRITE "${_dependency_install_stub}"
      "# Install disabled for dependency ${NAME} inside LapackWrapper.\n"
    )
  endif()

  set( CMAKE_SKIP_INSTALL_RULES "${_saved_skip_install_rules}" )
endfunction()

# ============================================================================
# BLAS/LAPACK backend selection
# ============================================================================
# The wrapper can be built against several numerical backends.  This section
# converts LAPACK_WRAPPER_BACKEND into three pieces of target data used later:
#
#   LAPACK_WRAPPER_LAPACK_LIBRARIES      libraries or imported targets to link.
#   LAPACK_WRAPPER_COMPILE_DEFINITIONS   preprocessor flags for C++ sources.
#   LAPACK_WRAPPER_CONFIG_DEFINE         text emitted into lapack_wrapper_config.hh.
#
# Accelerate is restricted to Apple platforms.  OpenBLAS is located with
# find_library() to support system package managers and custom library paths.
# MKL uses find_package(MKL) and expects an imported MKL::MKL target.  Netlib
# uses CMake's standard BLAS and LAPACK packages.
# ============================================================================
set( LAPACK_WRAPPER_LAPACK_LIBRARIES "" )
set( LAPACK_WRAPPER_COMPILE_DEFINITIONS "" )
set( LAPACK_WRAPPER_CONFIG_DEFINE "" )

# Select exactly one backend.  Each branch populates the common variables that
# are later attached to HSL and LapackWrapper as usage requirements.
if ( LAPACK_WRAPPER_BACKEND STREQUAL "Accelerate" )

  if ( NOT APPLE )
    message( FATAL_ERROR "Accelerate framework is only available on macOS" )
  endif()

  set( LAPACK_WRAPPER_LAPACK_LIBRARIES "-framework Accelerate" )
  list( APPEND LAPACK_WRAPPER_COMPILE_DEFINITIONS LAPACK_WRAPPER_USE_ACCELERATE )
  set( LAPACK_WRAPPER_CONFIG_DEFINE "#define LAPACK_WRAPPER_USE_ACCELERATE 1" )

elseif ( LAPACK_WRAPPER_BACKEND STREQUAL "OpenBLAS" )

  find_library( OPENBLAS_LIB NAMES openblas libopenblas )
  if ( NOT OPENBLAS_LIB )
    message( FATAL_ERROR "OpenBLAS library not found. Please install OpenBLAS or use -DLAPACK_WRAPPER_BACKEND=Accelerate on macOS." )
  endif()

  find_library( LAPACKE_LIB NAMES lapacke liblapacke )
  if ( LAPACKE_LIB )
    set( LAPACK_WRAPPER_LAPACK_LIBRARIES ${OPENBLAS_LIB} ${LAPACKE_LIB} )
  else()
    set( LAPACK_WRAPPER_LAPACK_LIBRARIES ${OPENBLAS_LIB} )
  endif()

  list( APPEND LAPACK_WRAPPER_COMPILE_DEFINITIONS LAPACK_WRAPPER_USE_OPENBLAS )
  set( LAPACK_WRAPPER_CONFIG_DEFINE "#define LAPACK_WRAPPER_USE_OPENBLAS 1" )

elseif ( LAPACK_WRAPPER_BACKEND STREQUAL "MKL" )

  find_package( MKL REQUIRED )
  if ( NOT TARGET MKL::MKL )
    message( FATAL_ERROR "MKL not found. Please set MKL_ROOT or install MKL." )
  endif()

  set( LAPACK_WRAPPER_LAPACK_LIBRARIES MKL::MKL )
  list( APPEND LAPACK_WRAPPER_COMPILE_DEFINITIONS LAPACK_WRAPPER_USE_MKL )
  set( LAPACK_WRAPPER_CONFIG_DEFINE "#define LAPACK_WRAPPER_USE_MKL 1" )

elseif ( LAPACK_WRAPPER_BACKEND STREQUAL "Netlib" )

  find_package( BLAS REQUIRED )
  find_package( LAPACK REQUIRED )
  set( LAPACK_WRAPPER_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  list( APPEND LAPACK_WRAPPER_COMPILE_DEFINITIONS LAPACK_WRAPPER_USE_NETLIB )
  set( LAPACK_WRAPPER_CONFIG_DEFINE "#define LAPACK_WRAPPER_USE_NETLIB 1" )

else()

  message( FATAL_ERROR "Unknown LAPACK_WRAPPER_BACKEND: ${LAPACK_WRAPPER_BACKEND}" )

endif()

# ============================================================================
# UtilsLite dependency resolution
# ============================================================================
# LapackWrapper depends on UtilsLite headers and, when available, on a UtilsLite
# library target or a prebuilt UtilsLite archive.  This section supports three
# workflows without requiring the old cmake_utils infrastructure:
#
#   1. Use a prebuilt sibling UtilsLite library from ../UtilsLite/lib/lib.
#   2. Add a local sibling UtilsLite checkout via FetchContent source override.
#   3. Fetch UtilsLite from GitHub when no local copy is available.
#
# UtilsLite install rules and tests are disabled before the dependency is added,
# so a parent LapackWrapper install does not accidentally install UtilsLite as a
# separate package and a dependency build does not spend time compiling tests.
#
# If no concrete UtilsLite target can be found after resolution, the fallback is
# an INTERFACE target that contributes only include directories.  This keeps the
# build usable for header-only portions of UtilsLite or partial local setups.
# ============================================================================
set( UTILS_INSTALL_LOCAL OFF CACHE BOOL "Disable UtilsLite install rules inside LapackWrapper" FORCE )
set( UTILS_ENABLE_TESTS OFF CACHE BOOL "Disable UtilsLite tests inside LapackWrapper" FORCE )

set( LAPACK_WRAPPER_UTILS_TARGET "" )
set( LAPACK_WRAPPER_UTILSLITE_SOURCE "not found" )

get_filename_component(
  LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS
  "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR}"
  ABSOLUTE
  BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Prefer an already built sibling UtilsLite library when present.  This mirrors
# the old large-project workflow without requiring cmake_utils functions.  The
# lookup uses the same platform suffix convention defined above, so static and
# shared UtilsLite artifacts can be discovered under the local lib directories.
if ( EXISTS "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}" )
  find_library(
    LAPACK_WRAPPER_UTILSLITE_IMPORTED_LIB
    NAMES
      UtilsLite_${LAPACK_WRAPPER_PLATFORM_SUFFIX}_static
      UtilsLite_${LAPACK_WRAPPER_PLATFORM_SUFFIX}
      UtilsLite
    PATHS
      "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/lib/lib"
      "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/lib"
    NO_DEFAULT_PATH
  )
endif()

if ( LAPACK_WRAPPER_UTILSLITE_IMPORTED_LIB )

  message( STATUS "Using prebuilt UtilsLite: ${LAPACK_WRAPPER_UTILSLITE_IMPORTED_LIB}" )
  set( LAPACK_WRAPPER_UTILSLITE_SOURCE "prebuilt: ${LAPACK_WRAPPER_UTILSLITE_IMPORTED_LIB}" )

  add_library( lapack_wrapper_UtilsLite_imported UNKNOWN IMPORTED GLOBAL )
  set_target_properties(
    lapack_wrapper_UtilsLite_imported PROPERTIES
    IMPORTED_LOCATION "${LAPACK_WRAPPER_UTILSLITE_IMPORTED_LIB}"
  )

  target_include_directories(
    lapack_wrapper_UtilsLite_imported SYSTEM INTERFACE
    "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/src"
    "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/src/Utils/3rd"
    "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/lib/include"
  )

  set( LAPACK_WRAPPER_UTILS_TARGET lapack_wrapper_UtilsLite_imported )
  set( UtilsLite_SOURCE_DIR "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}" )

else()

  FetchContent_Declare(
    UtilsLite
    GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
    GIT_TAG        Utils2.0
    GIT_SHALLOW    TRUE
  )
  lapack_wrapper_use_local_or_fetch( UtilsLite "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR}" )
  set( LAPACK_WRAPPER_UTILSLITE_SOURCE "${UtilsLite_SOURCE}" )

  foreach( _utils_candidate
    utils::UtilsLite_Static
    utils_UtilsLite_Static
    UtilsLite::UtilsLite
    UtilsLite
    utils_UtilsLite
  )
    if ( TARGET ${_utils_candidate} AND NOT LAPACK_WRAPPER_UTILS_TARGET )
      set( LAPACK_WRAPPER_UTILS_TARGET ${_utils_candidate} )
    endif()
  endforeach()

endif()

if ( NOT LAPACK_WRAPPER_UTILS_TARGET )
  message( WARNING "No UtilsLite library target found. Building with include directories only." )
  add_library( lapack_wrapper_UtilsLite_headers INTERFACE )
  if ( EXISTS "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/src" )
    target_include_directories(
      lapack_wrapper_UtilsLite_headers SYSTEM INTERFACE
      "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/src"
      "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/src/Utils/3rd"
      "${LAPACK_WRAPPER_UTILSLITE_LOCAL_DIR_ABS}/lib/include"
    )
  endif()
  set( LAPACK_WRAPPER_UTILS_TARGET lapack_wrapper_UtilsLite_headers )
endif()

foreach( _utils_tgt
  ${LAPACK_WRAPPER_UTILS_TARGET}
  utils_UtilsLite_Static
  utils_UtilsLite
)
  if ( TARGET ${_utils_tgt} )
    get_target_property( _utils_inc ${_utils_tgt} INTERFACE_INCLUDE_DIRECTORIES )
    if ( _utils_inc )
      set_target_properties(
        ${_utils_tgt} PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_utils_inc}"
      )
    endif()
  endif()
endforeach()

# ============================================================================
# Eigen dependency, optional and header-only
# ============================================================================
# Eigen is consumed as the imported INTERFACE target Eigen3::Eigen.  Because
# Eigen is header-only for this use case, creating an imported interface target
# is sufficient when a local or bundled include tree is found.
#
# Resolution order:
#   1. Reuse an existing Eigen3::Eigen target provided by a parent project.
#   2. Use LAPACK_WRAPPER_EIGEN_LOCAL_DIR when it contains Eigen/Core.
#   3. Use the bundled lib3rd/include copy when present.
#   4. Fetch the requested Eigen archive with FetchContent.
#
# The final summary records which path was chosen so CI logs and user reports
# make the dependency provenance clear.
# ============================================================================
set( LAPACK_WRAPPER_EIGEN_SOURCE "not used" )

if ( NOT TARGET Eigen3::Eigen )
  get_filename_component(
    LAPACK_WRAPPER_EIGEN_LOCAL_DIR_ABS
    "${LAPACK_WRAPPER_EIGEN_LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  if ( EXISTS "${LAPACK_WRAPPER_EIGEN_LOCAL_DIR_ABS}/Eigen/Core" )
    message( STATUS "Using local Eigen checkout: ${LAPACK_WRAPPER_EIGEN_LOCAL_DIR_ABS}" )
    set( LAPACK_WRAPPER_EIGEN_SOURCE "local: ${LAPACK_WRAPPER_EIGEN_LOCAL_DIR_ABS}" )
    add_library( Eigen3::Eigen INTERFACE IMPORTED GLOBAL )
    target_include_directories( Eigen3::Eigen SYSTEM INTERFACE "${LAPACK_WRAPPER_EIGEN_LOCAL_DIR_ABS}" )
  elseif ( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include/Eigen/Core" )
    message( STATUS "Using bundled Eigen: ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include" )
    set( LAPACK_WRAPPER_EIGEN_SOURCE "bundled: lib3rd/include" )
    add_library( Eigen3::Eigen INTERFACE IMPORTED GLOBAL )
    target_include_directories( Eigen3::Eigen SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include" )
  else()
    message( STATUS "Fetching Eigen ${LAPACK_WRAPPER_EIGEN_VERSION}" )
    set( LAPACK_WRAPPER_EIGEN_SOURCE "FetchContent: ${LAPACK_WRAPPER_EIGEN_VERSION}" )
    FetchContent_Declare(
      eigen3
      URL "https://gitlab.com/libeigen/eigen/-/archive/${LAPACK_WRAPPER_EIGEN_VERSION}/eigen-${LAPACK_WRAPPER_EIGEN_VERSION}.tar.gz"
      SOURCE_SUBDIR cmake
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable( eigen3 )
    if ( NOT TARGET Eigen3::Eigen )
      add_library( Eigen3::Eigen INTERFACE IMPORTED GLOBAL )
      target_include_directories( Eigen3::Eigen SYSTEM INTERFACE "${eigen3_SOURCE_DIR}" )
    endif()
  endif()
else()
  message( STATUS "Using existing Eigen3::Eigen target" )
  set( LAPACK_WRAPPER_EIGEN_SOURCE "existing target: Eigen3::Eigen" )
endif()

# ============================================================================
# Source and public-header collection
# ============================================================================
# The project source tree is intentionally discovered with GLOB_RECURSE and
# CONFIGURE_DEPENDS.  This keeps the CMake file short while still asking CMake
# to reconfigure when matching source/header files are added or removed.
#
# hsl_fake.cc is excluded from the main LapackWrapper source list because it is
# built as the separate HSL target below.  Keeping HSL as its own target makes
# the link structure explicit and preserves the public LapackWrapper::HSL alias.
#
# Public headers include .cxx files because some template or inline
# implementation files are intended to be installed with the headers.
# ============================================================================
set( SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )

file(
  GLOB_RECURSE LAPACK_WRAPPER_SOURCES
  CONFIGURE_DEPENDS
  "${SOURCE_DIR}/*.cc"
)
list( FILTER LAPACK_WRAPPER_SOURCES EXCLUDE REGEX ".*hsl_fake\\.cc$" )

set( LAPACK_WRAPPER_HSL_SOURCE "${SOURCE_DIR}/HSL/hsl_fake.cc" )
if ( NOT EXISTS "${LAPACK_WRAPPER_HSL_SOURCE}" )
  message( FATAL_ERROR "Missing HSL fake source: ${LAPACK_WRAPPER_HSL_SOURCE}" )
endif()

file(
  GLOB_RECURSE LAPACK_WRAPPER_PUBLIC_HEADERS
  CONFIGURE_DEPENDS
  "${SOURCE_DIR}/*.h"
  "${SOURCE_DIR}/*.hh"
  "${SOURCE_DIR}/*.hpp"
  "${SOURCE_DIR}/*.hxx"
  "${SOURCE_DIR}/*.cxx"
)
list( FILTER LAPACK_WRAPPER_PUBLIC_HEADERS EXCLUDE REGEX ".*\\.tmpl$" )

# ============================================================================
# Generated backend configuration header
# ============================================================================
# lapack_wrapper_config.hh exposes the selected BLAS/LAPACK backend to the C++
# sources.  When a template exists, placeholder tokens are replaced with the
# selected backend define.  If the template is absent but no generated header is
# available, a minimal header containing only the active backend define is
# written so the build can still proceed from lean source snapshots.
#
# The OpenBLAS NOSYSTEM placeholder is currently replaced with an empty string;
# the hook remains in place for compatibility with older templates.
# ============================================================================
set( LAPACK_WRAPPER_CONFIG_TEMPLATE "${SOURCE_DIR}/lapack_wrapper_config.hh.tmpl" )
set( LAPACK_WRAPPER_CONFIG_OUTPUT   "${SOURCE_DIR}/lapack_wrapper_config.hh" )

if ( EXISTS "${LAPACK_WRAPPER_CONFIG_TEMPLATE}" )
  file( READ "${LAPACK_WRAPPER_CONFIG_TEMPLATE}" LAPACK_WRAPPER_CONFIG_CONTENT )
  string(
    REPLACE
    "@@LAPACK_WRAPPER_USE@@"
    "${LAPACK_WRAPPER_CONFIG_DEFINE}"
    LAPACK_WRAPPER_CONFIG_CONTENT
    "${LAPACK_WRAPPER_CONFIG_CONTENT}"
  )
  string(
    REPLACE
    "@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@"
    ""
    LAPACK_WRAPPER_CONFIG_CONTENT
    "${LAPACK_WRAPPER_CONFIG_CONTENT}"
  )
  file( WRITE "${LAPACK_WRAPPER_CONFIG_OUTPUT}" "${LAPACK_WRAPPER_CONFIG_CONTENT}" )
else()
  message( STATUS "Config template not found; keeping existing lapack_wrapper_config.hh if present." )
  if ( NOT EXISTS "${LAPACK_WRAPPER_CONFIG_OUTPUT}" )
    file( WRITE "${LAPACK_WRAPPER_CONFIG_OUTPUT}" "${LAPACK_WRAPPER_CONFIG_DEFINE}\n" )
  endif()
endif()

# ============================================================================
# Output-directory helpers
# ============================================================================
# The main project writes its artifacts into the historical local layout:
#   libraries -> <source>/lib/lib
#   executables -> <source>/bin
#
# Dependencies that are built as part of this project are redirected to
# lib3rd/lib so they do not mix with the primary LapackWrapper and HSL artifacts.
# Alias targets are resolved to their real targets before properties are set.
# INTERFACE libraries are skipped because they do not produce binary artifacts.
# ============================================================================
function( lapack_wrapper_set_target_output_dirs TARGET_NAME LIB_DIR )
  # TARGET_NAME is a real binary target owned by this project.  LIB_DIR is a
  # source-relative directory such as lib/lib.  Runtime output is kept in bin so
  # local test and example commands can use a predictable executable path.
  set_target_properties(
    ${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

function( lapack_wrapper_set_dependency_output_dirs LIB_DIR )
  # ARGN contains candidate dependency targets.  Some entries may be aliases or
  # may not exist in the current configuration, so each entry is checked before
  # output-directory properties are applied.
  foreach( tgt IN LISTS ARGN )
    if ( TARGET ${tgt} )
      get_target_property( real_tgt ${tgt} ALIASED_TARGET )
      if ( NOT real_tgt )
        set( real_tgt ${tgt} )
      endif()
      if ( TARGET ${real_tgt} )
        get_target_property( tgt_type ${real_tgt} TYPE )
        if ( NOT tgt_type STREQUAL "INTERFACE_LIBRARY" )
          set_target_properties(
            ${real_tgt} PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
          )
        endif()
      endif()
    endif()
  endforeach()
endfunction()

# ============================================================================
# HSL compatibility target
# ============================================================================
# HSL is built from the local hsl_fake.cc source and exported as
# LapackWrapper::HSL.  This provides a stable link target for code that expects
# the HSL component while still allowing this repository to build without the
# proprietary or external HSL implementation.
#
# The target receives the same backend compile definitions as LapackWrapper, and
# it links UtilsLite because the HSL shim can depend on shared utility headers or
# implementation details.  POSITION_INDEPENDENT_CODE is enabled so the target is
# usable from both static and shared builds.
# ============================================================================
# HSL is intentionally separate from LapackWrapper so consumers can link or
# inspect the compatibility component independently.
add_library( HSL ${LAPACK_WRAPPER_HSL_SOURCE} )
add_library( LapackWrapper::HSL ALIAS HSL )

lapack_wrapper_set_target_output_dirs( HSL lib/lib )

target_compile_features( HSL PUBLIC cxx_std_20 )
target_compile_definitions( HSL PUBLIC _USE_MATH_DEFINES ${LAPACK_WRAPPER_COMPILE_DEFINITIONS} )
if ( LAPACK_WRAPPER_FMT_HEADER_ONLY )
  target_compile_definitions( HSL PUBLIC FMT_HEADER_ONLY )
endif()

target_include_directories(
  HSL PUBLIC
  "$<BUILD_INTERFACE:${SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(
  HSL SYSTEM PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include>"
)

target_link_libraries( HSL PUBLIC ${LAPACK_WRAPPER_UTILS_TARGET} )

set_target_properties(
  HSL PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME HSL
  OUTPUT_NAME HSL
  POSITION_INDEPENDENT_CODE ON
)

# ============================================================================
# Main LapackWrapper library target
# ============================================================================
# This is the primary public target and is also exposed through the namespace
# alias LapackWrapper::LapackWrapper.  Because add_library() does not specify
# STATIC or SHARED, the target type follows BUILD_SHARED_LIBS.
#
# Public usage requirements include the source include directory during the
# build and the install include directory for installed packages.  Third-party
# bundled headers are marked SYSTEM to suppress warnings from external code.
#
# The link interface is PUBLIC because consumers of LapackWrapper may need the
# same BLAS/LAPACK, Eigen, HSL, and UtilsLite usage requirements when compiling
# or linking their own targets.
# ============================================================================
# The main library type follows BUILD_SHARED_LIBS because no STATIC/SHARED
# keyword is supplied here.
add_library( LapackWrapper ${LAPACK_WRAPPER_SOURCES} )
add_library( LapackWrapper::LapackWrapper ALIAS LapackWrapper )

lapack_wrapper_set_target_output_dirs( LapackWrapper lib/lib )

target_compile_features( LapackWrapper PUBLIC cxx_std_20 )
target_compile_definitions(
  LapackWrapper PUBLIC
  _USE_MATH_DEFINES
  ${LAPACK_WRAPPER_COMPILE_DEFINITIONS}
)
if ( LAPACK_WRAPPER_FMT_HEADER_ONLY )
  target_compile_definitions( LapackWrapper PUBLIC FMT_HEADER_ONLY )
endif()

target_include_directories(
  LapackWrapper PUBLIC
  "$<BUILD_INTERFACE:${SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(
  LapackWrapper SYSTEM PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include>"
)

target_link_libraries(
  LapackWrapper PUBLIC
  LapackWrapper::HSL
  ${LAPACK_WRAPPER_UTILS_TARGET}
  Eigen3::Eigen
  ${LAPACK_WRAPPER_LAPACK_LIBRARIES}
)

set_target_properties(
  LapackWrapper PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME LapackWrapper
  OUTPUT_NAME lapack_wrapper
  POSITION_INDEPENDENT_CODE ON
)

# ============================================================================
# Compiler warnings and platform-specific compile/link flags
# ============================================================================
# Warning policy is target-scoped and applied only to the project-owned binary
# targets.  Normal builds use -Wall -Wextra on GCC/Clang-family compilers.
# LAPACK_WRAPPER_STRICT_WARNINGS adds conversion and shadow diagnostics for
# development builds where stricter code hygiene is desired.
#
# Clang builds explicitly use libc++ for both compilation and linking.  MSVC
# receives UTF-8 and big-object flags to support large generated or template-heavy
# translation units.
# ============================================================================
foreach( _target LapackWrapper HSL )
  if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU" )
    target_compile_options( ${_target} PRIVATE -Wall -Wextra )
    if ( LAPACK_WRAPPER_STRICT_WARNINGS )
      target_compile_options(
        ${_target} PRIVATE
        -Wconversion
        -Wsign-conversion
        -Wshadow
        -Wdouble-promotion
      )
    endif()
  endif()

  if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    target_compile_options( ${_target} PUBLIC -stdlib=libc++ )
    target_link_options( ${_target} PUBLIC -stdlib=libc++ )
  endif()

  if ( MSVC )
    target_compile_options( ${_target} PUBLIC /utf-8 /bigobj )
  endif()
endforeach()

# Put dependency artifacts, if any, under lib3rd/lib.  This avoids mixing
# third-party outputs with the public LapackWrapper/HSL artifacts under lib/lib.
lapack_wrapper_set_dependency_output_dirs(
  lib3rd/lib
  ${LAPACK_WRAPPER_UTILS_TARGET}
)

# ============================================================================
# Platform-qualified library names and compatibility aliases
# ============================================================================
# The real binary artifact receives a platform suffix, for example:
#   liblapack_wrapper_osx_static.a
#   libHSL_linux_static.a
#
# For backward compatibility, an unsuffixed file name is also created:
#   liblapack_wrapper.a -> liblapack_wrapper_osx_static.a
#
# On POSIX platforms the compatibility name is a symlink.  On Windows it is a
# copied file because symlink creation is less portable and may require elevated
# privileges.  The same alias operation is repeated at install time so installed
# trees expose both the explicit and legacy file names.
# ============================================================================
function( lapack_wrapper_add_platform_library_alias TARGET_NAME BASE_NAME )
  # TARGET_NAME is the CMake target to rename.  BASE_NAME is the unsuffixed stem
  # used for the compatibility alias.  The helper is applied to both the main
  # LapackWrapper library and the HSL compatibility library.
  if ( NOT TARGET ${TARGET_NAME} )
    return()
  endif()

  # The real library keeps the platform/static suffix, e.g.
  #   liblapack_wrapper_osx_static.a
  # and the unsuffixed name is kept only as a compatibility alias, e.g.
  #   liblapack_wrapper.a -> liblapack_wrapper_osx_static.a
  if ( BUILD_SHARED_LIBS )
    set( OUTPUT_NAME "${BASE_NAME}_${LAPACK_WRAPPER_PLATFORM_SUFFIX}" )
    set( ALIAS_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
    if ( WIN32 )
      set( REAL_FILENAME "lib${OUTPUT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
      set( OUTPUT_NAME "lib${OUTPUT_NAME}" )
    else()
      set( REAL_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${OUTPUT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
    endif()
  else()
    set( OUTPUT_NAME "${BASE_NAME}_${LAPACK_WRAPPER_PLATFORM_SUFFIX}_static" )
    set( ALIAS_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    if ( WIN32 )
      set( REAL_FILENAME "lib${OUTPUT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
      set( OUTPUT_NAME "lib${OUTPUT_NAME}" )
    else()
      set( REAL_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${OUTPUT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    endif()
  endif()

  set_target_properties(
    ${TARGET_NAME} PROPERTIES
    OUTPUT_NAME         "${OUTPUT_NAME}"
    ARCHIVE_OUTPUT_NAME "${OUTPUT_NAME}"
    LIBRARY_OUTPUT_NAME "${OUTPUT_NAME}"
    RUNTIME_OUTPUT_NAME "${OUTPUT_NAME}"
  )

  if ( WIN32 )
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "$<TARGET_FILE:${TARGET_NAME}>"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
      VERBATIM
    )

    if ( LAPACK_WRAPPER_INSTALL )
      install(
        CODE
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")\n
        execute_process(\n
          COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different\n
                  \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${REAL_FILENAME}\"\n
                  \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${ALIAS_FILENAME}\"\n
        )"
      )
    endif()
  else()
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "$<TARGET_FILE_NAME:${TARGET_NAME}>"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
      VERBATIM
    )

    if ( LAPACK_WRAPPER_INSTALL )
      install(
        CODE
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")\n
        execute_process(\n
          COMMAND \"${CMAKE_COMMAND}\" -E rm -f\n
                  \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${ALIAS_FILENAME}\"\n
        )\n
        execute_process(\n
          COMMAND \"${CMAKE_COMMAND}\" -E create_symlink\n
                  \"${REAL_FILENAME}\"\n
                  \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${ALIAS_FILENAME}\"\n
        )"
      )
    endif()
  endif()
endfunction()


# Apply the final artifact naming policy after the targets have been fully
# declared.  This ensures both the real platform-qualified name and the legacy
# alias are created for every produced library.
lapack_wrapper_add_platform_library_alias( LapackWrapper lapack_wrapper )
lapack_wrapper_add_platform_library_alias( HSL HSL )

# ============================================================================
# Test executables
# ============================================================================
# Tests are built only when BUILD_TESTING is enabled and this project is the top
# level project.  This prevents a parent project that includes LapackWrapper via
# FetchContent/add_subdirectory from unexpectedly compiling the entire local test
# suite.
#
# Test sources are discovered from src_tests/test*.cc, linked against the public
# LapackWrapper target, emitted into the source-level bin directory, and
# registered with CTest.  A convenience aggregate target named
# <project>_all_tests is created for local workflows.
# ============================================================================
include( CTest )

if ( BUILD_TESTING AND PROJECT_IS_TOP_LEVEL )
  file(
    GLOB LAPACK_WRAPPER_TEST_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/test*.cc"
  )
  list( SORT LAPACK_WRAPPER_TEST_SOURCES )

  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )

  foreach( TEST_SOURCE ${LAPACK_WRAPPER_TEST_SOURCES} )
    get_filename_component( TEST_NAME "${TEST_SOURCE}" NAME_WE )

    add_executable( ${TEST_NAME} "${TEST_SOURCE}" )
    target_compile_features( ${TEST_NAME} PRIVATE cxx_std_20 )

    if ( LAPACK_WRAPPER_FMT_HEADER_ONLY )
      target_compile_definitions( ${TEST_NAME} PRIVATE FMT_HEADER_ONLY )
    endif()

    target_link_libraries(
      ${TEST_NAME} PRIVATE
      LapackWrapper::LapackWrapper
      ${CMAKE_DL_LIBS}
    )

    set_target_properties(
      ${TEST_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
      SUFFIX ".exe"
    )

    add_test(
      NAME "${TEST_NAME}"
      COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${TEST_NAME}.exe"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    add_dependencies( "${PROJECT_NAME}_all_tests" ${TEST_NAME} )
  endforeach()
endif()

# ============================================================================
# Install rules and package archive setup
# ============================================================================
# Installation is optional and controlled by LAPACK_WRAPPER_INSTALL.  The install
# target exports the two project libraries and installs all public headers while
# excluding template files that are only needed to generate configure-time output.
#
# CPack is configured inside this block so packaging is available only when
# install rules are enabled.  Both binary and source archives are generated in
# TGZ and ZIP formats, while build directories and .git metadata are excluded
# from source packages.
# ============================================================================
if ( LAPACK_WRAPPER_INSTALL )
  install(
    TARGETS LapackWrapper HSL
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(
    DIRECTORY "${SOURCE_DIR}/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hh"
      PATTERN "*.hpp"
      PATTERN "*.hxx"
      PATTERN "*.cxx"
      PATTERN "*.tmpl" EXCLUDE
  )

  set( CPACK_PACKAGE_NAME "${PROJECT_NAME}" )
  set( CPACK_PACKAGE_VERSION "${PROJECT_VERSION}" )
  set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}" )
  set( CPACK_GENERATOR "TGZ;ZIP" )
  set( CPACK_SOURCE_GENERATOR "TGZ;ZIP" )
  set( CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/" )
  include( CPack )
endif()

# ============================================================================
# Configuration summary
# ============================================================================
# Print the decisions made during configuration.  This is deliberately concise
# but includes the backend, resolved libraries, dependency sources, target names,
# build type toggles, and language mode.  The summary is useful in CI logs and
# when diagnosing user reports because it records the effective dependency and
# option state after all cache variables and auto-detection have been applied.
# ============================================================================
message( STATUS "" )
message( STATUS "LapackWrapper configuration summary:" )
message( STATUS "  Version           : ${PROJECT_VERSION}" )
message( STATUS "  Backend           : ${LAPACK_WRAPPER_BACKEND}" )
message( STATUS "  Libraries         : ${LAPACK_WRAPPER_LAPACK_LIBRARIES}" )
message( STATUS "  UtilsLite         : ${LAPACK_WRAPPER_UTILSLITE_SOURCE}" )
message( STATUS "  UtilsLite target  : ${LAPACK_WRAPPER_UTILS_TARGET}" )
message( STATUS "  Eigen             : ${LAPACK_WRAPPER_EIGEN_SOURCE}" )
message( STATUS "  Shared libs       : ${BUILD_SHARED_LIBS}" )
message( STATUS "  Tests             : ${BUILD_TESTING}" )
message( STATUS "  Install           : ${LAPACK_WRAPPER_INSTALL}" )
message( STATUS "  C++ standard      : ${CMAKE_CXX_STANDARD}" )
message( STATUS "  fmt header-only   : ${LAPACK_WRAPPER_FMT_HEADER_ONLY}" )
message( STATUS "" )
