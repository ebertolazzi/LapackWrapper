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

# ============================================================================
# Version from Git tag, if available
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

project(
  LapackWrapper
  VERSION ${GIT_VERSION}
  DESCRIPTION "A C++ wrapper for BLAS and LAPACK"
  HOMEPAGE_URL "https://ebertolazzi.github.io/LapackWrapper/"
  LANGUAGES CXX C
)

# ============================================================================
# Platform suffix
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
# Options
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
  "5.0.0"
  CACHE STRING "Eigen version to fetch if no local checkout is found"
)

# Default install layout:
#   <source>/lib/include
#   <source>/lib/lib
if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE )
endif()
set( CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE )
set( CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE )

# C++20.  This replaces the old cmake_utils flags that forced -std=c++17.
set( CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard" FORCE )
set( CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require selected C++ standard" FORCE )
set( CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Disable compiler-specific C++ extensions" FORCE )

include( GNUInstallDirs )
include( FetchContent )

if ( APPLE )
  execute_process(
    COMMAND uname -m
    OUTPUT_VARIABLE ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set( CMAKE_OSX_ARCHITECTURES ${ARCHITECTURE} )
endif()

# ============================================================================
# Dependency helper, modeled after Splines but local to LapackWrapper
# ============================================================================
function( lapack_wrapper_use_local_or_fetch NAME LOCAL_DIR )
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
# BLAS/LAPACK backend
# ============================================================================
set( LAPACK_WRAPPER_LAPACK_LIBRARIES "" )
set( LAPACK_WRAPPER_COMPILE_DEFINITIONS "" )
set( LAPACK_WRAPPER_CONFIG_DEFINE "" )

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
# UtilsLite
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
# the old large-project workflow without requiring cmake_utils functions.
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
# Eigen, optional.  Prefer an existing target, then a local checkout, then fetch.
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
# Sources and public headers
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
# Generate lapack_wrapper_config.hh when a template is present
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
# Output directory helper
# ============================================================================
function( lapack_wrapper_set_target_output_dirs TARGET_NAME LIB_DIR )
  set_target_properties(
    ${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

function( lapack_wrapper_set_dependency_output_dirs LIB_DIR )
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
# HSL fake target
# ============================================================================
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
# Main library
# ============================================================================
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
# Compiler warnings and platform flags
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

# Put dependency artifacts, if any, under lib3rd/lib.
lapack_wrapper_set_dependency_output_dirs(
  lib3rd/lib
  ${LAPACK_WRAPPER_UTILS_TARGET}
)

# ============================================================================
# Platform alias files, following the Splines convention
# ============================================================================
function( lapack_wrapper_add_platform_library_alias TARGET_NAME BASE_NAME )
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


lapack_wrapper_add_platform_library_alias( LapackWrapper lapack_wrapper )
lapack_wrapper_add_platform_library_alias( HSL HSL )

# ============================================================================
# Tests
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
# Install
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
# Summary
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
