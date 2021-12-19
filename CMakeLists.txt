############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@ing.unitn.it                               #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.14)

include(./cmake/CMakeLists-common.txt)

project(
  lapack_wrapper
  VERSION ${EB_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/Utils/"
)

if ( EB_INSTALL_LOCAL )
  set( CMAKE_INSTALL_PREFIX ${CMAKE_HOME_DIRECTORY} )
endif()

message("LAPACK_WRAPPER current version = v${EB_PROJECT_VERSION}")

include(./cmake/CMakeLists-cflags.txt)
include(./cmake/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
ebertolazzi_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if(BUILD_SHARED)
  ebertolazzi_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#
set( SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )

set( SOURCES )
file( GLOB S ${SOURCE_DIR}/*.cc ${SOURCE_DIR}/**/*.cc )
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})
list(FILTER SOURCES EXCLUDE REGEX ".*hsl_fake.cc$")

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB S
  ${SOURCE_DIR}/*.h*
  ${SOURCE_DIR}/*/*.h*
  ${SOURCE_DIR}/*/*/*.h*
  ${SOURCE_DIR}/*/*.cxx
  ${SOURCE_DIR}/*/*/*.cxx
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

set( PRIVATE_HEADERS )
# no private headers

list(FILTER PUBLIC_HEADERS EXCLUDE REGEX ".*.tmpl$")
list(FILTER PRIVATE_HEADERS EXCLUDE REGEX ".*.tmpl$")

# extra include directory for compilation
include_directories( ./lib3rd/include )

#       _                           _                 _
#    __| | ___ _ __   ___ _ __   __| | ___ _ __   ___(_) ___  ___
#   / _` |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \ / __| |/ _ \/ __|
#  | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (__| |  __/\__ \
#   \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___|_|\___||___/
#             |_|
#
find_library(
  SUPERLU superlu${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
)

# DEPENDENCIES (TRY TO USE ALREADY COMPILED SUBMODULE IN LARGE PROJECT)
set(UPDIR ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(DEPEND_TARGETS)
find_library(
  UTILS Utils${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${UPDIR}/Utils/lib/lib
  NO_DEFAULT_PATH
)
if( UTILS )
  message("Found Utils${ARTIFACTS_STATIC_SUFFIX}")
  include_directories(${UPDIR}/Utils/lib/include)
else()
  message("NOT Found Utils${ARTIFACTS_STATIC_SUFFIX}")
  add_subdirectory(./submodules/Utils)
  set( UTILS ebertolazziUtilsStatic )
  set( DEPEND_TARGETS ${DEPEND_TARGETS} ${UTILS} )
endif()

# definisce lapackblas_libraries
include( ./cmake/CMakeLists-blas-lapack.txt )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#
# define target
ebertolazzi_setup_target(
  ${PROJECT_NAME}
  TRUE # build static
  "${SOURCE_DIR}"
  "${SOURCES}"
  "${PUBLIC_HEADERS}"
  "${PRIVATE_HEADERS}"
)

ebertolazzi_setup_target(
  "HSL"
  TRUE # build static
  "${SOURCE_DIR}"
  "HSL/hsl_fake.cc" # must be relative path
  "${PUBLIC_HEADERS}"
  "${PRIVATE_HEADERS}"
)

if ( BUILD_SHARED )
  ebertolazzi_setup_target(
    ${PROJECT_NAME}
    FALSE # build dynamic
    "${SOURCE_DIR}"
    "${SOURCES}"
    "${PUBLIC_HEADERS}"
    "${PRIVATE_HEADERS}"
  )
  ebertolazzi_setup_target(
    "HSL"
    FALSE # build static
    "${SOURCE_DIR}"
    "HSL/hsl_fake.cc" # must be relative path
    "${PUBLIC_HEADERS}"
    "${PRIVATE_HEADERS}"
  )
endif()

#   _____       _
#  |_   _|__ __| |_ ___
#    | |/ -_|_-<  _(_-<
#    |_|\___/__/\__/__/
#
if ( BUILD_EXECUTABLE OR ENABLE_CTEST )

  enable_testing()

  set( EXELISTBIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin" )
  set( EXELISTCPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src_tests" )

  if(APPLE OR UNIX)
    find_package( ZLIB REQUIRED )
  else()
    if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
      find_library( ZLIB_LIBRARIES zlibstaticd HINTS ./ThirdParties/zlib/lib )
    else()
      find_library( ZLIB_LIBRARIES zlibstatic HINTS ./ThirdParties/zlib/lib )
    endif()
    set( ZLIB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParties/zlib/include ")
  endif()

  file( MAKE_DIRECTORY ${EXELISTBIN_DIR} )

  set(
    EXELISTCPP
    test1-small-factorization
    test2-Timing
    test3-BandedMatrix
    test4-BFGS
    test5-BLOCKTRID
    test6-EIGS
    test7-SparseTool
    test8-SparseToolComplex
    test9-SparseToolTridiagonal
    test10-SparseToolTiming
    test11-SparseToolVector
    test12-SparseToolMatrix
    test13-SparseTool1
    test14-SparseTool2
  )

  add_custom_target( all_tests ALL )

  set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  set( LIBS ${EB_NAMESPACE}_${PROJECT_NAME}_Static ${EB_NAMESPACE}_HSL_Static ${lapackblas_libraries} ${UTILS} ${SUPERLU} ${ZLIB_LIBRARIES} ${CMAKE_DL_LIBS} )
  foreach( S ${EXELISTCPP} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc )
    target_link_libraries( ${S} ${LIBS} )
    add_test(
      NAME "${S}"
      COMMAND ./bin/${S}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    add_dependencies( all_tests ${S} )
  endforeach()

endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  # EXPORT (da definire)
  LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
  PUBLIC_HEADER DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION}
)

# copy Eigen header only library
#install( DIRECTORY "${SOURCE_DIR}/Eigen" DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION} )

#if ( WIN32 )
#  set( DIR3RD ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )
#  install( DIRECTORY "${DIR3RD}/include" DESTINATION "${INSTALL_PREFIX}/" )
#  install( DIRECTORY "${DIR3RD}/dll"     DESTINATION "${INSTALL_PREFIX}/" )
#  install( DIRECTORY "${DIR3RD}/lib"     DESTINATION "${INSTALL_PREFIX}/" )
#endif()

#   ___         _
#  | _ \__ _ __| |____ _ __ _ ___
#  |  _/ _` / _| / / _` / _` / -_)
#  |_| \__,_\__|_\_\__,_\__, \___|
#                       |___/
#

include( ./cmake/CMakeLists-cpack.txt )

ebertolazzi_final_messages()
