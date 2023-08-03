############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.14)

project( lapack_wrapper )
include(./cmake_utils/CMakeLists-common.txt)

project(
  ${PROJECT_NAME}
  VERSION ${UTILS_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/LapackWrapper/"
)

include(./cmake_utils/CMakeLists-cflags.txt)
include(./cmake_utils/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
utils_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( UTILS_BUILD_SHARED )
  utils_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

if ( APPLE )
  set( CMAKE_OSX_ARCHITECTURES arm64;x86_64 )
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
  file( RELATIVE_PATH RF ${SOURCE_DIR} ${F} )
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
  file( RELATIVE_PATH RF ${SOURCE_DIR} ${F} )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

# extra include directory for compilation
include_directories( ./lib3rd/include )

#       _                           _                 _
#    __| | ___ _ __   ___ _ __   __| | ___ _ __   ___(_) ___  ___
#   / _` |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \ / __| |/ _ \/ __|
#  | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (__| |  __/\__ \
#   \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___|_|\___||___/
#             |_|
#
set(DEPEND_TARGETS)

find_library(
  SUPERLU superlu${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  NO_DEFAULT_PATH
  REQUIRED
)

# DEPENDENCIES (TRY TO USE ALREADY COMPILED SUBMODULE IN LARGE PROJECT)
set(UPDIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
find_library(
  UTILS Utils${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${UPDIR}/Utils/lib/lib
  NO_DEFAULT_PATH
)
if( UTILS )
  message(STATUS "Found Utils${ARTIFACTS_STATIC_SUFFIX}")
  utils_copy_directory(
    "${UPDIR}/Utils/lib"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd"
  )
else()
  message(STATUS "NOT Found Utils${ARTIFACTS_STATIC_SUFFIX} use submodule")
  add_subdirectory(./submodules/Utils)
  include_directories(./submodules/Utils/src)
  set( UTILS utils_Utils_Static )
  set( DEPEND_TARGETS ${DEPEND_TARGETS} ${UTILS} )
  #add_custom_command(
  #  TARGET ${UTILS}
  #  PRE_BUILD
  #  COMMAND "rake install"
  #  WORKING_DIRECTORY "./submodules/Utils/ThirdParties"
  #)
endif()

# definisce lapackblas_libraries
include( ./cmake_utils/CMakeLists-blas-lapack.txt )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#
# define target
utils_setup_target(
  "${PROJECT_NAME}"
  TRUE # build static
  "${SOURCE_DIR}"
  "${SOURCES}"
  "${SOURCE_DIR}"
  "${PUBLIC_HEADERS}"
)
utils_setup_target(
  "HSL"
  TRUE # build static
  "${SOURCE_DIR}"
  "HSL/hsl_fake.cc" # must be relative path
  "${SOURCE_DIR}"
  "${PUBLIC_HEADERS}"
)
if ( UTILS_BUILD_SHARED )
  utils_setup_target(
    "${PROJECT_NAME}"
    FALSE # build dynamic
    "${SOURCE_DIR}"
    "${SOURCES}"
    "${SOURCE_DIR}"
    "${PUBLIC_HEADERS}"
  )
  utils_setup_target(
    "HSL"
    FALSE # build dynamic
    "${SOURCE_DIR}"
    "HSL/hsl_fake.cc" # must be relative path
    "${SOURCE_DIR}"
    "${PUBLIC_HEADERS}"
  )
  target_link_libraries(
    ${UTILS_NAMESPACE}_${PROJECT_NAME}
    ${SUPERLU} ${lapackblas_libraries} ${UTILS_NAMESPACE}_HSL_Static ${UTILS}
  )
endif()

if ( DEPEND_TARGETS )
  add_dependencies( ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static ${DEPEND_TARGETS} )
endif()

#   _____       _
#  |_   _|__ __| |_ ___
#    | |/ -_|_-<  _(_-<
#    |_|\___/__/\__/__/
#
if ( UTILS_ENABLE_TESTS )

  enable_testing()

  if( APPLE OR UNIX OR MINGW )
    find_package( ZLIB REQUIRED )
  else()
    if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
      find_library(
        ZLIB_LIBRARIES zlibstaticd
        HINTS ./submodules/Utils/ThirdParties/zlib-ng/lib/win/lib
              ${UPDIR}/Utils/lib3rd/lib ${UPDIR}/Utils/lib3rd/dll
      )
    else()
      find_library(
        ZLIB_LIBRARIES zlibstatic
        HINTS ./submodules/Utils/ThirdParties/zlib-ng/lib/win/lib
               ${UPDIR}/Utils/lib3rd/lib ${UPDIR}/Utils/lib3rd/dll
      )
    endif()
    set( ZLIB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParties/zlib-ng/win/include " "${UPDIR}/Utils/lib3rd/include ")
  endif()

  set( EXELISTBIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin" )
  set( EXELISTCPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src_tests" )

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

  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )

  set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  set(
    LIBS
    ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static
    ${UTILS_NAMESPACE}_HSL_Static
    ${UTILS}
    ${SUPERLU}
    ${ZLIB_LIBRARIES}
    ${lapackblas_libraries}
    ${CMAKE_DL_LIBS}
  )
  foreach( S ${EXELISTCPP} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc )
    target_link_libraries( ${S} ${LIBS} )
    set_target_properties( ${S} PROPERTIES SUFFIX ".exe" )
    add_test( NAME "${S}" COMMAND ./bin/${S}.exe WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( "${PROJECT_NAME}_all_tests" ${S} )
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

#   ___         _
#  | _ \__ _ __| |____ _ __ _ ___
#  |  _/ _` / _| / / _` / _` / -_)
#  |_| \__,_\__|_\_\__,_\__, \___|
#                       |___/
#

# do not work
include( ./cmake_utils/CMakeLists-cpack.txt )

utils_final_messages()
