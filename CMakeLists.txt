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

cmake_minimum_required( VERSION 3.0.0 )
set( CMAKE_VERBOSE_MAKEFILE OFF )
set( CMAKE_INSTALL_MESSAGE NEVER )

# get version from the git tag
execute_process(COMMAND git describe --tags --abbrev=0 OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "\n$" "" VERSION "${VERSION}")

project(
  lapack_wrapper
  VERSION ${VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
)

enable_testing()

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  set( DEBUG_SUFFIX "_debug" )
else()
  set( DEBUG_SUFFIX )
endif()

if( APPLE )
elseif( UNIX OR LINUX )
else()
  set( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
endif()

include(./CMakeLists-cflags.txt)
include(./CMakeUtilities.txt)

# Evaluating the suffix to append to target built by this CMake file
ebertolazzi_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if(BUILD_SHARED)
  ebertolazzi_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

set( TARGETS    lapack_wrapper${ARTIFACTS_STATIC_SUFFIX} )
set( TARGETHSLS HSL${ARTIFACTS_STATIC_SUFFIX} )

find_library(
  SUPERLU superlu${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
)

if ( BUILD_SHARED )
  set( TARGET    lapack_wrapper${ARTIFACTS_DYNAMIC_SUFFIX} )
  set( TARGETHSL HSL${ARTIFACTS_DYNAMIC_SUFFIX} )
endif()

# DEPENDENCIES (TRY TO USE ALREADY COMPILED SUBMODULE IN LARGE PROJECT)
set(UPDIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../..)
set(DEPEND_TARGETS)
if( EXISTS "${UPDIR}/Utils/lib" )
  include_directories(${UPDIR}/Utils/lib/include)
  find_library( UTILS Utils${ARTIFACTS_STATIC_SUFFIX} HINTS ${UPDIR}/Utils/lib/lib )
else()
  add_subdirectory(./submodules/Utils)
  set( UTILS UtilsStatic )
  set( DEPEND_TARGETS ${DEPEND_TARGETS} UtilsStatic )
endif()

set( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" )

# definisce lapackblas_libraries
include( ./CMakeLists-blas-lapack.txt )

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")

set( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
set( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

set( SOURCES )
file( GLOB S ./src/*/*.cc )
foreach (F ${S})
file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  if ( F MATCHES "hsl_fake.cc")
    message( STATUS "skip ${F}" )
  else()
    message( STATUS "add ${F}" )
    list( APPEND SOURCES ${RF} )
  endif()
endforeach (F ${S})

set( HEADERS )
file( GLOB S ./src/*/*.h* )
foreach( F ${S} )
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
endforeach (F ${S})
file( GLOB S ./src/*/*/*.h* )
foreach( F ${S} )
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
endforeach (F ${S})

add_library( ${TARGETS}    STATIC ${SOURCES} ${HEADERS} )
add_library( ${TARGETHSLS} STATIC src/HSL/hsl_fake.cc )

if ( BUILD_SHARED )
  add_library( ${TARGET}    SHARED ${SOURCES} ${HEADERS} )
  add_library( ${TARGETHSL} SHARED src/HSL/hsl_fake.cc )
  target_link_libraries( ${TARGET} ${TARGETHSL} ${lapackblas_libraries} ${UTILS} )
endif()

add_dependencies( ${TARGETS}    ${DEPEND_TARGETS} )
add_dependencies( ${TARGETHSLS} ${DEPEND_TARGETS} )

set( LLPATH ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )

include_directories( src )
include_directories( SYSTEM ${LLPATH}/include )

if ( APPLE OR UNIX )
  link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
else()
  link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin/$(BITS)
  )
endif()

if ( BUILD_EXECUTABLE )

  if (APPLE OR UNIX)
    find_package( ZLIB REQUIRED )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
  else()
    if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
      find_library( ZLIB_LIBRARIES zlibstaticd HINTS ${LLPATH}/lib ${LLPATH}/dll )
    else()
      find_library( ZLIB_LIBRARIES zlibstatic HINTS ${LLPATH}/lib ${LLPATH}/dll )
    endif()
  endif()

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
  ##test9-RRQR

  add_custom_target( all_tests ALL )

  set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  file( MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  foreach( S ${EXELISTCPP} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc ${HEADERS} )
    target_link_libraries( ${S} ${TARGETS} ${lapackblas_libraries} ${UTILS} ${SUPERLU} ${ZLIB_LIBRARIES} ${CMAKE_DL_LIBS} )
    add_test( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( all_tests ${S} )
  endforeach( S $(EXELIST) )

  foreach( S ${EXELISTC} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.c ${HEADERS} )
    target_link_libraries( ${S} ${TARGETS} ${lapackblas_libraries} ${UTILS} ${SUPERLU} ${ZLIB_LIBRARIES} ${CMAKE_DL_LIBS} )
    add_test( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( all_tests ${S} )
  endforeach( S $(EXELIST) )

  add_custom_command(
    TARGET all_tests
    COMMENT "Run tests"
    POST_BUILD COMMAND ctest ARGS --output-on-failure
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

endif()

file( GLOB HDR       ./src/*.h* )
file( GLOB HDR_HSL   ./src/HSL/*.h* )
file( GLOB HDR_LW    ./src/lapack_wrapper/*.h* )
file( GLOB HDR_LW1   ./src/lapack_wrapper/code/*.h* )
file( GLOB HDR_LW2   ./src/lapack_wrapper/code++/*.h* )
file( GLOB HDR_ST    ./src/sparse_tool/*.h* )
file( GLOB HDR_ST1   ./src/sparse_tool/interfaces/*.h* )
file( GLOB HDR_ST2   ./src/sparse_tool/iterative/*.h* )
file( GLOB HDR_ST3   ./src/sparse_tool/preconditioner/*.h* )
file( GLOB HDR_OBx64 ./src/openblas/x64/*.h* )
file( GLOB HDR_OBx86 ./src/openblas/x66/*.h* )

if ( BUILD_SHARED )
  install(
    TARGETS ${TARGETHSL} ${TARGETHSLS} ${TARGETS} ${TARGET}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
else()
  install(
    TARGETS ${TARGETHSLS} ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
endif()

install( FILES ${HDR}       DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
install( FILES ${HDR_HSL}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/HSL )
install( FILES ${HDR_LW}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper )
install( FILES ${HDR_LW1}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/code )
install( FILES ${HDR_LW2}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/code++ )
install( FILES ${HDR_ST}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool )
install( FILES ${HDR_ST1}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/interfaces )
install( FILES ${HDR_ST2}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/iterative )
install( FILES ${HDR_ST3}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/preconditioner )
install( FILES ${HDR_OBx64} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/openblas/x64 )
install( FILES ${HDR_OBx86} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/openblas/x86 )

message( STATUS "Using ${SSE_FLAGS} extensions")
message( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
message( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
message( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
message( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
message( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
message( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
message( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
message( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
message( STATUS "EXTRA_LIBS                  = ${EXTRA_LIBS}" )

if( APPLE )
  message( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  message( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  message( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
endif()

message( STATUS "lapackblas_libraries        = ${lapackblas_libraries}" )
message( STATUS "ZLIB_LIBRARIES              = ${ZLIB_LIBRARIES}" )
message( STATUS "ZLIB_INCLUDE_DIRS           = ${ZLIB_INCLUDE_DIRS}" )
message( STATUS "Utils                       = ${UTILS}" )
