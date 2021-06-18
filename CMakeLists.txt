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

CMAKE_MINIMUM_REQUIRED( VERSION 2.8.12 )
SET( CMAKE_VERBOSE_MAKEFILE OFF )
SET( CMAKE_INSTALL_MESSAGE NEVER )

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

PROJECT( lapack_wrapper CXX C )
ENABLE_TESTING()

IF( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  SET( DEBUG_SUFFIX "_debug" )
ELSE()
  SET( DEBUG_SUFFIX )
ENDIF()

IF( APPLE )
ELSEIF( UNIX OR LINUX )
ELSE()
  SET( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
ENDIF()

IF( APPLE )
  SET( TARGETS    lapack_wrapper_osx_static${DEBUG_SUFFIX} )
  SET( TARGETHSLS HSL_osx_static${DEBUG_SUFFIX} )
ELSEIF( UNIX )
  SET( TARGETS    lapack_wrapper_linux_static${DEBUG_SUFFIX} )
  SET( TARGETHSLS HSL_linux_static${DEBUG_SUFFIX} )
ELSE()
  SET( TARGETS    lapack_wrapper_win_${BITS}_static${DEBUG_SUFFIX} )
  SET( TARGETHSLS HSL_win_${BITS}_static${DEBUG_SUFFIX} )
ENDIF()

IF( APPLE )
  FIND_LIBRARY(
    UTILS Utils_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ELSEIF ( UNIX OR LINUX )
  FIND_LIBRARY(
    UTILS Utils_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ELSE()
  FIND_LIBRARY(
    UTILS Utils_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ENDIF()

IF ( BUILD_SHARED )
  IF(APPLE)
    SET( TARGET    lapack_wrapper_osx${DEBUG_SUFFIX} )
    SET( TARGETHSL HSL_osx${DEBUG_SUFFIX} )
  ELSEIF( UNIX )
    SET( TARGET    lapack_wrapper_linux${DEBUG_SUFFIX} )
    SET( TARGETHSL HSL_linux${DEBUG_SUFFIX} )
  ELSE()
    SET( TARGET    lapack_wrapper_win_${BITS}${DEBUG_SUFFIX} )
    SET( TARGETHSL HSL_win_${BITS}${DEBUG_SUFFIX} )
  ENDIF()
ENDIF()

SET( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" )

INCLUDE( ./CMakeLists-cflags.txt )

MESSAGE( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")

SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./src/*/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  IF ( F MATCHES "hsl_fake.cc")
    MESSAGE( STATUS "skip ${F}" )
  ELSE()
    MESSAGE( STATUS "add ${F}" )
    LIST( APPEND SOURCES ${RF} )
  ENDIF()
ENDFOREACH (F ${S})

SET( HEADERS )
FILE( GLOB S ./src/*/*.h* )
FOREACH( F ${S} )
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})
FILE( GLOB S ./src/*/*/*.h* )
FOREACH( F ${S} )
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

# definisce lapackblas_libraries
INCLUDE( ./CMakeLists-blas-lapack.txt )

ADD_LIBRARY( ${TARGETS}    STATIC ${SOURCES} ${HEADERS} )
ADD_LIBRARY( ${TARGETHSLS} STATIC src/HSL/hsl_fake.cc )

IF ( BUILD_SHARED )
  ADD_LIBRARY( ${TARGET}    SHARED ${SOURCES} ${HEADERS} )
  ADD_LIBRARY( ${TARGETHSL} SHARED src/HSL/hsl_fake.cc )
  TARGET_LINK_LIBRARIES( ${TARGET} ${TARGETHSL} ${lapackblas_libraries} ${UTILS} )
ENDIF()

MESSAGE( STATUS "YEAR = ${YEAR}" )
MESSAGE( STATUS "BITS = ${BITS}" )

SET( LLPATH ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )

INCLUDE_DIRECTORIES( src )
INCLUDE_DIRECTORIES( SYSTEM ${LLPATH}/include )

IF( APPLE OR UNIX )
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
ELSE()
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin/$(BITS)
  )
ENDIF()

IF( BUILD_EXECUTABLE )

  IF(APPLE OR UNIX)
    FIND_PACKAGE( ZLIB REQUIRED )
    INCLUDE_DIRECTORIES( ${ZLIB_INCLUDE_DIRS} )
  ELSE()
    IF ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
      FIND_LIBRARY( ZLIB_LIBRARIES zlibstaticd HINTS ${LLPATH}/lib ${LLPATH}/dll )
    ELSE()
      FIND_LIBRARY( ZLIB_LIBRARIES zlibstatic HINTS ${LLPATH}/lib ${LLPATH}/dll )
    ENDIF()
  ENDIF()

  SET(
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

  ADD_CUSTOM_TARGET( all_tests ALL )

  SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FILE( MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FOREACH( S ${EXELISTCPP} )
  	ADD_EXECUTABLE( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc ${HEADERS} )
    TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${lapackblas_libraries} ${UTILS} ${SUPERLU} ${ZLIB_LIBRARIES} )
    IF ( UNIX AND NOT APPLE )
      TARGET_LINK_LIBRARIES( ${S} ${DL} )
    ENDIF()
  	ADD_TEST( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    ADD_DEPENDENCIES( all_tests ${S} )
  ENDFOREACH( S $(EXELIST) )

  FOREACH( S ${EXELISTC} )
  	ADD_EXECUTABLE( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.c ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${lapackblas_libraries} ${UTILS} ${SUPERLU} ${ZLIB_LIBRARIES} )
    IF ( UNIX AND NOT APPLE )
      TARGET_LINK_LIBRARIES( ${S} ${DL} )
    ENDIF()
  	ADD_TEST( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    ADD_DEPENDENCIES( all_tests ${S} )
  ENDFOREACH( S $(EXELIST) )

  ADD_CUSTOM_COMMAND(
    TARGET all_tests
    COMMENT "Run tests"
    POST_BUILD COMMAND ctest ARGS --output-on-failure
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

ENDIF()

FILE( GLOB HDR       ./src/*.h* )
FILE( GLOB HDR_HSL   ./src/HSL/*.h* )
FILE( GLOB HDR_LW    ./src/lapack_wrapper/*.h* )
FILE( GLOB HDR_LW1   ./src/lapack_wrapper/code/*.h* )
FILE( GLOB HDR_LW2   ./src/lapack_wrapper/code++/*.h* )
FILE( GLOB HDR_ST    ./src/sparse_tool/*.h* )
FILE( GLOB HDR_ST1   ./src/sparse_tool/interfaces/*.h* )
FILE( GLOB HDR_ST2   ./src/sparse_tool/iterative/*.h* )
FILE( GLOB HDR_ST3   ./src/sparse_tool/preconditioner/*.h* )
FILE( GLOB HDR_OBx64 ./src/openblas/x64/*.h* )
FILE( GLOB HDR_OBx86 ./src/openblas/x66/*.h* )

IF ( BUILD_SHARED )
  INSTALL(
    TARGETS ${TARGETHSL} ${TARGETHSLS} ${TARGETS} ${TARGET}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ELSE()
  INSTALL(
    TARGETS ${TARGETHSLS} ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ENDIF()

INSTALL( FILES ${HDR}       DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
INSTALL( FILES ${HDR_HSL}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/HSL )
INSTALL( FILES ${HDR_LW}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper )
INSTALL( FILES ${HDR_LW1}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/code )
INSTALL( FILES ${HDR_LW2}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/code++ )
INSTALL( FILES ${HDR_ST}    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool )
INSTALL( FILES ${HDR_ST1}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/interfaces )
INSTALL( FILES ${HDR_ST2}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/iterative )
INSTALL( FILES ${HDR_ST3}   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/preconditioner )
INSTALL( FILES ${HDR_OBx64} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/openblas/x64 )
INSTALL( FILES ${HDR_OBx86} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/openblas/x86 )
INSTALL( DIRECTORY ./src/Eigen DESTINATION ${CMAKE_INSTALL_PREFIX}/include )

MESSAGE( STATUS "Using ${SSE_FLAGS} extensions")
MESSAGE( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
MESSAGE( STATUS "EXTRA_LIBS                  = ${EXTRA_LIBS}" )
IF( APPLE )
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()

MESSAGE( STATUS "lapackblas_libraries        = ${lapackblas_libraries}" )
MESSAGE( STATUS "Utils                       = ${UTILS}" )
