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

cmake_minimum_required(VERSION 3.5)

IF ( (CMAKE_BUILD_TYPE EQ "Debug") OR (CMAKE_BUILD_TYPE EQ "DEBUG") )
  IF(APPLE)
    SET( TARGET     lapack_wrapper_osx_debug )
    SET( TARGETS    lapack_wrapper_osx_static_debug )
    SET( TARGETHSL  HSL_osx_debug )
    SET( TARGETHSLS HSL_osx_static_debug )
  ELSEIF( UNIX )
    SET( TARGET     lapack_wrapper_linux_debug )
    SET( TARGETS    lapack_wrapper_linux_static_debug )
    SET( TARGETHSL  HSL_linux_debug )
    SET( TARGETHSLS HSL_linux_static_debug )
  ELSE()
    SET( TARGET     lapack_wrapper_win_${BITS}_debug )
    SET( TARGETS    lapack_wrapper_win_${BITS}_static_debug )
    SET( TARGETHSL  HSL_win_${BITS} )
    SET( TARGETHSLS HSL_win_${BITS}_static )
    SET( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
  ENDIF()
else()
  IF(APPLE)
    SET( TARGET     lapack_wrapper_osx )
    SET( TARGETS    lapack_wrapper_osx_static )
    SET( TARGETHSL  HSL_osx )
    SET( TARGETHSLS HSL_osx_static )
  ELSEIF( UNIX )
    SET( TARGET     lapack_wrapper_linux )
    SET( TARGETS    lapack_wrapper_linux_static )
    SET( TARGETHSL  HSL_linux )
    SET( TARGETHSLS HSL_linux_static )
  ELSE()
    SET( TARGET     lapack_wrapper_win_${BITS} )
    SET( TARGETS    lapack_wrapper_win_${BITS}_static )
    SET( TARGETHSL  HSL_win_${BITS} )
    SET( TARGETHSLS HSL_win_${BITS}_static )
    SET( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
  ENDIF()
endif()

PROJECT( ${TARGET} CXX C )
enable_testing()

SET(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/lib")

FILE( GLOB S ./CMakeLists-cflags.txt )
IF(EXISTS ${S})
  INCLUDE( ${S} )
ELSE()
  SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
  SET( CMAKE_VERBOSE_MAKEFILE            TRUE )
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    SET( CMAKE_CXX_FLAGS "-std=c++11 " )
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -Wall -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -funroll-loops -O3 -g0 " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -Wall -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -O0 -g3 -DMECHATRONIX_DEBUG " )
  ELSEIF( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    SET( CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ " )
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -funroll-loops -O3 -g0 -fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -O0 -gfull -DMECHATRONIX_DEBUG -fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation " )
  ELSEIF( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    SET( CMAKE_CXX_FLAGS_RELEASE "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_WINDOWS /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC /O2 /MD " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_WINDOWS /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC /Od /Ob0 /MDd /Zi /RTC1 /D_DEBUG " )
  ELSE()
    MESSAGE( FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
  ENDIF()
ENDIF()

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
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})
FILE( GLOB S ./src/*/*/*.h* )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

IF( APPLE )
  IF ( LAPACK_WRAPPER_USE_OPENBLAS )
    execute_process( COMMAND gfortran -print-libgcc-file-name OUTPUT_VARIABLE GFORTRANDIR )
    get_filename_component( GFORTRANDIR1 ${GFORTRANDIR} DIRECTORY )
    MESSAGE( STATUS ${GFORTRANDIR1} )
    FIND_LIBRARY( GFORTRAN libgfortran.dylib HINTS ${GFORTRANDIR1}/../../.. )
    FIND_LIBRARY( OPENBLAS libopenblas.dylib HINTS /usr/local/opt/openblas/lib )
    SET( lapackblas_libraries ${OPENBLAS} ${GFORTRAN} )
    INCLUDE_DIRECTORIES( /usr/local/opt/openblas/include )
  ELSEIF( LAPACK_WRAPPER_USE_MKL )
    FIND_LIBRARY( CORE libmkl_core.dylib HINTS /opt/intel/mkl/lib/ )
    FIND_LIBRARY( INTEL libmkl_intel.dylib HINTS /opt/intel/mkl/lib/ )
    FIND_LIBRARY( SEQUENTIAL libmkl_sequential.dylib HINTS /opt/intel/mkl/lib/ )
    SET( lapackblas_libraries ${CORE} ${INTEL} ${SEQUENTIAL} )
    INCLUDE_DIRECTORIES( /opt/intel/mkl/include/ )
  ELSE()
    FIND_PACKAGE(BLAS)
    FIND_PACKAGE(LAPACK)
    IF ( LAPACK_FOUND AND BLAS_FOUND )
      SET( lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
    ENDIF()
  ENDIF()
ELSEIF( UNIX )
  FIND_PACKAGE(BLAS)
  FIND_PACKAGE(LAPACK)
  IF ( LAPACK_FOUND AND BLAS_FOUND )
    SET( lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
  ENDIF()
ELSE()
  IF ( BITS MATCHES "x64")
    FIND_LIBRARY(
      lapackblas_libraries libopenblas_x64.lib
      HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    )
  ELSE()
    FIND_LIBRARY(
      lapackblas_libraries libopenblas.lib
      HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    )
  ENDIF()
ENDIF()

add_library( ${TARGET}     SHARED ${SOURCES} ${HEADERS} )
add_library( ${TARGETS}    STATIC ${SOURCES} ${HEADERS} )
add_library( ${TARGETHSL}  SHARED src/HSL/hsl_fake.cc )
add_library( ${TARGETHSLS} STATIC src/HSL/hsl_fake.cc )
TARGET_LINK_LIBRARIES( ${TARGET} ${TARGETHSL} ${lapackblas_libraries} )

SET(
  EXELISTCPP
  test1-small-factorization
  test2-Timing
  test3-BandedMatrix
  test4-BFGS
  test5-BLOCKTRID
  test6-EIGS
  #test7-SparseTool
  #test6-SparseToolComplex
)

MESSAGE( STATUS "YEAR = ${YEAR}" )
MESSAGE( STATUS "BITS = ${BITS}" )

SET(LLPATH ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )

INCLUDE_DIRECTORIES( src ${LLPATH}/include )

IF( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
  INCLUDE_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include" )
  IF ( ${BITS} STREQUAL "x64" )
    set(WINXXBITS Win64)
    LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64" )
  ELSE()
    set(WINXXBITS Win32)
    LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/ia32" )
  ENDIF()
ENDIF()

IF(APPLE)
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
ELSEIF( UNIX )
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

  #FIND_PACKAGE(BLAS)
  #FIND_PACKAGE(LAPACK)
  #if(LAPACK_FOUND AND BLAS_FOUND)
  #  SET(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
  #endif()

  ADD_CUSTOM_TARGET( all_tests ALL )

  SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  FILE(MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FOREACH( S ${EXELISTCPP} )
  	ADD_EXECUTABLE( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${lapackblas_libraries} )
  	ADD_TEST( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    ADD_DEPENDENCIES( all_tests ${S} )
  ENDFOREACH( S $(EXELIST) )

  FOREACH( S ${EXELISTC} )
  	ADD_EXECUTABLE( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.c ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${lapackblas_libraries} )
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

FILE( GLOB HDR_HSL  ./src/HSL/*.h* )
FILE( GLOB HDR_LW   ./src/lapack_wrapper/*.h* )
FILE( GLOB HDR_LW1  ./src/lapack_wrapper/code/*.h* )
FILE( GLOB HDR_LW2  ./src/lapack_wrapper/code++/*.h* )
FILE( GLOB HDR_ST   ./src/sparse_tool/*.h* )
FILE( GLOB HDR_ST1  ./src/sparse_tool/interfaces/*.h* )
FILE( GLOB HDR_ST2  ./src/sparse_tool/iterative/*.h* )
FILE( GLOB HDR_ST3  ./src/sparse_tool/preconditioner/*.h* )
FILE( GLOB HDR_ZS   ./src/zstream/*.h* )

INSTALL(
  TARGETS ${TARGETHSL} ${TARGETHSLS} ${TARGET} ${TARGETS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

INSTALL(FILES ${HDR_HSL} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/HSL )
INSTALL(FILES ${HDR_LW}  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper )
INSTALL(FILES ${HDR_LW1} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/code )
INSTALL(FILES ${HDR_LW2} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/code++ )
INSTALL(FILES ${HDR_ST}  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool )
INSTALL(FILES ${HDR_ST1} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/interfaces )
INSTALL(FILES ${HDR_ST2} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/iterative )
INSTALL(FILES ${HDR_ST3} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/preconditioner )
INSTALL(FILES ${HDR_ZS}  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/zstream )

# do the copying
foreach( file_glob ${HDR_HSL})
  get_filename_component( file_loc "${file_glob}" NAME)
  add_custom_command(
    TARGET ${TARGET}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${file_glob}
            ${CMAKE_INSTALL_PREFIX}/include/HSL/${file_loc}
)
endforeach( file_glob )

foreach( file_glob ${HDR_LW})
  get_filename_component( file_loc "${file_glob}" NAME)
  add_custom_command(
    TARGET ${TARGET}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${file_glob}
            ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/${file_loc}
)
endforeach( file_glob )

foreach( file_glob ${HDR_LW})
  get_filename_component( file_loc "${file_glob}" NAME)
  add_custom_command(
    TARGET ${TARGET}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${file_glob}
            ${CMAKE_INSTALL_PREFIX}/include/lapack_wrapper/${file_loc}
)
endforeach( file_glob )

foreach( file_glob ${HDR_ST})
  get_filename_component( file_loc "${file_glob}" NAME)
  add_custom_command(
    TARGET ${TARGET}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${file_glob}
            ${CMAKE_INSTALL_PREFIX}/include/sparse_tool/${file_loc}
)
endforeach( file_glob )

foreach( file_glob ${HDR_ZS})
  get_filename_component( file_loc "${file_glob}" NAME)
  add_custom_command(
    TARGET ${TARGET}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${file_glob}
            ${CMAKE_INSTALL_PREFIX}/include/zstream/${file_loc}
)
endforeach( file_glob )

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
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()

MESSAGE( STATUS "lapackblas_libraries        = ${lapackblas_libraries}" )
