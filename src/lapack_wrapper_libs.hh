/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: lapack_wrapper_libs.hh
///

#pragma once

#ifndef LAPACK_WRAPPER_LIBS_dot_HH
#define LAPACK_WRAPPER_LIBS_dot_HH

#include "lapack_wrapper_config.hh"

#if defined(LAPACK_WRAPPER_OS_WINDOWS) && defined(_MSC_VER)
  // automatically include windows libs in visual studio
  #pragma comment(lib, "shell32.lib")
  #pragma comment(lib, "kernel32.lib")
  #pragma comment(lib, "advapi32.lib")
  #pragma comment(lib, "ws2_32.lib")
  #pragma comment(lib, "IPHLPAPI.lib")
  #pragma comment(lib, "Shlwapi.lib")
  #pragma comment(lib, "iphlpapi.lib")
#endif

  // select LAPACK dll/lib for windows
  #if _MSC_VER >= 1900
    #pragma comment(lib, "legacy_stdio_definitions.lib")
  #endif

  #ifndef MINGW
    #if defined(LAPACK_WRAPPER_USE_LAPACK2)
      #if defined(_DEBUG) || defined(DEBUG)
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "cbia.lib.blas.dyn.dbg.x64.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.dbg.x64.12.lib")
        #else
          #pragma comment(lib, "cbia.lib.blas.dyn.dbg.x86.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.dbg.x86.12.lib")
        #endif
      #else
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "cbia.lib.blas.dyn.rel.x64.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.rel.x64.12.lib")
        #else
          #pragma comment(lib, "cbia.lib.blas.dyn.rel.x86.12.lib")
          #pragma comment(lib, "cbia.lib.lapack.dyn.rel.x86.12.lib")
        #endif
      #endif
    #elif defined(LAPACK_WRAPPER_USE_LAPACK)
      #if defined(_DEBUG) || defined(DEBUG)
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "blas_win64_MTd.lib")
          #pragma comment(lib, "lapack_win64_MTd.lib")
        #else
          #pragma comment(lib, "blas_win32_MTd.lib")
          #pragma comment(lib, "lapack_win32_MTd.lib")
        #endif
      #else
        #ifdef LAPACK_WRAPPER_ARCH64
          #pragma comment(lib, "blas_win64_MT.lib")
          #pragma comment(lib, "lapack_win64_MT.lib")
        #else
          #pragma comment(lib, "blas_win32_MT.lib")
          #pragma comment(lib, "lapack_win32_MT.lib")
        #endif
      #endif
    #elif defined(LAPACK_WRAPPER_USE_OPENBLAS) || \
          defined(LAPACK_WRAPPER_USE_BLASFEO)
      // no debug version
      #ifdef LAPACK_WRAPPER_ARCH64
        #pragma comment(lib, "libopenblas_x64.lib")
      #else
        #pragma comment(lib, "libopenblas_x86.lib")
      #endif
    #elif defined(LAPACK_WRAPPER_USE_MKL)
      #ifdef LAPACK_WRAPPER_ARCH64
        #pragma comment(lib, "mkl_intel_lp64_dll.lib")
      #else
        #pragma comment(lib, "mkl_intel_c_dll.lib")
      #endif
      #pragma comment(lib, "mkl_sequential_dll.lib")
      #pragma comment(lib, "mkl_core_dll.lib")
    #else
      #error "Only standard Lapack, Openblas, and MKL are supported for WINDOWS!"
    #endif
  #endif
#else
  #include <cstdint>
#endif

#endif

///
/// eof: lapack_wrapper_libs.hh
///
