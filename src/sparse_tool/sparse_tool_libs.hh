/*

  \file     sparse_tools_libs.hh
  \date     2011, July 21
  \version  1.1
  \note     based on SparseLib 11.0 (C 1999 -- 2008)

  \author Enrico Bertolazzi

  \par Affiliations:
       Dipartimento di Ingegneria Industriale<br>
       Università degli Studi di Trento <br>
       enrico.bertolazzi\@unitn.it

*/

#pragma once
#ifndef SPARSETOOL_LIBS_HH
#define SPARSETOOL_LIBS_HH

// automatic library inclusion
#if defined(UTILS_OS_WINDOWS) && !defined(MINGW)
  #if defined(_DEBUG) || defined(DEBUG)
    #ifdef UTILS_ARCH64
      #pragma comment(lib, "superlu_win_x64_static_debug.lib")
    #else
      #pragma comment(lib, "superlu_win_x86_static_debug.lib")
    #endif
  #else
    #ifdef UTILS_ARCH64
      #pragma comment(lib, "superlu_win_x64_static.lib")
    #else
      #pragma comment(lib, "superlu_win_x86_static.lib")
    #endif
  #endif
#endif

#endif
