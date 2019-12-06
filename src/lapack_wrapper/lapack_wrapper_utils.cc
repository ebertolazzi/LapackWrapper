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

#include "../lapack_wrapper_config.hh"

// workaround for windows
#ifdef _MSC_VER
  #ifdef max
    #undef max
  #endif
  #ifdef min
    #undef min
  #endif
#endif

#include <algorithm>
#include <string>
#include <sstream>

#ifndef LAPACK_WRAPPER_OS_WINDOWS
#include <execinfo.h> // for backtrace
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#endif

#include <thread>

namespace lapack_wrapper {

  /*\
  :|:  _
  :|: | |__  __ _ ___ ___ _ _  __ _ _ __  ___
  :|: | '_ \/ _` (_-</ -_) ' \/ _` | '  \/ -_)
  :|: |_.__/\__,_/__/\___|_||_\__,_|_|_|_\___|
  :|:
  \*/

  #ifdef LAPACK_WRAPPER_OS_WINDOWS
    std::string
    basename( char const path[] ) {
      static char drive[100];
      static char dir[1024];
      static char fname[256];
      static char ext[128];
      errno_t e = _splitpath_s(
        path,
        drive, 100,
        dir,   1024,
        fname, 256,
        ext,   128
      );
      LW_ASSERT0( e == 0, "lapack_wrapper, basename failed!\n" );
      return fname;
    }
  #else
    std::string
    basename( char const path[] ) {

      if ( path[0] == '\0' ) return std::string("");

      std::string filename(path);

      size_t len   = filename.length();
      size_t index = filename.find_last_of("/\\");

      if ( index == std::string::npos ) return filename;
      if ( index + 1 >= len ) {
        --len;
        index = filename.substr(0, len).find_last_of("/\\");

        if ( len   == 0 ) return filename;
        if ( index == 0 ) return filename.substr(1, len - 1);
        if ( index == std::string::npos ) return filename.substr(0, len);
        return filename.substr(index + 1, len - index - 1);
      }
      return filename.substr(index + 1, len - index);
    }
  #endif

  void
  Console::changeLevel( int new_level ) {
    LW_ASSERT(
      new_level >= -1 && new_level <= 4,
      "Console::changeLevel( new_level = {})\nnew_level must be in the range [-1,4]\n",
      new_level
    );
    this->level = new_level;
  }

  void
  Console::changeStream( ostream_type * new_p_stream  ) {
    this->p_stream = new_p_stream;
  }

  Console const &
  Console::black( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::black << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::red( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::red << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::green( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::green << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::yellow( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::yellow << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::blue( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::blue << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::magenta( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::magenta << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::cyan( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::cyan << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::gray( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream) << rang::fg::gray << msg << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::black_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::black << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::red_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::red << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::green_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::green << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::yellow_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::yellow << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::blue_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::blue << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::magenta_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::magenta << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::cyan_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::cyan << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::gray_reversed( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << rang::fg::gray << rang::style::reversed
        << msg << rang::style::reset << rang::fg::reset;
    return *this;
  }

  Console::Console( ostream_type * _p_stream, int _level )
  : p_stream(_p_stream)
  , level(_level)
  {
    this->message_style.s = rang::style::reset;
    this->message_style.f = rang::fg::reset;
    this->message_style.b = rang::bg::reset;

    this->warning_style.s = rang::style::reset;
    this->warning_style.f = rang::fg::yellow;
    this->warning_style.b = rang::bg::reset;

    this->error_style.s = rang::style::italic;
    this->error_style.f = rang::fg::red;
    this->error_style.b = rang::bg::reset;

    this->fatal_style.s = rang::style::underline;
    this->fatal_style.f = rang::fg::red;
    this->fatal_style.b = rang::bg::reset;
  }

  Console const &
  Console::semaphore(
    unsigned            rvg,
    std::string const & msg,
    int                 msg_level
  ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    static rang::fg rvg_color[3] = { rang::fg::red, rang::fg::yellow, rang::fg::green };
    if ( msg_level <= level )
      (*p_stream)
        << rang::style::reset
        << rang::bg::reset
        << rvg_color[rvg%3]
        << msg
        << rang::fg::reset;
    return *this;
  }

  Console const &
  Console::message( std::string const & msg, int msg_level ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( msg_level <= level )
      (*p_stream)
        << message_style.s
        << message_style.f
        << message_style.b
        << msg
        << rang::style::reset
        << rang::fg::reset
        << rang::bg::reset;
    return *this;
  }

  Console const &
  Console::warning( std::string const & msg ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( level >= 2 )
      (*p_stream)
        << warning_style.s
        << warning_style.f
        << warning_style.b
        << msg
        << rang::style::reset
        << rang::fg::reset
        << rang::bg::reset;
    return *this;
  }

  Console const &
  Console::error( std::string const & msg ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    if ( level >= 1 )
      (*p_stream)
        << error_style.s
        << error_style.f
        << error_style.b
        << msg
        << rang::style::reset
        << rang::fg::reset
        << rang::bg::reset;
    return *this;
  }

  Console const &
  Console::fatal( std::string const & msg ) const {
    std::lock_guard<std::mutex> lock_access(message_mutex);
    (*p_stream)
      << fatal_style.s
      << fatal_style.f
      << fatal_style.b
      << msg
      << rang::style::reset
      << rang::fg::reset
      << rang::bg::reset;
    return *this;
  }

  const char*
  Runtime_Error::what() const noexcept {
    return std::runtime_error::what();
  }

  const char*
  Runtime_TraceError::what() const noexcept {
    return std::runtime_error::what();
  }

  #ifdef LAPACK_WRAPPER_OS_WINDOWS

  void
  printTrace(
    int                 line,
    char        const   file[],
    std::string const & msg,
    ostream_type & stream
  ) {
    fmt::print( stream,
      "---------------------------------------------------------\n"
      "file: {}:{}\n{}\n"
      "---------------------------------------------------------\n",
      file, line, msg
    );
  }

  /*
    #include "StackWalker.h"
    static
    inline
    void
    printStackTrace( FILE *out = stderr ) {
      fprintf( out, "stack trace:\n" );
      StackWalker sw;
      sw.ShowCallstack();
    }
  */
  std::string
  Runtime_TraceError::grab_backtrace(
    std::string const & reason,
    char const          file[],
    int                 line
  ) const {
    return fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line );
  }

  #else

  static
  inline
  std::string
  demang( char const mangled_name[] ) {
    int status = 0 ;
    std::string retval = mangled_name;
    char * name = abi::__cxa_demangle( mangled_name, nullptr, nullptr, &status );
    if ( status == 0 ) retval = name;
    if ( name != nullptr ) std::free(name) ;
    return retval;
  }
  //! print a trace stack used in debug
  void
  printTrace(
    int                 line,
    char        const   file[],
    std::string const & reason,
    ostream_type      & stream
  ) {

    fmt::print(
      stream, "\n{}\nOn File:{}:{}\nprocess ID:{}, parent process ID:{}\nstack trace:\n",
      reason, basename(file), line, getpid(), getppid()
    );

    //  record stack trace upto 128 frames
    void *callstack[128] = {};

    // collect stack frames
    int frames = backtrace( callstack, 128);

    // get the human-readable symbols (mangled)
    char** strs = backtrace_symbols( callstack, frames );

    for ( int i = 1; i < frames; ++i) {
      #ifdef LAPACK_WRAPPER_OS_LINUX
        Dl_info dlinfo;
        if( !dladdr(callstack[i], &dlinfo) ) continue;
        fmt::print( stream, "{:2} {}\n", i, demang( dlinfo.dli_sname ) );
      #else
        char functionSymbol[1024] = {};
        char moduleName[1024]     = {};
        int  offset               = 0;
        char addr[48]             = {};
        // split the string, take out chunks out of stack trace
        // we are primarily interested in module, function and address
        sscanf(
          strs[i], "%*s %s %s %s %*s %d",
          moduleName, addr, functionSymbol, &offset
        );
        //  if this is a C++ library, symbol will be demangled
        //  on success function returns 0
        //
        fmt::print( stream, "{:2} {:30}  [{}] {} + {}\n",
          i, moduleName, addr, demang( functionSymbol ), offset
        );
      #endif
    }
    free(strs);
  }

  std::string
  Runtime_TraceError::grab_backtrace(
    std::string const & reason,
    char const          file[],
    int                 line
  ) const {
    std::ostringstream ost;
    printTrace( line, file, reason, ost );
    return ost.str();
  }
  #endif


  //============================================================================
  /*    __                       _ _   _       _   _
  //   / _| ___  _   _ _ __   __| | \ | | __ _| \ | |
  //  | |_ / _ \| | | | '_ \ / _` |  \| |/ _` |  \| |
  //  |  _| (_) | |_| | | | | (_| | |\  | (_| | |\  |
  //  |_|  \___/ \__,_|_| |_|\__,_|_| \_|\__,_|_| \_|
  */
  //! check if the vector `pv` os size `DIM` contains only regular floats
  bool
  foundNaN( double const pv[], int DIM ) {
    for ( int i = 0; i < DIM; ++i )
      if ( !isRegular(pv[i]) )
        return true;
    return false;
  }

  bool
  foundNaN( float const pv[], int DIM ) {
    for ( int i = 0; i < DIM; ++i )
      if ( !isRegular(pv[i]) )
        return true;
    return false;
  }

  /*       _               _    _   _       _   _
  //   ___| |__   ___  ___| | _| \ | | __ _| \ | |
  //  / __| '_ \ / _ \/ __| |/ /  \| |/ _` |  \| |
  // | (__| | | |  __/ (__|   <| |\  | (_| | |\  |
  //  \___|_| |_|\___|\___|_|\_\_| \_|\__,_|_| \_|
  */

  #define LINE_LINE_LINE_LINE "--------------------------------------------------------------------------------"

  //! check if the vector `pv` os size `DIM` contains only regular floats. If not an error is issued
  void
  checkNaN(
    double const pv[],
    char   const v_name[],
    int          DIM,
    int          line,
    char   const file[]
  ) {
    for ( int i = 0; i < DIM; ++i ) {
      if ( isInfinite(pv[i]) ) {
        LW_ERROR(
          "{}\n({}):{}) found Infinity at {}[{}]\n{}\n",
          LINE_LINE_LINE_LINE,
          basename(file), line, v_name, i,
          LINE_LINE_LINE_LINE
        );
      } else if ( isNaN(pv[i]) ) {
        LW_ERROR(
          "{}\n({}):{}) found NaN at {}[{}]\n{}\n",
          LINE_LINE_LINE_LINE,
          basename(file), line, v_name, i,
          LINE_LINE_LINE_LINE
        );
      }
    }
  }

  void
  checkNaN(
    float const pv[],
    char  const v_name[],
    int         DIM,
    int         line,
    char const  file[]
  ) {
    for ( int i = 0; i < DIM; ++i ) {
      if ( isInfinite(pv[i]) ) {
        LW_ERROR(
          "{}\n({}):{}) found Infinity at {}[{}]\n{}\n",
          LINE_LINE_LINE_LINE,
          basename(file), line, v_name, i,
          LINE_LINE_LINE_LINE
        );
      } else if ( isNaN(pv[i]) ) {
        LW_ERROR(
          "{}\n({}):{}) found NaN at {}[{}]\n{}\n",
          LINE_LINE_LINE_LINE,
          basename(file), line, v_name, i,
          LINE_LINE_LINE_LINE
        );
      }
    }
  }
} // end namespace lapack_wrapper

///
/// eof: lapack_wrapper_utils.cc
///

