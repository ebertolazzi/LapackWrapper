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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef TIC_TOC_HH
#define TIC_TOC_HH

#include "../lapack_wrapper_config.hh"

#ifdef LAPACK_WRAPPER_OS_WINDOWS
  #include <windows.h>
  class TicToc {

    typedef double real_type;
    LARGE_INTEGER m_frequency;   // ticks per second
    LARGE_INTEGER m_t1, m_t2;    // ticks
    real_type     m_elapsed_time;

    TicToc( TicToc const & ) = delete;
    TicToc const & operator = ( TicToc const & ) const = delete;

  public:

    TicToc()
    : m_elapsed_time(0)
    { QueryPerformanceFrequency(&m_frequency); tic(); }

    ~TicToc() {}

    void
    tic()
    { QueryPerformanceCounter(&m_t1); }

    void
    toc() {
      QueryPerformanceCounter(&t2);
      m_elapsed_time = (m_t2.QuadPart - m_t1.QuadPart) * 1000.0 / m_frequency.QuadPart;;
    }

    real_type
    elapsed_s() const
    { return 1e-3*m_elapsed_time; }

    real_type
    elapsed_ms() const
    { return m_elapsed_time; }

  };

  inline
  void
  sleep_for_seconds( unsigned s )
  { Sleep(DWORD(s)*1000); }

  inline
  void
  sleep_for_milliseconds( unsigned ms )
  { Sleep(DWORD(ms)); }

#else

  #include <chrono>
  #include <thread>

  class TicToc {

    typedef double real_type;
    #ifdef TIC_TOC_USE_HRC
    typedef std::chrono::high_resolution_clock clock;
    #else
    typedef std::chrono::steady_clock clock;
    #endif

    using elapsed_resolution = std::chrono::microseconds;

    clock::time_point m_start_time;
    clock::time_point m_stop_time;

    elapsed_resolution m_elapsed_time;

    TicToc( TicToc const & );
    TicToc const & operator = ( TicToc const & ) const;

   public:

    TicToc()
    : m_elapsed_time(0)
    { this->tic(); }

    ~TicToc() {}

    void
    tic()
    { m_start_time = clock::now(); }

    void
    toc() {
      m_stop_time    = clock::now();
      m_elapsed_time = std::chrono::duration_cast<elapsed_resolution>(m_stop_time - m_start_time);
    }

    real_type
    elapsed_s() const
    { return 1e-6*m_elapsed_time.count(); }

    real_type
    elapsed_ms() const
    { return 1e-3*m_elapsed_time.count(); }
  };

  inline
  void
  sleep_for_seconds( unsigned s )
  { std::this_thread::sleep_for(std::chrono::seconds(s)); }

  inline
  void
  sleep_for_milliseconds( unsigned ms )
  { std::this_thread::sleep_for(std::chrono::milliseconds(ms)); }

#endif

#endif
