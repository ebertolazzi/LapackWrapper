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
/// file: lapack_wrapper_utils.hh
///

#ifndef LAPACK_WRAPPER_UTILS_HH
#define LAPACK_WRAPPER_UTILS_HH

#include <mutex> // std::mutex

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#endif

#ifndef LW_ERROR0
  #define LW_ERROR0(MSG) \
  throw lapack_wrapper::Runtime_Error( MSG, __FILE__, __LINE__ )
#endif

#ifndef LW_ASSERT0
  #define LW_ASSERT0(COND,MSG) if ( !(COND) ) LW_ERROR0( MSG )
#endif

#ifndef LW_ERROR
  #define LW_ERROR(...) \
  throw lapack_wrapper::Runtime_Error( fmt::format(__VA_ARGS__), __FILE__, __LINE__ )
#endif

#ifndef LW_ASSERT
  #define LW_ASSERT(COND,...) if ( !(COND) ) LW_ERROR( __VA_ARGS__ )
#endif



#ifndef LW_ERROR_TRACE0
  #define LW_ERROR_TRACE0(MSG) \
  throw lapack_wrapper::Runtime_TraceError( MSG, __FILE__, __LINE__ )
#endif

#ifndef LW_ASSERT_TRACE0
  #define LW_ASSERT_TRACE0(COND,MSG) if ( !(COND) ) LW_ERROR_TRACE0( MSG )
#endif

#ifndef LW_ERROR_TRACE
  #define LW_ERROR_TRACE(...) \
  throw lapack_wrapper::Runtime_TraceError( fmt::format(__VA_ARGS__), __FILE__, __LINE__ )
#endif

#ifndef LW_ASSERT_TRACE
  #define LW_ASSERT_TRACE(COND,...) if ( !(COND) ) LW_ERROR_TRACE( __VA_ARGS__ )
#endif

#ifdef LAPACK_WRAPPER_NO_DEBUG
  #ifndef LW_ASSERT0_DEBUG
    #define LW_ASSERT0_DEBUG(COND,MSG)
  #endif
  #ifndef LW_ASSERT_DEBUG
    #define LW_ASSERT_DEBUG(COND,...)
  #endif
#else
  #ifndef LW_ASSERT0_DEBUG
    #define LW_ASSERT0_DEBUG(COND,MSG) LW_ASSERT0(COND,MSG)
  #endif
  #ifndef LW_ASSERT_DEBUG
    #define LW_ASSERT_DEBUG(COND,...) LW_ASSERT(COND,__VA_ARGS__)
  #endif
#endif

namespace lapack_wrapper {

  typedef std::basic_ostream<char> ostream_type;

  class Runtime_TraceError : public std::runtime_error {
  private:
    std::string
    grab_backtrace(
      std::string const & reason,
      char        const   file[],
      int                 line
    ) const;

  public:
    explicit
    Runtime_TraceError( std::string const & reason, char const file[], int line )
    : std::runtime_error( grab_backtrace( reason, file, line ) )
    { }

    virtual const char* what() const noexcept override;
  };

  class Runtime_Error : public std::runtime_error {
  public:
    explicit
    Runtime_Error( std::string const & reason, char const file[], int line )
    : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    { }

    explicit
    Runtime_Error( char const reason[], char const file[], int line )
    : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    { }

    virtual const char* what() const noexcept override;
  };

  class Console {

    mutable std::mutex message_mutex; // mutex for critical section

  public:
    class Console_style {
    public:
      rang::style s;
      rang::fg    f;
      rang::bg    b;
    };

  private:

    ostream_type * p_stream;

    // 0 only fatal, error
    // 1 + warning
    // 2
    int level;
    int n_thread;
    int max_n_thread;

    Console_style message_style;
    Console_style warning_style;
    Console_style error_style;
    Console_style fatal_style;

    Console() = delete;
    Console( Console const & ) = delete;

  public:

    Console( ostream_type * p_stream = &std::cout, int level = 4 );

    void changeLevel( int new_level );
    void changeStream( ostream_type * new_p_stream );
    void changeNthread( int new_n_thread );

    int getLevel() const { return level; }
    int getNthread() const { return n_thread; }
    int getMaxNthread() const { return max_n_thread; }

    ostream_type * getStream() const { return p_stream; }

    Console const &
    message( std::string const & msg, int msg_level = 4 ) const;

    Console const &
    semaphore( unsigned ryg, std::string const & msg, int msg_level = 0 ) const;

    Console const &
    warning( std::string const & msg ) const; // level >= 2

    Console const &
    error( std::string const & msg ) const; // level >= 1

    Console const &
    fatal( std::string const & msg ) const; // level >= 0

    Console const & black   ( std::string const & msg, int msg_level = 0 ) const;
    Console const & red     ( std::string const & msg, int msg_level = 0 ) const;
    Console const & green   ( std::string const & msg, int msg_level = 0 ) const;
    Console const & yellow  ( std::string const & msg, int msg_level = 0 ) const;
    Console const & blue    ( std::string const & msg, int msg_level = 0 ) const;
    Console const & magenta ( std::string const & msg, int msg_level = 0 ) const;
    Console const & cyan    ( std::string const & msg, int msg_level = 0 ) const;
    Console const & gray    ( std::string const & msg, int msg_level = 0 ) const;

    void
    setMessageStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      message_style.s = s;
      message_style.f = f;
      message_style.b = b;
    }

    void
    setWarningStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      warning_style.s = s;
      warning_style.f = f;
      warning_style.b = b;
    }

    void
    setErrorStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      error_style.s = s;
      error_style.f = f;
      error_style.b = b;
    }

    void
    setFatalStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      fatal_style.s = s;
      fatal_style.f = f;
      fatal_style.b = b;
    }

  };

  /*
  //    ____                _              _
  //   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
  //  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
  //  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
  //   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
  */

  /// Not a number constant
  template <typename T> T NaN();
  template <> inline float NaN()
  { return std::numeric_limits<float>::quiet_NaN(); }
  template <> inline double NaN()
  { return std::numeric_limits<double>::quiet_NaN(); }

  /// machine epsilon
  template <typename T> T machineEps();

  template <> inline float
  machineEps() { return std::numeric_limits<float>::epsilon(); }

  template <> inline double
  machineEps() { return std::numeric_limits<double>::epsilon(); }

  /// square root of machine epsilon
  template <typename T> T sqrtMachineEps();

  template <> inline float
  sqrtMachineEps() { return std::sqrt(std::numeric_limits<float>::epsilon()); }

  template <> inline double
  sqrtMachineEps() { return std::sqrt(std::numeric_limits<double>::epsilon()); }

  /// maximum representable value
  template <typename T> T maximumValue();

  template <> inline float
  maximumValue() { return std::sqrt(std::numeric_limits<float>::max()); }

  template <> inline double
  maximumValue() { return std::sqrt(std::numeric_limits<double>::max()); }

  /// minimum representable value
  template <typename T> T minimumValue();

  template <> inline float
  minimumValue() { return std::sqrt(std::numeric_limits<float>::min()); }

  template <> inline double
  minimumValue() { return std::sqrt(std::numeric_limits<double>::min()); }

  static
  inline
  bool isZero( double x ) { return FP_ZERO == std::fpclassify(x); }

  static
  inline
  bool isZero( float x ) { return FP_ZERO == std::fpclassify(x); }

  static
  inline
  bool isInfinite( double x ) { return FP_INFINITE == std::fpclassify(x); }

  static
  inline
  bool isInfinite( float x ) { return FP_INFINITE == std::fpclassify(x); }

  static
  inline
  bool isNaN( double x ) { return FP_NAN == std::fpclassify(x); }

  static
  inline
  bool isNaN( float x ) { return FP_NAN == std::fpclassify(x); }

  static
  inline
  bool isRegular( double x )
  { return !( FP_INFINITE == std::fpclassify(x) ||
              FP_NAN      == std::fpclassify(x) ); }

  static
  inline
  bool isRegular( float x )
  { return !( FP_INFINITE == std::fpclassify(x) ||
              FP_NAN      == std::fpclassify(x) ); }

  static
  inline
  bool isInteger( double x )
  { return isZero( x-static_cast<long>(std::floor(x)) ); }

  static
  inline
  bool isInteger( float x )
  { return isZero( x-static_cast<long>(std::floor(x)) ); }

  static
  inline
  bool isUnsigned( double x ) { return isInteger(x) && x >= 0; }

  static
  inline
  bool isUnsigned( float x ) { return isInteger(x) && x >= 0; }

  //============================================================================

  bool
  foundNaN( double const pv[], int DIM );

  bool
  foundNaN( float const pv[], int DIM );

  void
  checkNaN(
    double const pv[],
    char   const v_name[],
    int          DIM,
    int          line,
    char   const file[]
  );

  void
  checkNaN(
    float const pv[],
    char  const v_name[],
    int         DIM,
    int         line,
    char  const file[]
  );

  //============================================================================

  //! `m_e` the value of \f$ e \f$.
  static double const m_e = 2.718281828459045235360287471352662497757;

  //! `m_pi` the value of \f$ \pi \f$.
  static double const m_pi = 3.141592653589793238462643383279502884197;

  //! `m_2pi` the value of \f$ 2\pi \f$.
  static double const m_2pi = 6.283185307179586476925286766559005768394;

  //! `m_pi_2` the value of \f$ \pi/2 \f$.
  static double const m_pi_2 = 1.570796326794896619231321691639751442098;

  //! `m_pi_4` the value of \f$ \pi/4 \f$.
  static double const m_pi_4 = 0.7853981633974483096156608458198757210492;

  //! `m_1_pi` the value of \f$ 1/\pi \f$.
  static double const m_1_pi = 0.3183098861837906715377675267450287240689;

  //! `m_2_pi` the value of \f$ 2/\pi \f$.
  static double const m_2_pi = 0.6366197723675813430755350534900574481378;

  //! `m_sqrtpi` the value of \f$ \sqrt{\pi} \f$.
  static double const m_sqrtpi = 1.772453850905516027298167483341145182798;

  //! `m_2_sqrtpi` the value of \f$ 2/\sqrt{\pi} \f$.
  static double const m_2_sqrtpi = 1.128379167095512573896158903121545171688;

  //! `m_sqrt2` the value of \f$ \sqrt{2} \f$.
  static double const m_sqrt2 = 1.414213562373095048801688724209698078570;

  //! `m_1_sqrt2` the value of \f$ 1/\sqrt{2} \f$.
  static double const m_1_sqrt2 = 0.7071067811865475244008443621048490392850;

}

#endif

///
/// eof: lapack_wrapper_utils.hh
///
