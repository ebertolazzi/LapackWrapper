/*

  \file     sparse_tool.hh
  \mainpage Sparse_tool: a Sparse Matrix Manager
  \date     2011, July 21
  \version  1.1
  \note     based on SparseLib 11.0 (C 1999 -- 2008)

  \author Enrico Bertolazzi

  Affiliations:
  Dipartimento di Ingegneria Industriale
  Universita` degli Studi di Trento
  email: enrico.bertolazzi@unitn.it
*/

#pragma once

#ifndef SPARSETOOL_dot_HH
#define SPARSETOOL_dot_HH

#include <cmath>
#include <complex>
#include <string>
#include <string_view>
#include <cstdint>

// standard includes I/O
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utility>

// STL lib
#include <vector>
#include <algorithm>
#include <functional>


#ifdef NO_SYSTEM_UTILS
  #include "Utils_eigen.hh"
#else
  #include <Utils_eigen.hh>
#endif

#include "lapack_wrapper/lapack_wrapper.hh"

// workaround for windows macros
#ifdef max
  #undef max
#endif

#ifdef min
  #undef min
#endif

/*
// #    #    ##     ####   #####    ####    ####
// ##  ##   #  #   #    #  #    #  #    #  #
// # ## #  #    #  #       #    #  #    #   ####
// #    #  ######  #       #####   #    #       #
// #    #  #    #  #    #  #   #   #    #  #    #
// #    #  #    #   ####   #    #   ####    ####
*/

//!
//! Default number of pre-allocated nonzeros.
//!
#ifndef SPARSETOOL_DEFAULT_NNZ
  #define SPARSETOOL_DEFAULT_NNZ 100
#endif

//!
//! Issue an error message.
//!
#define SPARSETOOL_ERR(W)                             \
  { using namespace ::std;                            \
    cerr << "\n" << W  << "\nin file `" << __FILE__   \
         << "' on line " << __LINE__ << "\n";;        \
    exit(0); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define SPARSETOOL_ASSERT(X,W) if ( !(X) ) SPARSETOOL_ERR(W)

#ifdef SPARSETOOL_DEBUG
  #define SPARSETOOL_INLINE
#else
  #define SPARSETOOL_INLINE inline
#endif

// loops
#define SPARSELIB_LOOP(N,DO)   { for ( integer i = 0; i < N; ++i ) { DO; } }
#define SPARSELIB_V1LOOP(DO)   SPARSELIB_LOOP(this->size(),DO)
#define SPARSELIB_V2LOOP(V,DO) SPARSELIB_LOOP(minIndex(V.size(),this->size()),DO)

#endif

//!
//! The namespace with the Sparse_tool toolkit.
//!
namespace Sparse_tool {

  using Utils::ostream_type;
  using Utils::istream_type;

  using std::vector;
  using std::string;
  using std::string_view;
  using std::istream;
  using std::ostream;
  using std::cin;
  using std::cout;
  using std::cerr;
  using std::setw;
  using std::greater;
  using std::less;

  // FROM BLITZ++

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  //!
  //! This class is copied from `BLITZ++`.
  //! The C++ compiler supports partial specialization, so type promotion
  //! can be done the elegant way.
  //! This implementation is after ideas by Jean-Louis Leroy.
  //!

  template<typename T>
  struct precision_trait {
    static unsigned const precisionRank     = 0;     //!< precision rank
    static bool     const knowPrecisionRank = false; //!< is rank known ?
  };

  #define SPARSELIB_DECLARE_PRECISION(T,rank)         \
    template<>                                        \
    struct precision_trait< T > {                     \
      static unsigned const precisionRank     = rank; \
      static bool     const knowPrecisionRank = true; \
    };

  SPARSELIB_DECLARE_PRECISION(int,100)
  SPARSELIB_DECLARE_PRECISION(unsigned int,200)
  SPARSELIB_DECLARE_PRECISION(long,300)
  SPARSELIB_DECLARE_PRECISION(unsigned long,400)
  SPARSELIB_DECLARE_PRECISION(float,500)
  SPARSELIB_DECLARE_PRECISION(double,600)
  SPARSELIB_DECLARE_PRECISION(long double,700)
  SPARSELIB_DECLARE_PRECISION(std::complex<float> ,800)
  SPARSELIB_DECLARE_PRECISION(std::complex<double> ,900)
  SPARSELIB_DECLARE_PRECISION(std::complex<long double> ,1000)

  #endif

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  template<typename T>
  struct return_trait {
    using real_type = T; //!< type of the trait
  };

  // Partial specialization for `complex` type.
  template<typename T> struct return_trait<std::complex<T> > { using real_type = T; };

  //! \internal Class for type promotion
  template<typename T>
  struct autopromote_trait { using T_numtype = T; };

  //! promote `bool` to `int`
  template<> struct autopromote_trait<bool>
  { using T_numtype = int; };

  //! promote `char` to `int`
  template<> struct autopromote_trait<char>
  { using T_numtype = int; };

  //! promote `bool` to `int`
  template<> struct autopromote_trait<unsigned char>
  { using T_numtype = int; };

  //! promote `short` to `int`
  template<> struct autopromote_trait<short>
  { using T_numtype = int; };

  // promote `unsigned` `short` to `unsigned`
  template<> struct autopromote_trait<unsigned short>
  { using T_numtype = unsigned; };

  #endif

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <typename T1, typename T2>
  class promote_trait {
  public:
    //! promoted type
    typedef decltype( std::declval<T1>() * std::declval<T2>() ) T_promote;
  };

  //! Class for promotion with three types
  template<typename T1, typename T2, typename T3>
  class promote_trait3 {
  public:
    //! promoted type
    typedef decltype( std::declval<T1>() * std::declval<T2>() * std::declval<T3>() ) T_promote;
  };

  #ifndef SPARSELIB_FUNCTION_TYPE
  #define SPARSELIB_FUNCTION_TYPE float
  #endif

  #ifndef SPARSELIB_ACCUMULATOR_TYPE
  #define SPARSELIB_ACCUMULATOR_TYPE double
  #endif

  template <typename T> class Vector;

  #endif

  //!
  //! Use `lapack_wrapper::integer` as the type for indexing vector and matrices.
  //!
  using lapack_wrapper::integer;

  //!
  //! Return minimum value between `A` and `b`.
  //!
  inline integer
  minIndex(integer a, integer b)
  { return a < b ? a : b; }

  //!
  //! Return maximum value between `A` and `b`.
  //!
  inline integer
  maxIndex(integer a, integer b)
  { return a > b ? a : b; }

  //!
  //! \defgroup Comparator Comparator for Sparse Matrix
  //! costructor and conversion
  //!
  //@{

  //!
  //! Comparator class for selecting **all** elements
  //!
  struct all_ok {
    //! return always true
    bool operator () ( integer /* i */, integer /* j */ ) const { return true; }
  };

  //!
  //! Comparator class for selecting;lements **under** the diagonal
  //!
  struct lower_ok {
    //! select lower diagonal indexes
    bool operator () ( integer i, integer j ) const { return i > j; }
  };

  //!
  //! Comparator class for selecting elements **under** the diagonal.
  //!
  struct lowereq_ok {
    //! select lower diagonal indexes
    bool operator () ( integer i, integer j ) const { return i >= j; }
  };

  //!
  //! Comparator class for selec;ng elements **over** the diagonal.
  //!
  struct upper_ok {
    //! select upper diagonal indexes
    bool operator () ( integer i, integer j ) const { return i < j; }
  };

  //!
  //! Comparator class for selecting elements **over** the diagonal.
  //!
  struct uppereq_ok {
    //! select upper diagonal indexes
    bool operator () ( integer i, integer j ) const { return i <= j; }
  };

  //!
  //! Comparator class for selecting elements **on** the diagonal.
  //!
  struct diag_ok {
    //! select diagonal indexes
    bool operator () ( integer i, integer j ) const { return i==j; }
  };

  //@}

  /*
  //  ####   #####     ##    #####    ####   ######
  // #       #    #   #  #   #    #  #       #
  //  ####   #    #  #    #  #    #   ####   #####
  //      #  #####   ######  #####        #  #
  // #    #  #       #    #  #   #   #    #  #
  //  ####   #       #    #  #    #   ####   ######
  //
  //  #####  #####   ######  ######
  //    #    #    #  #       #
  //    #    #    #  #####   #####
  //    #    #####   #       #
  //    #    #   #   #       #
  //    #    #    #  ######  ######
  */

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  //!
  //! \defgroup MVStructures Structures Storing Matrix/Vector operations.
  //!
  //@{

  //!
  //! Structure storing the operation `a/M`
  //! (vector-matrix division i.e. solve `A*x=a`)
  //!
  template <typename VA, typename MATRIX>
  struct Vector_V_div_M {
    VA     const & a;
    MATRIX const & M;
    Vector_V_div_M( VA const & m_a, MATRIX const & m_M )
    : a(m_a)
    , M(m_M)
    {}
  };

  //!
  //! Structure storing the operation `M*a` (matrix-vector multiply)
  //!
  template <typename MATRIX, typename VA>
  struct Vector_M_mul_V {
    MATRIX const & M;
    VA     const & a;
    Vector_M_mul_V( MATRIX const & m_M, VA const & m_a )
    : M(m_M)
    , a(m_a)
    {}
  };

  //!
  //! Structure storing the operation `s*M` (scalar-matrix multiply)
  //!
  template <typename S, typename MATRIX>
  struct Vector_S_mul_M {
    S      const & s;
    MATRIX const & M;
    Vector_S_mul_M( S const & m_s, MATRIX const & m_M )
    : s(m_s)
    , M(m_M)
    {}
  };

  //!
  //! Structure storing the operation `M^a`
  //! (transpose matrix-vector multiply).
  //!
  template <typename MATRIX, typename VA>
  struct Vector_Mt_mul_V {
    MATRIX const & M;
    VA     const & a;
    Vector_Mt_mul_V( MATRIX const & m_M, VA const & m_a )
    : M(m_M)
    , a(m_a)
    {}
  };

  //!
  //! Structure storing the operation `s*(M*a)`
  //! (scalar-matrix-vector multiply).
  //!
  template <typename SCALAR, typename MATRIX, typename VA>
  struct Vector_S_mul_M_mul_V {
    SCALAR const & s;
    MATRIX const & M;
    VA     const & a;
    Vector_S_mul_M_mul_V(
      SCALAR const & m_s,
      MATRIX const & m_M,
      VA     const & m_a
    )
    : s(m_s)
    , M(m_M)
    , a(m_a)
    {}
  };

  //!
  //! Structure storing the operation `s*(M^a)`
  //! (transpose scalar-matrix-vector multiply).
  //!
  template <typename SCALAR, typename MATRIX, typename VA>
  struct Vector_S_mul_Mt_mul_V {
    SCALAR const & s;
    MATRIX const & M;
    VA     const & a;
    Vector_S_mul_Mt_mul_V(
      SCALAR const & m_s,
      MATRIX const & m_M,
      VA     const & m_a
    )
    : s(m_s)
    , M(m_M)
    , a(m_a)
    {}
  };

  //!
  //! Structure storing the operation `a+M*b`
  //!
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sum_M_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_M_mul_V(
      VA     const & m_a,
      MATRIX const & m_M,
      VB     const & m_b
    )
    : a(m_a)
    , M(m_M)
    , b(m_b)
    {}
  };

  //! structure storing the operation `a+M^b`
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sum_Mt_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_Mt_mul_V(
      VA     const & m_a,
      MATRIX const & m_M,
      VB     const & m_b
    )
    : a(m_a)
    , M(m_M)
    , b(m_b)
    {}
  };

  //!
  //! Structure storing the operation `a-M*b`.
  //!
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sub_M_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sub_M_mul_V(
      VA     const & m_a,
      MATRIX const & m_M,
      VB     const & m_b
    )
    : a(m_a)
    , M(m_M)
    , b(m_b)
    {}
  };

  //!
  //! Structure storing the operation `a+M^b`.
  //!
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sub_Mt_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sub_Mt_mul_V(
      VA     const & m_a,
      MATRIX const & m_M,
      VB     const & m_b
    )
    : a(m_a)
    , M(m_M)
    , b(m_b)
    {}
  };

  //!
  //! Structure storing the operation `a+s*M*b`.
  //!
  template <typename VA, typename SCALAR, typename MATRIX, typename VB>
  struct Vector_V_sum_S_mul_M_mul_V {
    VA     const & a;
    SCALAR const & s;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_S_mul_M_mul_V(
      VA     const & m_a,
      SCALAR const & m_s,
      MATRIX const & m_M,
      VB     const & m_b
    )
    : a(m_a)
    , s(m_s)
    , M(m_M)
    , b(m_b)
    {}
  };

  //!
  //! Structure storing the operation `a+s*M^b`.
  //!
  template <typename VA, typename SCALAR, typename MATRIX, typename VB>
  struct Vector_V_sum_S_mul_Mt_mul_V {
    VA     const & a;
    SCALAR const & s;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_S_mul_Mt_mul_V(
      VA     const & m_a,
      SCALAR const & m_s,
      MATRIX const & m_M,
      VB     const & m_b
    )
    : a(m_a)
    , s(m_s)
    , M(m_M)
    , b(m_b)
    {}
  };

  //!
  //! Structure storing the operation `a/P`,
  //! i.e. the application of Preconditioner.
  //!
  template <typename VECTOR, typename PRECO>
  struct Vector_V_div_P {
    VECTOR const & a;
    PRECO  const & P;
    Vector_V_div_P( VECTOR const & m_a, PRECO const & m_P )
    : a(m_a)
    , P(m_P)
    {}
  };

  //@}

  #endif

  template <typename T>
  class Vector : public Eigen::Matrix<T,Eigen::Dynamic,1> {
  public:

    using V_base    = Eigen::Matrix<T,Eigen::Dynamic,1>;
    using real_type = T;

    Vector() : Eigen::Matrix<T,Eigen::Dynamic,1>(0) {}
    explicit Vector( integer dim ) : Eigen::Matrix<T,Eigen::Dynamic,1>(dim) {}

    Vector( Vector<T> const & v ) {
      this->to_eigen() = v.to_eigen();
    }

    V_base       & to_eigen()       { return *static_cast<V_base*>(this); }
    V_base const & to_eigen() const { return *static_cast<V_base const *>(this); }

    Vector<T> const &
    operator = ( T const & v )
    { this->to_eigen().array() = v; return *this; }

    Vector<T> const &
    operator += ( T const & v )
    { this->to_eigen().array() += v; return *this; }

    Vector<T> const &
    operator -= ( T const & v )
    { this->to_eigen().array() -= v; return *this; }

    Vector<T> const &
    operator *= ( T const & v )
    { this->to_eigen().array() *= v; return *this; }

    Vector<T> const &
    operator /= ( T const & v )
    { this->to_eigen().array() /= v; return *this; }

    Vector<T> const &
    operator = ( Vector<T> const & v )
    { this->to_eigen() = v.to_eigen(); return *this; }

    Vector<T> const &
    operator += ( Vector<T> const & v )
    { this->to_eigen() += v.to_eigen(); return *this; }

    Vector<T> const &
    operator -= ( Vector<T> const & v )
    { this->to_eigen() -= v.to_eigen(); return *this; }

    Vector<T> const &
    operator *= ( Vector<T> const & v )
    { this->to_eigen().array() *= v.to_eigen().array(); return *this; }

    Vector<T> const &
    operator /= ( Vector<T> const & v )
    { this->to_eigen().array() /= v.to_eigen().array(); return *this; }

    Vector<T> const &
    operator = ( V_base const & v )
    { this->to_eigen() = v; return *this; }

    Vector<T> const &
    operator += ( V_base const & v )
    { this->to_eigen() += v; return *this; }

    Vector<T> const &
    operator -= ( V_base const & v )
    { this->to_eigen() -= v; return *this; }

    Vector<T> const &
    operator *= ( V_base const & v )
    { this->to_eigen().array() *= v.array(); return *this; }

    Vector<T> const &
    operator /= ( V_base const & v )
    { this->to_eigen().array() /= v.array(); return *this; }

    //! \name MM_COMPLEX MATRIX OPERATIONS
    ///@{

    //!
    //! Assign to `*this`  the evaluation of the expression `a/M`
    //! contained in `op` of type `Vector_V_div_M`
    //!
    template <typename MATRIX, typename VA>
    Vector<T> const &
    operator = ( Vector_V_div_M<VA,MATRIX> const & op )
    { op.M.ass_V_div_M( *this, op.a ); return *this; }

    //!
    //! Assign to `*this`  the evaluation of the expression `M*a`
    //! contained in `op` of type `Vector_M_mul_V`
    //!
    template <typename MATRIX, typename VA>
    Vector<T> const &
    operator = ( Vector_M_mul_V<MATRIX,VA> const & op ) {
      this->setZero();
      op.M.add_S_mul_M_mul_V( *this, real_type(+1), op.a );
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `a/M`
    //! contained in `op` of type `Vector_V_div_M`
    //!
    template <typename MATRIX, typename VA>
    Vector<T> const &
    operator += ( Vector_M_mul_V<MATRIX,VA> const & op ) {
      op.M.add_S_mul_M_mul_V( *this, real_type(+1), op.a );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `M*a`
    //! contained in `op` of type `Vector_M_mul_V`
    //!
    template <typename MATRIX, typename VA>
    Vector<T> const &
    operator -= ( Vector_M_mul_V<MATRIX,VA> const & op ) {
      op.M.add_S_mul_M_mul_V( *this, real_type(-1), op.a );
      return *this;
    }

    //!
    //! Assign to `*this` the evaluation of the expression `M^T*a`
    //! contained in `op` of type `Vector_Mt_mul_V`
    //!
    template <typename MATRIX, typename VA>
    Vector<T> const &
    operator = ( Vector_Mt_mul_V<MATRIX,VA> const & op ) {
      this->setZero();
      op.M.add_S_mul_Mt_mul_V( *this, real_type(+1), op.a );
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `M^T*a`
    //! contained in `op` of type `Vector_Mt_mul_V`
    //!
    template <typename MATRIX, typename VA>
    Vector<T> const &
    operator += ( Vector_Mt_mul_V<MATRIX,VA> const & op ) {
      op.M.add_S_mul_Mt_mul_V( *this, real_type(+1), op.a );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `M^T*a`
    //! contained in `op` of type `Vector_Mt_mul_V`
    //!
    template <typename MATRIX, typename VA>
    Vector<T> const &
    operator -= (Vector_Mt_mul_V<MATRIX,VA> const & op) {
      op.M.add_S_mul_Mt_mul_V(*this, real_type(-1), op.a);
      return *this;
    }

    //!
    //! Assign to `*this` the evaluation of the expression `s*(M*a)`
    //! contained in `op` of type `Vector_S_mul_M_mul_V`
    //!
    template <typename SCALAR, typename MATRIX, typename VA>
    Vector<T> const &
    operator = ( Vector_S_mul_M_mul_V<SCALAR,MATRIX,VA> const & op ) {
      this->setZero();
      op.M.add_S_mul_M_mul_V( *this, op.s, op.a );
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `s*(M*a)`
    //! contained in `op` of type `Vector_S_mul_M_mul_V`
    //!
    template <typename SCALAR, typename MATRIX, typename VA>
    Vector<T> const &
    operator += ( Vector_S_mul_M_mul_V<SCALAR,MATRIX,VA> const & op ) {
      op.M.add_S_mul_M_mul_V( *this, op.s, op.a );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `s*(M*a)`
    //! contained in `op` of type `Vector_S_mul_M_mul_V`
    //!
    template <typename SCALAR, typename MATRIX, typename VA>
    Vector<T> const &
    operator -= ( Vector_S_mul_M_mul_V<SCALAR,MATRIX,VA> const & op ) {
      op.M.add_S_mul_M_mul_V( *this, -op.s, op.a );
      return *this;
    }

    //!
    //! Assign to `*this` the evaluation of the expression `s*(M^T*a)`
    //! contained in `op` of type `Vector_S_mul_Mt_mul_V`
    //!
    template <typename SCALAR, typename MATRIX, typename VA>
    Vector<T> const &
    operator = ( Vector_S_mul_Mt_mul_V<SCALAR,MATRIX,VA> const & op ) {
      this->setZero();
      op.M.add_S_mul_Mt_mul_V( *this, op.s, op.a );
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `s*(M^T*a)`
    //! contained in `op` of type `Vector_S_mul_Mt_mul_V`
    //!
    template <typename SCALAR, typename MATRIX, typename VA>
    Vector<T> const &
    operator += ( Vector_S_mul_Mt_mul_V<SCALAR,MATRIX,VA> const & op ) {
      op.M.add_S_mul_Mt_mul_V( *this, op.s, op.a );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `s*(M^T*a)`
    //! contained in `op` of type `Vector_S_mul_Mt_mul_V`
    //!
    template <typename SCALAR, typename MATRIX, typename VA>
    Vector<T> const &
    operator -= ( Vector_S_mul_Mt_mul_V<SCALAR,MATRIX,VA> const & op ) {
      op.M.add_S_mul_Mt_mul_V( *this, -op.s, op.a );
      return *this;
    }

    //!
    //! Assign to `*this` the evaluation of the expression `a+M*b`
    //! contained in `op` of type `Vector_V_sum_M_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator = ( Vector_V_sum_M_mul_V<VA,MATRIX,VB> const & op ) {
      *this = op.a;
      op.M.add_S_mul_M_mul_V( *this, real_type(+1), op.b );
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `a+M*b`
    //! contained in `op` of type `Vector_V_sum_M_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator += ( Vector_V_sum_M_mul_V<VA,MATRIX,VB> const & op ) {
      *this += op.a;
      op.M.add_S_mul_M_mul_V( *this, real_type(+1), op.b );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `a+M*b`
    //! contained in `op` of type `Vector_V_sum_M_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator -= ( Vector_V_sum_M_mul_V<VA,MATRIX,VB> const & op ) {
      *this -= op.a;
      op.M.add_S_mul_M_mul_V( *this, real_type(-1), op.b );
      return *this;
    }

    //!
    //! Assign to `*this` the evaluation of the expression `a+M^T*b`
    //! contained in `op` of type `Vector_V_sum_Mt_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator = ( Vector_V_sum_Mt_mul_V<VA,MATRIX,VB> const & op ) {
      *this = op.a;
      op.M.add_S_mul_Mt_mul_V( *this, real_type(+1), op.b );
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `a+M^T*b`
    //! contained in `op` of type `Vector_V_sum_Mt_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator += ( Vector_V_sum_Mt_mul_V<VA,MATRIX,VB> const & op ) {
      *this += op.a;
      op.M.add_S_mul_Mt_mul_V( *this, real_type(+1), op.b );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `a+M^T*b`
    //! contained in `op` of type `Vector_V_sum_Mt_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator -= ( Vector_V_sum_Mt_mul_V<VA,MATRIX,VB> const & op ) {
      *this -= op.a;
      op.M.add_S_mul_Mt_mul_V( *this, real_type(-1), op.b );
      return *this;
    }

    //!
    //! Assign to `*this` the evaluation of the expression `a-M*b`
    //! contained in `op` of type `Vector_V_sub_M_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator = ( Vector_V_sub_M_mul_V<VA,MATRIX,VB> const & op ) {
      *this = op.a;
      op.M.add_S_mul_M_mul_V( *this, real_type(-1), op.b );
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `a-M*b`
    //! contained in `op` of type `Vector_V_sub_M_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator += ( Vector_V_sub_M_mul_V<VA,MATRIX,VB> const & op ) {
      *this += op.a;
      op.M.add_S_mul_M_mul_V( *this, real_type(-1), op.b );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `a-M*b`
    //! contained in `op` of type `Vector_V_sub_M_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator -= ( Vector_V_sub_M_mul_V<VA,MATRIX,VB> const & op ) {
      *this -= op.a;
      op.M.add_S_mul_M_mul_V( *this, real_type(+1), op.b );
      return *this;
    }

    //!
    //! Assign to `*this` the evaluation of the expression `a-M^T*b`
    //! contained in `op` of type `Vector_V_sub_Mt_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator = ( Vector_V_sub_Mt_mul_V<VA,MATRIX,VB> const & op ) {
      *this = op.a;
      op.M.add_S_mul_Mt_mul_V( *this, real_type(-1), op.b);
      return *this;
    }

    //!
    //! Add to `*this` the evaluation of the expression `a-M^T*b`
    //! contained in `op` of type `Vector_V_sub_Mt_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator += ( Vector_V_sub_Mt_mul_V<VA,MATRIX,VB> const & op ) {
      *this += op.a;
      op.M.add_S_mul_Mt_mul_V( *this, real_type(-1), op.b );
      return *this;
    }

    //!
    //! Subtract to `*this`  the evaluation of the expression `a-M^T*b`
    //! contained in `op` of type `Vector_V_sub_Mt_mul_V`
    //!
    template <typename VA, typename MATRIX, typename VB>
    Vector<T> const &
    operator -= ( Vector_V_sub_Mt_mul_V<VA,MATRIX,VB> const & op ) {
      *this -= op.a;
      op.M.add_S_mul_Mt_mul_V( *this, real_type(+1), op.b );
      return *this;
    }

    //!
    //! Assign to `*this`  the evaluation of the expression `a+s*(M*b)`
    //! contained in `op` of type `Vector_V_sum_S_mul_M_mul_V`
    //!
    template <typename VA, typename SCALAR, typename MATRIX, typename VB>
    Vector<T> const &
    operator = (Vector_V_sum_S_mul_M_mul_V<VA,SCALAR,MATRIX,VB> const & op) {
      *this = op.a;
      op.M.add_S_mul_M_mul_V( *this, op.s, op.b );
      return *this;
    }

    //!
    //! Add to `*this`  the evaluation of the expression `a+s*(M*b)`
    //! contained in `op` of type `Vector_V_sum_S_mul_M_mul_V`
    //!
    template <typename VA, typename SCALAR, typename MATRIX, typename VB>
    Vector<T> const &
    operator += ( Vector_V_sum_S_mul_M_mul_V<VA,SCALAR,MATRIX,VB> const & op ) {
      *this += op.a;
      op.M.add_S_mul_M_mul_V( *this, op.s, op.b );
      return *this;
    }

    //!
    //! Subtract to `*this` the evaluation of the expression `a+s*(M*b)`
    //! contained in `op` of type `Vector_V_sum_S_mul_M_mul_V`
    //!
    template <typename VA, typename SCALAR, typename MATRIX, typename VB>
    Vector<T> const &
    operator -= ( Vector_V_sum_S_mul_M_mul_V<VA,SCALAR,MATRIX,VB> const & op ) {
      *this -= op.a;
      op.M.add_S_mul_M_mul_V( *this, -op.s, op.b );
      return *this;
    }

    //!
    //! Assign to `*this`  the evaluation of the expression `a+s*(M^T*b)`
    //! contained in `op` of type `Vector_V_sum_S_mul_M_mul_V`
    //!
    template <typename VA, typename SCALAR, typename MATRIX, typename VB>
    Vector<T> const &
    operator = ( Vector_V_sum_S_mul_Mt_mul_V<VA,SCALAR,MATRIX,VB> const & op ) {
      *this = op.a;
      op.M.add_S_mul_Mt_mul_V( *this, op.s, op.b );
      return *this;
    }

    //!
    //! Add to `*this`  the evaluation of the expression `a+s*(M^T*b)`
    //! contained in `op` of type `Vector_V_sum_S_mul_M_mul_V`
    //!
    template <typename VA, typename SCALAR, typename MATRIX, typename VB>
    Vector<T> const &
    operator += ( Vector_V_sum_S_mul_Mt_mul_V<VA,SCALAR,MATRIX,VB> const & op ) {
      *this += op.a;
      op.M.add_S_mul_Mt_mul_V( *this, op.s, op.b );
      return *this;
    }

    //!
    //! Subtract to `*this`  the evaluation of the expression `a+s*(M^T*b)`
    //! contained in `op` of type `Vector_V_sum_S_mul_M_mul_V`
    //!
    template <typename VA, typename SCALAR, typename MATRIX, typename VB>
    Vector<T> const &
    operator -= ( Vector_V_sum_S_mul_Mt_mul_V<VA,SCALAR,MATRIX,VB> const & op ) {
      *this -= op.a;
      op.M.add_S_mul_Mt_mul_V( *this, -op.s, op.b );
      return *this;
    }

    //!
    //! Assign to `*this`  the evaluation of the preconditioner `a/P`
    //! contained in `op` of type `Vector_V_div_P`
    //!
    template <typename VEC, typename PRECO>
    Vector<T> const &
    operator = ( Vector_V_div_P<VEC,PRECO> const & op ) {
      op.P.ass_preco( *this, op.a );
      return *this;
    }

    //!
    //! Assign to `*this`  the evaluation of the preconditioner `*this/P`
    //!
    template <typename PRECO>
    Vector<T> const &
    operator /= ( PRECO const & P ) {
      P.ass_preco( *this );
      return *this;
    }

    ///@}
  };

  /*
  //  #####  #     #  ###   #####  #    #        #####  ####### ######  #######
  // #     # #     #   #   #     # #   #        #     # #     # #     #    #
  // #     # #     #   #   #       #  #         #       #     # #     #    #
  // #     # #     #   #   #       ###           #####  #     # ######     #
  // #   # # #     #   #   #       #  #               # #     # #   #      #
  // #    #  #     #   #   #     # #   #        #     # #     # #    #     #
  //  #### #  #####   ###   #####  #    #        #####  ####### #     #    #
  */

  //!
  //! \name Sorting Structures.
  //!
  //@{

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  using stack_node = struct { integer lo, hi; };

  static inline
  bool GT( integer i, integer i1 )
  { return i>i1; }

  static inline
  bool GT( integer i,  integer j,
           integer i1, integer j1 )
  { return j==j1 ? i>i1 : j>j1; }

  #endif

  //!
  //! Function to sort two vector of indexes and values
  //! ascending respect to index vector `I`.
  //!
  //! \param I           C-array of row index vector
  //! \param A           C-array of values
  //! \param total_elems total number of element to sort
  //! \param MAX_THRESH  threshold value to swich to insert sort
  //!
  template <typename I_type, typename T_type>
  static
  void
  QuickSortI(
    I_type  I[],
    T_type  A[],
    integer total_elems,
    integer MAX_THRESH = 4
  ) {

    using ::std::swap;

    if ( total_elems <= 1 ) return;
    if ( total_elems <= MAX_THRESH ) goto insert_sort;
    {
      integer lo = 0;
      integer hi = total_elems - 1;

      stack_node stack[128];
      stack_node *top = stack;

      while ( top >= stack ) { // Stack not empty

        /* Select median value from among LO, MID, and HI. Rearrange *\
         * LO and HI so the three values are sorted. This lowers the *
         * probability of picking a pathological pivot value and     *
        \* skips a comparison for both the LEFT_PTR and RIGHT_PTR.   */

        I_type & I_hi  = I[hi];
        I_type & I_lo  = I[lo];
        I_type & I_mid = I[ (hi + lo) / 2 ];

        T_type & A_hi  = A[hi];
        T_type & A_lo  = A[lo];
        T_type & A_mid = A[ (hi + lo) / 2 ];

        if ( GT(I_lo, I_mid) ) {
          swap(I_mid, I_lo);
          swap(A_mid, A_lo);
        }
        if ( GT(I_mid, I_hi) ) {
          swap(I_hi, I_mid);
          swap(A_hi, A_mid);
          if ( GT(I_lo,I_mid) ) {
            swap(I_mid, I_lo);
            swap(A_mid, A_lo);
          }
        }

        integer Pivot     = I_mid;
        integer left_ptr  = lo + 1;
        integer right_ptr = hi - 1;

        /* Here's the famous ``collapse the walls'' section of quicksort. *\
         * Gotta like those tight inner loops!  They are the main reason  *
        \* that this algorithm runs much faster than others.              */

        do {

          while ( GT(Pivot, I[left_ptr])  ) ++left_ptr;
          while ( GT(I[right_ptr], Pivot) ) --right_ptr;

          if ( left_ptr < right_ptr ) {
            swap(I[left_ptr], I[right_ptr]);
            swap(A[left_ptr], A[right_ptr]);
            ++left_ptr;
            --right_ptr;
          } else if ( left_ptr == right_ptr ) {
            ++left_ptr;
            --right_ptr;
            break;
          }

        } while ( left_ptr <= right_ptr );

        /* Set up pointers for next iteration.  First determine whether   *\
         * left and right partitions are below the threshold size. If so, *
         * ignore one or both.  Otherwise, push the larger partition's    *
        \* bounds on the stack and continue sorting the smaller one.      */

        if ( (right_ptr - lo) <= MAX_THRESH ) {
          if ((hi - left_ptr) <= MAX_THRESH ) {
            --top;
            if ( top < stack ) break;
            lo = top -> lo;
            hi = top -> hi;
          } else {
            lo = left_ptr;
          }
        } else if ((hi - left_ptr) <= MAX_THRESH) {
          hi = right_ptr;
        } else if ((right_ptr - lo) > (hi - left_ptr)) {
          top -> lo = lo;
          top -> hi = right_ptr;
          ++top;
          lo = left_ptr;
        } else {
          top -> lo = left_ptr;
          top -> hi = hi;
          ++top;
          hi = right_ptr;
        }
      }

      /* Once the BASE_PTR array is partially sorted by quicksort the rest  *\
       * is completely sorted using insertion sort, since this is efficient *
      \* for partitions below MAX_THRESH size.                              */
    }

  insert_sort:

    for ( integer i = 1; i < total_elems; ++i ) {
      for ( integer j = i; j > 0; --j ) {
        if ( GT(I[j], I[j-1]) ) break;
        swap( I[j], I[j-1] );
        swap( A[j], A[j-1] );
      }
    }
  }

  //!
  //! Function to sort two vector of indexes
  //! ascending respect to index vector `I` and `J`.
  //!
  //! \param I           C-array of row index vector
  //! \param J           C-array of column index vector
  //! \param total_elems total number of element to sort
  //! \param MAX_THRESH  threshold value to swich to insert sort
  //!
  template <typename I_type>
  static
  void
  QuickSortIJ(
    I_type  I[],
    I_type  J[],
    integer total_elems,
    integer MAX_THRESH = 4
  ) {

    using ::std::swap;

    if ( total_elems <= 1 ) return;
    if ( total_elems <= MAX_THRESH ) goto insert_sort;

    {
      integer lo = 0;
      integer hi = total_elems - 1;

      stack_node stack[128];
      stack_node *top = stack;

      while ( top >= stack ) { // Stack not empty

        /* Select median value from among LO, MID, and HI. Rearrange *\
         * LO and HI so the three values are sorted. This lowers the *
         * probability of picking a pathological pivot value and     *
        \* skips a comparison for both the LEFT_PTR and RIGHT_PTR.   */

        I_type & I_hi  = I[hi];
        I_type & I_lo  = I[lo];
        I_type & I_mid = I[ (hi + lo) / 2 ];

        I_type & J_hi  = J[hi];
        I_type & J_lo  = J[lo];
        I_type & J_mid = J[ (hi + lo) / 2 ];

        if ( GT(I_lo, J_lo, I_mid, J_mid) ) {
          swap(I_mid, I_lo);
          swap(J_mid, J_lo);
        }
        if ( GT(I_mid, J_mid, I_hi, J_hi) ) {
          swap(I_hi, I_mid);
          swap(J_hi, J_mid);
          if ( GT(I_lo, J_lo, I_mid, J_mid) ) {
            swap(I_mid, I_lo);
            swap(J_mid, J_lo);
          }
        }

        integer IPivot    = I_mid;
        integer JPivot    = J_mid;
        integer left_ptr  = lo + 1;
        integer right_ptr = hi - 1;

        /* Here's the famous ``collapse the walls'' section of quicksort. *\
         * Gotta like those tight inner loops!  They are the main reason  *
        \* that this algorithm runs much faster than others.              */

        do {

          while ( GT(IPivot, JPivot, I[left_ptr], J[left_ptr] )  ) ++left_ptr;
          while ( GT(I[right_ptr], J[right_ptr], IPivot, JPivot) ) --right_ptr;

          if ( left_ptr < right_ptr ) {
            swap(I[left_ptr], I[right_ptr]);
            swap(J[left_ptr], J[right_ptr]);
            ++left_ptr;
            --right_ptr;
          } else if ( left_ptr == right_ptr ) {
            ++left_ptr;
            --right_ptr;
            break;
          }

        } while ( left_ptr <= right_ptr );

        /* Set up pointers for next iteration.  First determine whether   *\
         * left and right partitions are below the threshold size. If so, *
         * ignore one or both.  Otherwise, push the larger partition's    *
        \* bounds on the stack and continue sorting the smaller one.      */

        if ( (right_ptr - lo) <= MAX_THRESH ) {
          if ((hi - left_ptr) <= MAX_THRESH ) {
            --top;
            if ( top < stack ) break;
            lo = top -> lo;
            hi = top -> hi;
          } else {
            lo = left_ptr;
          }
        } else if ((hi - left_ptr) <= MAX_THRESH) {
          hi = right_ptr;
        } else if ((right_ptr - lo) > (hi - left_ptr)) {
          top -> lo = lo;
          top -> hi = right_ptr;
          ++top;
          lo = left_ptr;
        } else {
          top -> lo = left_ptr;
          top -> hi = hi;
          ++top;
          hi = right_ptr;
        }
      }

      /* Once the BASE_PTR array is partially sorted by quicksort the rest   *\
       * is completely sorted using insertion sort, since this is efficient  *
       * for partitions below MAX_THRESH size.                               */
    }

  insert_sort:

    for ( integer i = 1; i < total_elems; ++i ) {
      for ( integer j = i; j > 0; --j ) {
        if ( GT(I[j], J[j], I[j-1], J[j-1]) ) break;
        swap( I[j], I[j-1] );
        swap( J[j], J[j-1] );
      }
    }
  }

  //!
  //! Function to sort three vector of indexes and values
  //! ascending respect to index vector `I` and `J`.
  //!
  //! \param I           C-array of row index vector
  //! \param J           C-array of column index vector
  //! \param A           C-array of values
  //! \param total_elems total number of element to sort
  //! \param MAX_THRESH  threshold value to swich to insert sort
  //!
  template <typename I_type, typename T_type>
  static
  void
  QuickSortIJ2(
    I_type  I[],
    I_type  J[],
    T_type  A[],
    integer total_elems,
    integer MAX_THRESH = 4
  ) {

    using ::std::swap;

    if ( total_elems <= 1 ) return;
    if ( total_elems <= MAX_THRESH ) goto insert_sort;

    {
      integer lo = 0;
      integer hi = total_elems - 1;

      stack_node stack[128];
      stack_node *top = stack;

      while ( top >= stack ) { // Stack not empty

        /* Select median value from among LO, MID, and HI. Rearrange *\
         * LO and HI so the three values are sorted. This lowers the *
         * probability of picking a pathological pivot value and     *
        \* skips a comparison for both the LEFT_PTR and RIGHT_PTR.   */

        I_type & I_hi  = I[hi];
        I_type & I_lo  = I[lo];
        I_type & I_mid = I[ (hi + lo) / 2 ];

        I_type & J_hi  = J[hi];
        I_type & J_lo  = J[lo];
        I_type & J_mid = J[ (hi + lo) / 2 ];

        T_type & A_hi  = A[hi];
        T_type & A_lo  = A[lo];
        T_type & A_mid = A[ (hi + lo) / 2 ];

        if ( GT(I_lo, J_lo, I_mid, J_mid) ) {
          swap(I_mid, I_lo);
          swap(J_mid, J_lo);
          swap(A_mid, A_lo);
        }
        if ( GT(I_mid, J_mid, I_hi, J_hi) ) {
          swap(I_hi, I_mid);
          swap(J_hi, J_mid);
          swap(A_hi, A_mid);
          if ( GT(I_lo, J_lo, I_mid, J_mid) ) {
            swap(I_mid, I_lo);
            swap(J_mid, J_lo);
            swap(A_mid, A_lo);
          }
        }

        integer IPivot    = I_mid;
        integer JPivot    = J_mid;
        integer left_ptr  = lo + 1;
        integer right_ptr = hi - 1;

        /* Here's the famous ``collapse the walls'' section of quicksort. *\
         * Gotta like those tight inner loops!  They are the main reason  *
        \* that this algorithm runs much faster than others.              */

        do {

          while ( GT(IPivot, JPivot, I[left_ptr], J[left_ptr] )  ) ++left_ptr;
          while ( GT(I[right_ptr], J[right_ptr], IPivot, JPivot) ) --right_ptr;

          if ( left_ptr < right_ptr ) {
            swap(I[left_ptr], I[right_ptr]);
            swap(J[left_ptr], J[right_ptr]);
            swap(A[left_ptr], A[right_ptr]);
            ++left_ptr;
            --right_ptr;
          } else if ( left_ptr == right_ptr ) {
            ++left_ptr;
            --right_ptr;
            break;
          }

        } while ( left_ptr <= right_ptr );

        /* Set up pointers for next iteration.  First determine whether   *\
         * left and right partitions are below the threshold size. If so, *
         * ignore one or both.  Otherwise, push the larger partition's    *
        \* bounds on the stack and continue sorting the smaller one.      */

        if ( (right_ptr - lo) <= MAX_THRESH ) {
          if ((hi - left_ptr) <= MAX_THRESH ) {
            --top;
            if ( top < stack ) break;
            lo = top -> lo;
            hi = top -> hi;
          } else {
            lo = left_ptr;
          }
        } else if ((hi - left_ptr) <= MAX_THRESH) {
          hi = right_ptr;
        } else if ((right_ptr - lo) > (hi - left_ptr)) {
          top -> lo = lo;
          top -> hi = right_ptr;
          ++top;
          lo = left_ptr;
        } else {
          top -> lo = left_ptr;
          top -> hi = hi;
          ++top;
          hi = right_ptr;
        }
      }

      /* Once the BASE_PTR array is partially sorted by quicksort the rest  *\
       * is completely sorted using insertion sort, since this is efficient *
       * for partitions below MAX_THRESH size.                              */
    }

  insert_sort:

    for ( integer i = 1; i < total_elems; ++i ) {
      for ( integer j = i; j > 0; --j ) {
        if ( GT(I[j], J[j], I[j-1], J[j-1]) ) break;
        swap( I[j], I[j-1] );
        swap( J[j], J[j-1] );
        swap( A[j], A[j-1] );
      }
    }
  }
  //@}

  /*
  //  ####   #####     ##    #####    ####   ######
  // #       #    #   #  #   #    #  #       #
  //  ####   #    #  #    #  #    #   ####   #####
  //      #  #####   ######  #####        #  #
  // #    #  #       #    #  #   #   #    #  #
  //  ####   #       #    #  #    #   ####   ######
  */

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  //!
  //! Macro to define the operator overloading of
  //! standard matrix-vector expressions.
  //!
  #define SPARSELIB_MUL_STRUCTURES(MATRIX)                                    \
  template <typename TM, typename VEC> inline                                 \
  Vector_M_mul_V<MATRIX<TM>,VEC>                                              \
  operator * (MATRIX<TM> const & M, VEC const & a) {                          \
    return Vector_M_mul_V<MATRIX<TM>,VEC>(M,a);                               \
  }                                                                           \
                                                                              \
  template <typename S, typename TM> inline                                   \
  Vector_S_mul_M<S,MATRIX<TM> >                                               \
  operator * (S const & s, MATRIX<TM> const & M) {                            \
    return Vector_S_mul_M<S,MATRIX<TM> >(s,M);                                \
  }                                                                           \
                                                                              \
  template <typename S, typename TM, typename VEC> inline                     \
  Vector_S_mul_M_mul_V<S,MATRIX<TM>, VEC>                                     \
  operator * (Vector_S_mul_M<S,MATRIX<TM> > const & sM,                       \
              VEC                           const & v) {                      \
    return Vector_S_mul_M_mul_V<S,MATRIX<TM>,VEC>(sM.s,sM.M,v);               \
  }                                                                           \
                                                                              \
  template <typename S, typename TM, typename VEC> inline                     \
  Vector_S_mul_M_mul_V<S,MATRIX<TM>,VEC>                                      \
  operator * (S const & s, Vector_M_mul_V<MATRIX<TM>,VEC> const & Mv) {       \
    return Vector_S_mul_M_mul_V<S,MATRIX<TM>,VEC>(s,Mv.M,Mv.a);               \
  }                                                                           \
                                                                              \
  template <typename TM, typename VA, typename VB> inline                     \
  Vector_V_sum_M_mul_V<VA,MATRIX<TM>,VB>                                      \
  operator + (VA const & a, Vector_M_mul_V<MATRIX<TM>,VB> const & Mv ) {      \
    return Vector_V_sum_M_mul_V<VA,MATRIX<TM>,VB>(a,Mv.M,Mv.a);               \
  }                                                                           \
                                                                              \
  template <typename TM, typename VA, typename VB> inline                     \
  Vector_V_sub_M_mul_V<VA,MATRIX<TM>,VB>                                      \
  operator - ( VA const & a, Vector_M_mul_V<MATRIX<TM>, VB> const & Mv ) {    \
    return Vector_V_sub_M_mul_V<VA,MATRIX<TM>,VB>(a,Mv.M,Mv.a);               \
  }                                                                           \
                                                                              \
  template <typename VA, typename S, typename TM, typename VB> inline         \
  Vector_V_sum_S_mul_M_mul_V<VA,S,MATRIX<TM>,VB>                              \
  operator + ( VA const & a,                                                  \
               Vector_S_mul_M_mul_V<S,MATRIX<TM>,VB> const & sMv ) {          \
    return Vector_V_sum_S_mul_M_mul_V<VA,S,MATRIX<TM>,VB>(a,sMv.s,sMv.M,sMv.a); \
  }                                                                           \
                                                                              \
  template <typename TM, typename VEC> inline                                 \
  Vector_Mt_mul_V<MATRIX<TM>, VEC>                                            \
  operator ^ ( MATRIX<TM> const & M, VEC const & a ) {                        \
    return Vector_Mt_mul_V<MATRIX<TM>,VEC>(M,a);                              \
  }                                                                           \
                                                                              \
  template <typename S, typename TM, typename VEC> inline                     \
  Vector_S_mul_Mt_mul_V<S,MATRIX<TM>,VEC>                                     \
  operator ^ (Vector_S_mul_M<S,MATRIX<TM> > const & sM, VEC const & v) {      \
    return Vector_S_mul_Mt_mul_V<S,MATRIX<TM>,VEC>(sM.s,sM.M,v);              \
  }                                                                           \
                                                                              \
  template <typename S, typename TM, typename VEC> inline                     \
  Vector_S_mul_Mt_mul_V<S,MATRIX<TM>,VEC>                                     \
  operator * (S const & s,                                                    \
              Vector_Mt_mul_V<MATRIX<TM>,VEC> const & Mv ) {                  \
    return Vector_S_mul_Mt_mul_V<S,MATRIX<TM>,VEC>(s,Mv.M,Mv.a);              \
  }                                                                           \
                                                                              \
  template <typename VA, typename TM, typename VB> inline                     \
  Vector_V_sum_Mt_mul_V<VA,MATRIX<TM>,VB>                                     \
  operator + ( VA const & a,                                                  \
               Vector_Mt_mul_V<MATRIX<TM>,VB> const & Mv ) {                  \
    return Vector_V_sum_Mt_mul_V<VA,MATRIX<TM>,VB>(a,Mv.M,Mv.a);              \
  }                                                                           \
                                                                              \
  template <typename VA, typename TM, typename VB> inline                     \
  Vector_V_sub_Mt_mul_V<VA,MATRIX<TM>,VB>                                     \
  operator - ( VA const & a,                                                  \
               Vector_Mt_mul_V<MATRIX<TM>,VB> const & Mv ) {                  \
    return Vector_V_sub_Mt_mul_V<VA,MATRIX<TM>,VB>(a,Mv.M,Mv.a);              \
  }                                                                           \
                                                                              \
  template <typename VA, typename S, typename TM, typename VB> inline         \
  Vector_V_sum_S_mul_Mt_mul_V<VA,S,MATRIX<TM>,VB>                             \
  operator + ( VA const & a,                                                  \
               Vector_S_mul_Mt_mul_V<S,MATRIX<TM>,VB> const & sMv ) {         \
    return Vector_V_sum_S_mul_Mt_mul_V<VA,S,MATRIX<TM>,VB>(a,sMv.s,sMv.M,sMv.a); \
  }

  #endif

  //!
  //! This class in the base class for all the
  //! sparse matrix classes of `Sparse_tool`.
  //! This class is incomplete and is used as a pivot for internal operations.
  //!
  template <typename Matrix>
  class SparseBase {
  protected:

    integer sp_nrows;     //!< Number of rows of the derived class
    integer sp_ncols;     //!< Number of columns of the derived class
    integer sp_min_size;  //!< Minimum between `sp`_nrows and `sp`_ncols
    integer sp_max_size;  //!< Minimum between `sp`_nrows and `sp`_ncols
    integer sp_nnz;       //!< Total number of nonzeros of the derived sparse matrix
    integer sp_lower_nnz; //!< Total number of nonzeros under the main diagonal
    integer sp_upper_nnz; //!< Total number of nonzeros over the main diagonal
    integer sp_diag_nnz;  //!< Total number of nonzeros on the main diagonal
    bool      sp_is_ordered;
    /*!< \brief Some sparse matrix can be internally in a state not ordered.
         When in this state the random access to element is unpredictable
         and can result in runtime error.
         This method return `true` if the matrix is ordered. */

    //!
    //! Check the index `idx`. If out of the range `0..nnz-1`
    //! and error is issued.
    //!
    void
    test_nnz(integer idx) const {
      UTILS_ASSERT(
        idx < sp_nnz,
        "Sparse::operator [{}] index out of range\n", idx
      );
    }

    //!
    //! Check the indices `i` and `j`.
    //! If out of the range of the matrix and error is issued.
    //!
    void
    test_index(integer i, integer j) const {
      UTILS_ASSERT(
        i < sp_nrows && j < sp_ncols,
        "Sparse::test_index({},{}) index out of range\n", i, j
      );
    }

    //!
    //! Check the indices `i`.
    //! If out of the row range of the matrix and error is issued.
    //!
    void
    test_row(integer i) const {
      UTILS_ASSERT(
        i < sp_nrows,
        "Sparse::test_row({}) index out of range\n", i
      );
    }

    //!
    //! Check the indices `j`.
    //! If out of the column range of the matrix and error is issued.
    //!
    void
    test_col(integer j) const {
      UTILS_ASSERT(
        j < sp_ncols,
        "Sparse::test_col({}) index out of range\n", j
      );
    }

    //!
    //! Initialize the class `Sparse` with `nr` rows and `nc` columns.
    //!
    void
    setup( integer nr, integer nc ) {
      sp_nrows      = nr;
      sp_ncols      = nc;
      sp_min_size   = std::min(nr,nc);
      sp_max_size   = std::max(nr,nc);
      sp_nnz        = 0;
      sp_lower_nnz  = 0;
      sp_upper_nnz  = 0;
      sp_diag_nnz   = 0;
      sp_is_ordered = true;
    }

    //!
    //! Update the counter for the lower, upper and diagonal elements.
    //!
    void
    ldu_count(integer i, integer j) {
      if      ( j < i ) ++sp_lower_nnz;
      else if ( j > i ) ++sp_upper_nnz;
      else              ++sp_diag_nnz;
    }

  public:

    SparseBase(void) { setup(0, 0); }
    ~SparseBase(void) {}

    // common data
    integer nrows    (void) const { return sp_nrows; }    //!< return the number of rows of the `Sparse` object
    integer ncols    (void) const { return sp_ncols; }    //!< return the number of columns of the `Sparse` object
    integer min_size (void) const { return sp_min_size; } //!< return the minimum between rows and columns
    integer max_size (void) const { return sp_max_size; } //!< return the maximum between rows and columns
    integer nnz      (void) const { return sp_nnz; }      //!< return the total number of nonzeros

    // ALIAS
    #ifdef LAPACK_WRAPPER_USE_ALIAS
    integer numRows (void) const { return nrows(); }    //!< return the number of rows of the `Sparse` object
    integer numCols (void) const { return ncols(); }    //!< return the number of columns of the `Sparse` object
    integer minSize (void) const { return min_size(); } //!< return the minimum between rows and columns
    integer maxSize (void) const { return max_size(); } //!< return the maximum between rows and columns
    #endif

    //!
    //! Return the number of nonzeros of the derived sparse
    //! matrix under the main diagonal (**lower**)
    //! over the main diagonal (**upper**) and on the diagonal (**diag**)
    //!
    void
    nnz( integer & lower, integer & diag, integer & upper) const {
      lower = sp_lower_nnz;
      diag  = sp_diag_nnz;
      upper = sp_upper_nnz;
    }
    //!
    //! Return `true` if the internal data is ordered.
    //!
    bool is_ordered (void) const { return sp_is_ordered; }

    //!
    //! \name Iterators
    //!
    //! These methods are useful for accessing all the nonzero elements
    //! of the derived sparse matrix.  The methods
    //!
    //! - `void Begin()`
    //! - `void Next()`
    //! - `bool End()`
    //!
    //! permits to loops on all the elements, while the methods
    //!
    //! - `integer row()`
    //! - `integer column()`
    //! - `real_type value()`
    //!
    //! permits to access values of the actual elements
    //! pointed by the iterator. For example to print all the stored
    //! values of the `Sparse` object \c S we can do:
    //!
    //! \code
    //!   for ( S.Begin(); S.End(); S.Next() ) {
    //!     cout << " row   = " << S.row()
    //!          << " col   = " << S.column()
    //!          << " value = " << S.value()
    //!          << '\n';
    //!   }
    //! \endcode
    //!
    //@{
    //!
    //! Set the iterator at the begin of the loop.
    //!
    void Begin(void) const { static_cast<Matrix const *>(this) -> Begin(); }
    //!
    //! Go the next item.
    //!
    void Next(void) const { static_cast<Matrix const *>(this) -> Next(); }
    //!
    //! Return `true` unless we are at the end of the loop.
    //!
    bool End(void) const { return static_cast<Matrix const *>(this) -> End(); }

    //!
    //! The **row** of the pointed element.
    //!
    integer row(void) const { return static_cast<Matrix const *>(this) -> row(); }

    //!
    //! The **column** of the pointed element.
    //!
    integer column(void) const { return static_cast<Matrix const *>(this) -> column(); }

    //!
    //! Assign the pointed element to `rhs`.
    //!
    template <typename T>
    void assign( T & rhs ) const {
      Matrix const * This = static_cast<Matrix const *>(this);
      This -> assign(rhs);
    }

    //@}

    //!
    //! Return `true` if the element `(i,j)` is
    //! a nonzeros of the derived matrix.
    //!
    bool
    exists( integer i, integer j )
    { return static_cast<Matrix const *>(this) -> exists( i, j ); }

    //!
    //! Return the position of the `(i,j)` elements
    //! in the internal structure. If the element `(i,j)` do not
    //! exist it return `nnz()`
    //!
    integer
    position( integer i, integer j ) const
    { return static_cast<Matrix const *>(this) -> position( i, j ); }

  };

  //!
  //! This class in the base class for all the
  //! sparse matrix classes of `Sparse_tool`.
  //! This class is incomplete and is used
  //! as a pivot for internal operations.
  //!

  template <typename T, typename Matrix>
  class Sparse : public SparseBase<Matrix> {
    using SBASE = SparseBase<Matrix>;
  public:
    using real_type = T; //!< type of the element of the matrix

    Sparse(void) : SparseBase<Matrix>() { SBASE::setup(0, 0); }
    ~Sparse(void) {}

    //!
    //! The **value** of the pointed element.
    //!
    real_type value(void) const { return static_cast<Matrix const *>(this) -> value(); }

    //!
    //! Assign the pointed element to `rhs`.
    //!
    template <typename TS>
    void assign( TS & rhs ) const { rhs = this -> value(); }

    //!
    //! \name Random access to elements.
    //!
    //! The following methods permits random access to
    //! nonzeros elements of the derived class.
    //!
    //@{
    //!
    //! Access to the elements of the sparse matrix in
    //! a random way.  For example given a sparse matrix `A`  the code
    //! `A(i,j)` return the **value** of the elements of the
    //! matrix at the `i-th` row and `j`-th column.
    //! If at the `(i,j)` coordinate there are no elements
    //! an error is issued.
    //!
    real_type const &
    operator () ( integer i, integer j ) const
    { return static_cast<Matrix const *>(this) -> operator () (i,j); }

    //!
    //! Access to the elements of the sparse matrix in
    //! a random way.  For example given a sparse matrix `A`  the code
    //! `A(i,j)` return the **value** of the elements of the
    //! matrix at the `i-th` row and `j`-th column.
    //! If at the `(i,j)` coordinate there are no elements an error is issued.
    //!
    real_type &
    operator () ( integer i, integer j )
    { return static_cast<Matrix *>(this) -> operator () (i,j); }

    //!
    //! This method permits to access the elements of the sparse matrix in
    //! a random way.  For example given a sparse matrix `A`  the code
    //! `A(i,j)` return the **reference** of the elements of the
    //! matrix at the `i-th` row and `j-th` column.  If at the `(i,j)`
    //! coordinate there are no elements a reference to a value **0** is returned
    //!
    real_type const &
    value( integer i, integer j ) const
    { return static_cast<Matrix const *>(this) -> value(i,j); }

    //!
    //! Return `true` if the element `(i,j)`
    //! is a nonzeros of the derived matrix.
    //!
    bool
    exists( integer i, integer j )
    { return static_cast<Matrix const *>(this) -> exists( i, j ); }

    //!
    //! Return the position of the `(i,j)` elements
    //! in the internal structure. If the element `(i,j)` do not
    //! exist it return `nnz()`.
    //!
    integer
    position( integer i, integer j ) const
    { return static_cast<Matrix const *>(this) -> position( i, j ); }

    //!
    //! Return the reference of the `idx-th` element
    //! of the vector of stored values.
    //!
    real_type const &
    operator [] (integer idx) const
    { return static_cast<Matrix const *>(this) -> operator[] (idx); }

    //!
    //! Return the reference of the `idx-th`
    //! element of the vector of stored values.
    //!
    real_type &
    operator [] (integer idx)
    { return static_cast<Matrix *>(this) -> operator[] (idx); }

    //@}

    //!
    //! Set all the nonzeros of the sparse matrix to the value `0`.
    //!
    void
    setZero()
    { return static_cast<Matrix const *>(this) -> setZero(); }

    //!
    //! Multiply all the nonzeros of the sparse matrix by `s`.
    //!
    void
    scale_values( real_type const & s )
    { return static_cast<Matrix const *>(this) -> scale_values( s ); }

    //!
    //! Multiply all the nonzeros of the row `nr` matrix by `val`.
    //!
    void
    scaleRow( integer nr, real_type const & val )
    { return static_cast<Matrix const *>(this) -> scaleRow( nr, val ); }

    //!
    //! Multiply all the nonzeros of the column `nc` matrix by `val`.
    //!
    void
    scaleColumn( integer nc, real_type const & val )
    { return static_cast<Matrix const *>(this) -> scaleColumn( nc, val ); }

    //!
    //! \name Diagonal internal operation
    //!
    //@{
    //!
    //! Set `s` to all the components of the diagonal,
    //! and `0` the others elements.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    & 1  &    \\
    //!     1 & 2  &    &    & 5  \\
    //!       &    & 2  & -3 &    \\
    //!       &    & 3  & 2  & -4 \\
    //!       &    &    & 4  & 2
    //! \end{pmatrix},\qquad
    //! (\mathbf{A}=2.1) = \begin{pmatrix}
    //!    2.1 & 0   &     & 0   &    \\
    //!    0   & 2.1 &     &     & 0  \\
    //!        &     & 2.1 & 0   &    \\
    //!        &     & 0   & 2.1 & 0  \\
    //!        &     &     & 0   & 2.1
    //! \end{pmatrix}
    //! \f]
    //!
    Matrix &
    operator = ( real_type const & s ) {
      setZero();
      if ( s != 0 )
        for ( integer k = 0; k < SBASE::sp_min_size; ++k )
          (*this)(k,k) = s;
      return *this;
    }

    //!
    //! Add `s` to all the components of the diagonal.
    //! For example:
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //! 2 & -1 &    & 1  &    \\
    //! 1 & 2  &    &    & 5  \\
    //!   &    & 2  & -3 &    \\
    //!   &    & 3  & 2  & -4 \\
    //!   &    &    & 4  & 2
    //! \end{pmatrix},\qquad
    //! (\mathbf{A}\verb|+=| 2) = \begin{pmatrix}
    //! 4 & -1 &    & 1  &    \\
    //! 1 & 4  &    &    & 5  \\
    //!   &    & 4  & -3 &    \\
    //!   &    & 3  & 4  & -4 \\
    //!   &    &    & 4  & 4
    //! \end{pmatrix}
    //! \f]
    //!
    Matrix &
    operator += ( real_type const & s ) {
      for ( integer k = 0; k < SBASE::sp_min_size; ++k )
        (*this)(k,k) += s;
      return *this;
    }

    //!
    //! Subtract `s` to all the components of the diagonal.
    //! For example:
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //! 2 & -1 &    & 1  &    \\
    //! 1 & 2  &    &    & 5  \\
    //!   &    & 2  & -3 &    \\
    //!   &    & 3  & 2  & -4 \\
    //!   &    &    & 4  & 2
    //! \end{pmatrix},\qquad
    //! (\mathbf{A} \verb|-=| 2) = \begin{pmatrix}
    //! 0 & -1 &    & 1  &    \\
    //! 1 & 0  &    &    & 5  \\
    //!   &    & 0  & -3 &    \\
    //!   &    & 3  & 0  & -4 \\
    //!   &    &    & 4  & 0
    //! \end{pmatrix}
    //! \f]
    //!
    Matrix &
    operator -= ( real_type const & s ) {
      for ( integer k = 0; k < SBASE::sp_min_size; ++k )
        (*this)(k,k) -= s;
      return *this;
    }

    //!
    //! Multiply by `s` to all the components of the diagonal.
    //! For example:
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //! 2 & -1 &    & 1  &    \\
    //! 1 & 2  &    &    & 5  \\
    //!   &    & 2  & -3 &    \\
    //!   &    & 3  & 2  & -4 \\
    //!   &    &    & 4  & 2
    //! \end{pmatrix},\qquad
    //! (\mathbf{A} \verb|*=| 2) = \begin{pmatrix}
    //! 4 & -2 &    & 2  &    \\
    //! 2 &  4 &    &    & 10 \\
    //!   &    & 4  & -6 &    \\
    //!   &    & 6  & 4  & -8 \\
    //!   &    &    & 8  & 4
    //! \end{pmatrix}
    //! \f]
    //!
    Matrix &
    operator *= ( real_type const & s )
    { scale_values(s); return * this; }

    //!
    //! Divide all the nonzeros of the derived matrix by `s`.
    //! For example:
    //!
    //! \f[
    //! \mathbf{A}= \begin{pmatrix}
    //! 2 & -1 &    & 1  &    \\
    //! 1 & 2  &    &    & 5  \\
    //!   &    & 2  & -3 &    \\
    //!   &    & 3  & 2  & -4 \\
    //!   &    &    & 4  & 2
    //! \end{pmatrix},\qquad
    //! (\mathbf{A} \verb|/=| 10) = \begin{pmatrix}
    //! 0.2 & -0.1 &      & 0.1  &      \\
    //! 0.1 & 0.2  &      &      & 0.5  \\
    //!     &      & 0.2  & -0.3 &      \\
    //!     &      & 0.3  & 0.2  & -0.4 \\
    //!     &      &      & 0.4  & 0.2
    //! \end{pmatrix}
    //! \f]
    //!
    Matrix &
    operator /= ( real_type const & s)
    { real_type ss = 1/s; scale_values(ss); return * this; }

    // ------------------------
    //!
    //! Set the diagonal values of the sparse
    //! matrix to `v`. For example
    //!
    //! \f[
    //!     \mathbf{A}= \begin{pmatrix}
    //!     2 & -1 &    & 1  &    \\
    //!     1 & 2  &    &    & 5  \\
    //!       &    & 2  & -3 &    \\
    //!       &    & 3  & 2  & -4 \\
    //!       &    &    & 4  & 2
    //!     \end{pmatrix},\qquad
    //!     \mathbf{v} = \begin{pmatrix}
    //!     1 \\ 2\\ 3
    //!     \end{pmatrix},
    //!     \qquad
    //!     (\mathbf{A}=\mathbf{v}) = \begin{pmatrix}
    //!     1 & 0 &   & 0 &   \\
    //!     0 & 2 &   &   & 0 \\
    //!       &   & 3 & 0 &   \\
    //!       &   & 0 & 0 & 0 \\
    //!       &   &   & 0 & 0
    //!    \end{pmatrix}
    //! \f]
    //!
    //! Notice that only the first `v.size()` diagonal elements are
    //! set while the rest of the matrix is set to `0`.
    //!
    template <typename TV> inline
    Matrix &
    operator = ( Vector<TV> const & v ) {
      integer n = minIndex(v.size(), SBASE::sp_min_size );
      setZero();
      for ( integer k = 0; k < n; ++k )
        (*this)(k,k) = v(k);
      return *this;
    }

    //!
    //! Add the components of `v`  to the diagonal.
    //! For example:
    //!
    //! \f[
    //!   \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    & 1  &    \\
    //!     1 & 2  &    &    & 5  \\
    //!       &    & 2  & -3 &    \\
    //!       &    & 3  & 2  & -4 \\
    //!       &    &    & 4  & 2
    //!   \end{pmatrix},\qquad
    //!   \mathbf{v} = \begin{pmatrix}
    //!     1 \\ 2\\ 3
    //!   \end{pmatrix},
    //!   \qquad
    //!   (\mathbf{A}\verb|-=| \mathbf{v}) = \begin{pmatrix}
    //!     3 & -1 &    & 1  &    \\
    //!     1 & 4  &    &    & 5  \\
    //!       &    & 5  & -3 &    \\
    //!       &    & 3  & 2  & -4 \\
    //!       &    &    & 4  & 2
    //!   \end{pmatrix}
    //! \f]
    //!
    template <typename TV> inline
    Matrix &
    operator += ( Vector<TV> const & v ) {
      integer n = minIndex( v.size(), SBASE::sp_min_size );
      for ( integer k = 0; k < n; ++k )
        (*this)(k,k) += v(k);
      return *this;
    }

    //!
    //! Subtract the components of `v`  to the diagonal.
    //! For example:
    //!
    //! \f[
    //!  \mathbf{A} = \begin{pmatrix}
    //!    2 & -1 &    & 1  &    \\
    //!    1 & 2  &    &    & 5  \\
    //!      &    & 2  & -3 &    \\
    //!      &    & 3  & 2  & -4 \\
    //!      &    &    & 4  & 2
    //!  \end{pmatrix},\qquad
    //!  \mathbf{v} = \begin{pmatrix}
    //!    1 \\ 2\\ 3
    //!  \end{pmatrix},
    //!  \qquad
    //!  (\mathbf{A}\verb|-=| \mathbf{v}) = \begin{pmatrix}
    //!    1 & -1 &    & 1  &    \\
    //!    1 & 0  &    &    & 5  \\
    //!      &    & -1  & -3 &   \\
    //!      &    & 3  & 2  & -4 \\
    //!      &    &    & 4  & 2
    //!  \end{pmatrix}
    //! \f]
    //!
    template <typename TV> inline
    Matrix &
    operator -= ( Vector<TV> const & v ) {
      integer n = minIndex(v.size(),  SBASE::sp_min_size);
      for ( integer k = 0; k < n; ++k )
        (*this)(k,k) -= v(k);
      return *this;
    }

    //!
    //! Assign the vector expression `e` to the components of the diagonal,
    //! and `0` the others elements.  For example:
    //!
    //! \f[
    //!  \mathbf{A} = \begin{pmatrix}
    //!      2 & -1 &    & 1  &    \\
    //!      1 & 2  &    &    & 5  \\
    //!        &    & 2  & -3 &    \\
    //!        &    & 3  & 2  & -4 \\
    //!        &    &    & 4  & 2
    //!  \end{pmatrix},\qquad
    //!  \mathbf{a} = \begin{pmatrix}
    //!      1 \\ 2 \\ 3 \\ 4 \\ 5
    //!  \end{pmatrix},\qquad
    //!  \mathbf{b} = \begin{pmatrix}
    //!      1 \\ -1 \\ 2
    //!  \end{pmatrix},
    //! \f]
    //!
    //! \f[
    //!  (\mathbf{A}=\mathbf{a}+2\mathbf{b}) = \begin{pmatrix}
    //!     3 & 0 &   & 0  &    \\
    //!     0 & 0 &   &    & 0  \\
    //!       &   & 7 & 0  &    \\
    //!       &   & 0 & 0  & 0  \\
    //!       &   &   & 0  & 0
    //!   \end{pmatrix}
    //! \f]
    //!
    template <typename VEC_EXPR> inline
    Matrix &
    operator = ( VEC_EXPR const & e ) {
      integer n = minIndex(e.size(), SBASE::sp_min_size );
      setZero();
      for ( integer k = 0; k < n; ++k )
        (*this)(k,k) = e(k);
      return *this;
    }

    //!
    //! Add the vector expression `e` to the components of the diagonal,
    //! and `0` the others elements.  For example:
    //!
    //! \f[
    //!  \mathbf{A} = \begin{pmatrix}
    //!      2 & -1 &    & 1  &    \\
    //!      1 & 2  &    &    & 5  \\
    //!        &    & 2  & -3 &    \\
    //!        &    & 3  & 2  & -4 \\
    //!        &    &    & 4  & 2
    //!  \end{pmatrix},\qquad
    //!  \mathbf{a} = \begin{pmatrix}
    //!      1 \\ 2 \\ 3 \\ 4 \\ 5
    //!  \end{pmatrix},\qquad
    //!  \mathbf{b} = \begin{pmatrix}
    //!      1 \\ -1 \\ 2
    //!  \end{pmatrix},
    //! \f]
    //!
    //! \f[
    //!  (\mathbf{A}\verb|+=| \mathbf{a}+2\mathbf{b}) = \begin{pmatrix}
    //!      5 & -1 &    & 1  &    \\
    //!      1 & 2  &    &    & 5  \\
    //!        &    & 9  & -3 &    \\
    //!        &    & 3  & 2  & -4 \\
    //!        &    &    & 4  & 2
    //!  \end{pmatrix}
    //! \f]
    //!
    template <typename VEC_EXPR> inline
    Matrix &
    operator += ( VEC_EXPR const & e ) {
      integer n = minIndex(e.size(), SBASE::sp_min_size );
      for ( integer k = 0; k < n; ++k )
        (*this)(k,k) += e(k);
      return *this;
    }

    //!
    //! Add the vector expression `e` to the components of the diagonal,
    //! and `0` the others elements.  For example:
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    & 1  &    \\
    //!     1 & 2  &    &    & 5  \\
    //!       &    & 2  & -3 &    \\
    //!       &    & 3  & 2  & -4 \\
    //!       &    &    & 4  & 2
    //! \end{pmatrix},\qquad
    //! \mathbf{a} = \begin{pmatrix}
    //!     1 \\ 2 \\ 3 \\ 4 \\ 5
    //! \end{pmatrix},\qquad
    //! \mathbf{b} = \begin{pmatrix}
    //!     1 \\ -1 \\ 2
    //! \end{pmatrix},
    //! \f]
    //!
    //! \f[
    //! (\mathbf{A} \verb|-=| \mathbf{a}+2\mathbf{b}) = \begin{pmatrix}
    //!    -1 & -1 &    & 1  &    \\
    //!     1 & 2  &    &    & 5  \\
    //!       &    & -5 & -3 &    \\
    //!       &    & 3  & 2  & -4 \\
    //!       &    &    & 4  & 2
    //! \end{pmatrix}
    //! \f]
    //!
    template <typename VEC_EXPR>
    inline
    Matrix &
    operator -= ( VEC_EXPR const & e ) {
      integer n = minIndex(e.size(),  SBASE::sp_min_size );
      for ( integer k = 0; k < n; ++k )
        (*this)(k,k) -= e(k);
      return *this;
    }
    //@}

    //!
    //! \name Matrix-Vector multiplication
    //!
    //@{

    //!
    //! Perform the operation `res += s * (A * x)`
    //!
    template <typename VEC>
    void
    add_S_mul_M_mul_V(
      VEC       & res,
      T   const & s,
      VEC const & x
    ) const {
      UTILS_ASSERT0(
        std::addressof(res) != std::addressof(x),
        "Sparse_tool: [res += s * (A * x)]: add_S_mul_M_mul_V `res` and `x` cant be the same"
      );
      UTILS_ASSERT0(
        res.size() >= SBASE::sp_nrows,
        "Sparse_tool: [res += s * (A * x)] result vector too small"
      );
      for ( SBASE::Begin(); SBASE::End(); SBASE::Next() )
        res(SBASE::row()) += s * value() * x(SBASE::column());
    }

    //!
    //! perform the operation `res += s * (A^T * x)`
    //!
    template <typename VEC>
    void
    add_S_mul_Mt_mul_V(
      VEC       & res,
      T   const & s,
      VEC const & x
    ) const {
      UTILS_ASSERT0(
        std::addressof(res) != std::addressof(&x),
        "Sparse_tool: [res += s * (A^T * x)] add_S_mul_M_mul_V, `res` and `x` cant be the same"
      );
      UTILS_ASSERT0(
        res.size() >= SBASE::sp_ncols,
        "Sparse_tool: [res += s * (A^T * x)] result vector too small"
      );
      for ( SBASE::Begin(); SBASE::End(); SBASE::Next() )
        res(SBASE::column()) += s * value() * x(SBASE::row());
    }

    //!
    //! Perform the operation `res += s * (A * x)`
    //!
    template <typename VA, typename VB>
    void
    add_S_mul_M_mul_V(
      VA       & res,
      T  const & s,
      VB const & x
    ) const {
      UTILS_ASSERT0(
        res.size() >= SBASE::sp_nrows,
        "Sparse_tool: [res += s * (A * x)] result vector too small"
      );
      for ( SBASE::Begin(); SBASE::End(); SBASE::Next() )
        res(SBASE::row()) += s * value() * x(SBASE::column());
    }

    //!
    //! perform the operation `res += s * (A^T * x)`
    //!
    template <typename VA, typename VB>
    void
    add_S_mul_Mt_mul_V(
      VA       & res,
      T  const & s,
      VB const & x
    ) const {
      UTILS_ASSERT0(
        res.size() >= SBASE::sp_ncols,
        "Sparse_tool: [res += s * (A^T * x)] result vector too small"
      );
      for ( SBASE::Begin(); SBASE::End(); SBASE::Next() )
        res(SBASE::column()) += s * value() * x(SBASE::row());
    }

    //@}
  };

  //!
  //! This class in the base class for all the preconditioner
  //! of sparse matrix classes of `Sparse_tool`.
  //! This class is incomplete and is used as a
  //! pivot for internal operations.
  //!
  template <typename PRECO>
  class Preco {
  protected:
    //! The number of row/column of the preconditioner matrix.
    integer pr_size;
  public:
    //!
    //! Create an empty preconditioner.
    //!
    Preco(void) : pr_size(0) {}
    ~Preco(void) {}
    //!
    //! Return the number of row/column of the preconditioner matrix.
    //!
    integer size(void) const { return pr_size; }
  };

  /*
  //  #####
  // #     #  #####     ##    #####    ####   ######
  // #        #    #   #  #   #    #  #       #
  //  #####   #    #  #    #  #    #   ####   #####
  //       #  #####   ######  #####        #  #
  // #     #  #       #    #  #   #   #    #  #
  //  #####   #       #    #  #    #   ####   ######
  //
  // ######
  // #     #    ##    #####  #####  ######  #####   #    #
  // #     #   #  #     #      #    #       #    #  ##   #
  // ######   #    #    #      #    #####   #    #  # #  #
  // #        ######    #      #    #       #####   #  # #
  // #        #    #    #      #    #       #   #   #   ##
  // #        #    #    #      #    ######  #    #  #    #
  */

  //!
  //! Sparse Pattern in compressed coordinate
  //! ---------------------------------------
  //!
  //! The sparse pattern internal structure is essentially a sparse
  //! compressed coordinate ones.  It consists of two big vector of
  //! integer which contain the coordinate of nonzero elements.  We call
  //! `I` the vector that store the first coordinate, while we call
  //! `J` the vector that store the second coordinate.  For example the
  //! following `6 x 7` sparse matrix pattern
  //!
  //! \f[
  //!   \mathbf{S} = \begin{pmatrix}
  //!  * & * &   &   &   &   & * \\
  //!    & * & * &   &   &   &   \\
  //!    &   &   & * & * &   &   \\
  //!    & * &   &   & * &   &   \\
  //!  * &   & * &   &   & * &   \\
  //!    & * &   &   &   &   & *
  //!    \end{pmatrix}
  //! \f]
  //!
  //! \f[
  //!    \begin{array}{c|cccccccccccccc}
  //!      I & 0 & 0 & 0 & 1 & 1 & 2 & 2 & 3 & 3 & 4 & 4 & 4 & 5 & 5 \\
  //!      \hline
  //!      J & 0 & 1 & 6 & 1 & 2 & 3 & 4 & 1 & 4 & 0 & 2 & 5 & 1 & 6
  //!    \end{array}
  //! \f]
  //!
  //! Notice that the index follows the C convention, starting from **0**.
  //! The class `SparsePattern` try to manage such a structure in a
  //! simple way for the user.  To define a `SparsePattern` class you
  //! can use uno of the following scripture
  //!
  //! \code
  //! SparsePatter sp;              // define an empty SparsePattern class
  //! SparsePatter sp(nr, nc, nnz); // define a SparsePattern class on a pattern of nr rows and nc columns.
  //!                               // The number nnz is an integer which determine the preallocated
  //!                               // number of pattern elements stored in the class, in practice is the
  //!                               // initial dimension of vectors I and J.
  //! SparsePatter sp(sp1);         // define a SparsePattern object which is the copy of SparsePattern object sp1.
  //! SparsePatter sp(sobj);        // define a SparsePattern object which is the pattern of nonzero of the
  //!                               // Sparse object sobj. In this way it is possible to obtain the sparse
  //!                               // pattern of any object derived from Sparse.
  //! \endcode
  //!
  //! It is possible in any moment to change the sizes and the maximum
  //! number of nonzero by the `resize` methods:
  //!
  //! \code
  //! sp.resize(nr, nc, nnz);
  //! sp.resize(sp1);
  //! sp.resize(sobj);
  //! \endcode
  //!
  //! to complete the basic description of the `SparsePattern` a sintetic
  //! description of the remaining methods of the class are presented
  //!
  //! - `insert(i,j)` this method permit to insert an item of nonzero.
  //!   For example the first pattern can be constructed as
  //!
  //! \code
  //!   SparsePattern sp(6,7);
  //!   sp.insert(0, 0); sp.insert(0, 1); sp.insert(0, 6); sp.insert(1, 1);
  //!   sp.insert(1, 2); sp.insert(2, 3); sp.insert(2, 4); sp.insert(3, 1);
  //!   sp.insert(3, 4); sp.insert(4, 0); sp.insert(4, 2); sp.insert(4, 5);
  //!   sp.insert(5, 1); sp.insert(5, 6); sp.internal_order();
  //! \endcode
  //!
  //! - `internal_order`()
  //!   This methods reorder internally the nonzero elements of the sparse pattern
  //!   and remove all duplicated entries.
  //!   If `k1 <= k2` we have one of the two following cases
  //!
  //! - `J(k1) < J(k2)`
  //! - `J(k1) == J(k2)` and `I(k1) <= I(k2)`
  //!
  //! - `bool is_ordered()`
  //!   this methods return `true` if the elements inside the sparse
  //!   pattern are ordered, `false` otherwise.
  //!
  class SparsePattern : public SparseBase<SparsePattern> {
    using MATRIX = SparsePattern;
    using SPARSE = SparseBase<SparsePattern>;

    Vector<integer> I;    //!< Vector of row index
    Vector<integer> J;    //!< Vector of column index
    mutable integer ipos; //!< Actual position of iterator

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    template <typename MAT, typename Compare>
    void
    convert( SparseBase<MAT> const & M, Compare cmp ) {
      // count nonzero
      integer nz{0};
      for ( M.Begin(); M.End(); M.Next() )
        if ( cmp( M.row(), M.column() ) ) ++nz;

      resize( M.nrows(), M.ncols(), nz );

      for ( M.Begin(); M.End(); M.Next() ) {
        integer i = M.row();
        integer j = M.column();
        if ( cmp(i,j) ) insert(i,j);
      }
      internal_order();
    }
    #endif

  public:

    //!
    //! Initialize and empty sparse pattern of `0` rows and `0` columns.
    //!
    SparsePattern(void) : SparseBase<SparsePattern>()
    { }

    //!
    //! Initialize and empty sparse pattern of `nr` rows and `nc` columns
    //! with reserved room for `mnnz` nonzeros.
    //!
    SparsePattern(
      integer nr,
      integer nc,
      integer mnnz = SPARSETOOL_DEFAULT_NNZ
    ) : SparsePattern()
    { resize(nr, nc, mnnz); }

    //!
    //! Assign the pointed element to `rhs`.
    //!
    template <typename T>
    void assign( T & ) const { /* do nothing! */ }

    //!
    //! Insert the element of the sparse pattern `sp`
    //! which satify `cmp` to `*this`.
    //!
    //! \param sp  sparse pattern to be copied
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename Compare>
    SparsePattern( SparsePattern const & sp, Compare cmp )
    : SparsePattern()
    { convert(sp,cmp); }

    //!
    //! Copy the sparse pattern `sp`to `*this`.
    //!
    //! \param sp sparse pattern to be copied
    //!
    SparsePattern(SparsePattern const & sp) : SparsePattern()
    { convert(sp,all_ok()); }

    // convert => SparsePattern
    //!
    //! Insert the element of the pattern of `M`.
    //! which satify `cmp` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    SparsePattern( SparseBase<MAT> const & M, Compare cmp )
    : SparsePattern()
    { convert(M,cmp); }

    //!
    //! Copy the sparse pattern of `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    SparsePattern( SparseBase<MAT> const & M )
    : SparsePattern()
    { convert(M,all_ok()); }

    //!
    //! Copy the sparse pattern of `sp`to `*this`.
    //!
    //! \param SP sparse pattern object
    //!
    SparsePattern &
    operator = ( SparsePattern const & SP ) {
      if ( &SP != this ) { // avoid copy to itself
        SPARSE::setup( SP.nrows(), SP.ncols() );
        SPARSE::sp_nnz        = SP.nnz();
        SPARSE::sp_is_ordered = SP.is_ordered();
        I.resize( SPARSE::sp_nnz );
        J.resize( SPARSE::sp_nnz );
        I = SP.I.head(SPARSE::sp_nnz);
        J = SP.J.head(SPARSE::sp_nnz);
        if ( !SPARSE::sp_is_ordered ) internal_order();
      }
      return *this;
    }

    // convert => SparsePattern
    //!
    //! Copy the sparse pattern of `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    SparsePattern &
    operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return *this; }

    //!
    //! Initialize and empty sparse pattern of `nr` rows and `nc` columns
    //! with reserved room for `mnnz` nonzeros.
    //!
    void
    resize( integer nr, integer nc, integer mnnz = SPARSETOOL_DEFAULT_NNZ ) {
      SPARSE::setup( nr, nc );
      I.resize( mnnz );
      J.resize( mnnz );
    }

    //!
    //! Initialize the pattern and insert the element of the pattern of `M`.
    //! to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    void
    resize(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    //!
    //! Initialize the pattern and insert the elements M(i,j)
    //! which satify `cmp(i,j)` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    void
    resize(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Order the internal structure to permit random
    //! access to nonzeros.
    //!
    void
    internal_order() {
      QuickSortIJ<integer>( I.data(), J.data(), SPARSE::sp_nnz );
      // eliminate duplicate elements
      integer i1{0}, i{0};
      SPARSE::sp_lower_nnz = 0;
      SPARSE::sp_diag_nnz  = 0;
      SPARSE::sp_upper_nnz = 0;
      while ( i1 < SPARSE::sp_nnz ) {
        // copy i1 => i
        I(i) = I(i1); J(i) = J(i1);
        // setup statistic
        SPARSE::ldu_count(I(i),J(i));
        // find i1 != i
        for ( ++i1; i1 < SPARSE::sp_nnz && I(i1) == I(i) && J(i1) == J(i); ++i1 ) {}
        ++i;
      }
      SPARSE::sp_nnz        = i;
      SPARSE::sp_is_ordered = true;
      I.resize( SPARSE::sp_nnz );
      J.resize( SPARSE::sp_nnz );
    }

    //!
    //! Check if pattern is empty.
    //!
    bool empty(void) const { return SPARSE::sp_nnz == 0; }

    //!
    //! Insert `(i,j)` in the sparse pattern.
    //!
    SparsePattern &
    insert( integer i, integer j ) {
      SPARSE::test_index(i,j);
      if ( I.size() <= SPARSE::sp_nnz ) {
        integer newdim = SPARSE::sp_nnz + (SPARSE::sp_nnz>>3) + 100;
        I.conservativeResize( newdim );
        J.conservativeResize( newdim );
      }
      I(SPARSE::sp_nnz) = i;
      J(SPARSE::sp_nnz) = j;
      ++SPARSE::sp_nnz;
      SPARSE::sp_is_ordered = false;
      return *this;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    //! Return the row index vector.
    Vector<integer> const & getI(void) const { return I; }
    //! Return the column index vector.
    Vector<integer> const & getJ(void) const { return J; }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // ITERATORS
    //! Initialize iterator.
    void Begin(void) const { ipos = 0; }

    //! Go to the next element.
    void Next(void) const { ++ipos; }

    //! Check if iterator terminate the loop.
    bool End(void) const { return ipos < SPARSE::sp_nnz; }

    //! The row of the element pointed by iterator.
    integer row(void) const { return I(ipos); }

    //! The column of the element pointed by iterator.
    integer column(void) const { return J(ipos); }
  };

  /*
  //  #####   #####
  // #     # #     #   ####    ####   #####
  // #       #        #    #  #    #  #    #
  // #       #        #    #  #    #  #    #
  // #       #        #    #  #    #  #####
  // #     # #     #  #    #  #    #  #   #
  //  #####   #####    ####    ####   #    #
  //
  // #     #
  // ##   ##    ##    #####  #####   #  #    #
  // # # # #   #  #     #    #    #  #   #  #
  // #  #  #  #    #    #    #    #  #    ##
  // #     #  ######    #    #####   #    ##
  // #     #  #    #    #    #   #   #   #  #
  // #     #  #    #    #    #    #  #  #    #
  //
  */

  //!
  //! Compressed Coordinate Matrix Storage
  //! ------------------------------------
  //!
  //! The class `CCoorMatrix<T>` implement a
  //! `Compressed Coordinate` storage sparse scheme.
  //! It consists of two big vector of integer which contain
  //! the coordinate of nonzero elements and a big one of real number
  //! which contain the values.  We call `I` the vector that store the
  //! rows coordinate, `J` the vector that store the columns coordinate
  //! and `A`  the vector that store the nonzero values.
  //! For example the following **6** x **7** sparse matrix
  //!
  //! \f[
  //! \mathbf{A} = \begin{pmatrix}
  //!  1 & 2    &    &   &   &   & 9 \\
  //!    & -1   & 0  &   &   &   &   \\
  //!    &      &    & 3 & 4 &   &   \\
  //!    & 2    &    &   & 5 &   &   \\
  //!  2 &      & -2 &   &   & 1 &   \\
  //!    & -1.5 &    &   &   &   & -1
  //! \end{pmatrix}
  //! \f]
  //! \f[
  //! \begin{array}{c|cccccccccccccc}
  //!    I & 0 & 0 & 0 & 1 & 1 & 2 & 2 & 3 & 3 & 4 & 4 & 4 &  5  & 5 \\
  //!    \hline
  //!    J & 0 & 1 & 6 & 1 & 2 & 3 & 4 & 1 & 4 & 0 & 2 & 5 &  1  & 6 \\
  //!    \hline
  //!    A & 1 & 2 & 9 &-1 & 0 & 3 & 4 & 2 & 5 & 2 &-2 & 1 &-1.5 &-1
  //! \end{array}
  //! \f]
  //!
  //! Notice that the index follows the C convention, starting from
  //! `0`.  The class `CCoorMatrix<T>` try to manage such a
  //! structure in a simple way for the user.  To define a
  //! `CCoorMatrix<T>` class you can use one of the following scripture
  //!
  //! \code
  //! CCoorMatrix<double> ccoor; // define an empty CCoorMatrix<double> class
  //!
  //! // instance a CCoorMatrix<double> class of nr rows and nc columns.
  //! // The number nnz is an integer which determine the pre-allocated number of elements
  //! // stored in the class, in practice it is the initial dimension of vectors `I`, `J` and `A`.
  //! //  If needed the vectors are enlarged.
  //! CCoorMatrix<double> ccoor(nr, nc, nnz);
  //!
  //! // instance a CCoorMatrix<double> class with the sparsity pattern defined of the SparsePattern class `sp`.
  //! // nr will be sp.nrows(), ncol will be sp.ncols().
  //! CCoorMatrix<double> ccoor(sp);
  //!
  //! // instance a CCoorMatrix<double> class with which is the copy of the CCoorMatrix<double> object ccoor1.
  //! CCoorMatrix<double> ccoor(ccoor1);
  //!
  //! // instance a CCoorMatrix<double> class with which is the copy of the Sparse object sobj.
  //! CCoorMatrix<double> ccoor(sobj);
  //! \endcode
  //!
  //! It is possible in any moment to change the sizes and the maximum
  //! number of nonzero by the `resize` methods:
  //!
  //! \code
  //! ccoor.resize(nrow, ncol, nnz);
  //! ccoor.resize(sp);
  //! ccoor.resize(ccoor1);
  //! ccoor.resize(sobj);
  //! \endcode
  //!
  //! - `real_type & insert(i,j)`
  //!    this method permit to insert an item in the matrix.  For example
  //!    matrix `A`  previously defined can be constructed with
  //!
  //! \code
  //! CCoorMatrix<double> A(6,7);
  //! A.insert(0, 0) = 1;    A.insert(0, 1) = 2; A.insert(0, 6) = 9;  A.insert(1, 1) = -1;
  //! A.insert(1, 2) = 0;    A.insert(2, 3) = 3; A.insert(2, 4) = 4;  A.insert(3, 1) = 2;
  //! A.insert(3, 4) = 5;    A.insert(4, 0) = 2; A.insert(4, 2) = -2; A.insert(4, 5) = 1;
  //! A.insert(5, 1) = -1.5; A.insert(5, 6) = 1; A.internal_order();
  //! \endcode
  //!
  //! - `internal_order()`
  //!   This methods reorder internally the nonzero elements of
  //!   `CCoorMatrix<double>` in such a way if
  //!   `k1 < k2` we have one of the two following cases
  //!
  //!   -# `I(k1) <  I(k2)`
  //!   -# `I(k1) == I(k2)` and `J(k1) <= J(k2)`
  //!
  //!   Moreover all duplicated entries are added togheter.
  //!
  //! - `bool is_ordered()`
  //!   this methods return `true` if the elements inside the class
  //!   `CCoorMatrix<double>` are ordered, `false` otherwise.
  //!
  template <typename T>
  class CCoorMatrix : public Sparse<T,CCoorMatrix<T> > {
    using MATRIX = CCoorMatrix<T>;
    using SPARSE = Sparse<T,MATRIX>;
  public:
    using real_type = T; //!< the type of the elements of the matrix

  private:

    Vector<integer>   I, J, C;
    Vector<real_type> A;

    mutable integer ipos;

    template <typename MAT, typename Compare> inline
    void
    convert( SparseBase<MAT> const & M, Compare cmp ) {
      // count nonzero
      integer nz{0};
      for ( M.Begin(); M.End(); M.Next() )
        if ( cmp( M.row(), M.column() ) ) ++nz;

      resize( M.nrows(), M.ncols(), nz );

      for ( M.Begin(); M.End(); M.Next() ) {
        integer i = M.row();
        integer j = M.column();
        if ( cmp(i,j) ) M.assign( insert(i,j) );
      }
      internal_order();
    }

  public:

    //!
    //! Initialize and empty sparse compressed coordinate matrix
    //! of `0` rows and `0` columns.
    //!
    CCoorMatrix(void) {}

    //!
    //! Initialize and empty empty sparse compressed coordinate matrix
    //! of `nr` rows and `nc` columns with reserved room for `mnnz` nonzeros.
    //!
    CCoorMatrix( integer nr, integer nc, integer mnnz = SPARSETOOL_DEFAULT_NNZ )
    { resize(nr, nc, mnnz); }

    //!
    //! Insert the element of the sparse compressed coordinate matrix `M`.
    //! which satify `cmp` to `*this`.
    //!
    //! \param M   sparse compressed coordinate matrix to be copied
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename Compare>
    CCoorMatrix( CCoorMatrix<T> const & M, Compare cmp )
    { convert(M,cmp); }

    //!
    //! Copy the sparse compressed coordinate matrix `M` to `*this`.
    //!
    //! \param M sparse compressed coordinate matrix to be copied
    //!
    CCoorMatrix(CCoorMatrix<T> const & M)
    { convert(M,all_ok()); }

    //!
    //! Insert the element of the sparse compressed coordinate matrix `M`.
    //! which satify `cmp` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    CCoorMatrix(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Copy the sparse compressed coordinate matrix `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    CCoorMatrix(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    //!
    //! Copy the sparse compressed coordinate matrix `M` to `*this`.
    //!
    //! \param M sparse compressed coordinate matrix
    //!
    CCoorMatrix<T> &
    operator = (CCoorMatrix<T> const & M) {
      if ( &M == this ) return *this; // avoid copy to itself
      resize( M.nrows(), M.ncols(), M.nnz() );
      SPARSE::sp_nnz        = M.nnz();
      SPARSE::sp_is_ordered = M.is_ordered();

      A.resize(SPARSE::sp_nnz); A = M.A.head(SPARSE::sp_nnz);
      I.resize(SPARSE::sp_nnz); I = M.I.head(SPARSE::sp_nnz);
      J.resize(SPARSE::sp_nnz); J = M.J.head(SPARSE::sp_nnz);
      C.resize(SPARSE::sp_nnz); C = M.C.head(SPARSE::sp_nnz);
      if ( !SPARSE::sp_is_ordered ) internal_order();
      return *this;
    }

    // convert => CCoorMatrix
    //!
    //! Copy the sparse matrix `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    CCoorMatrix<T> & operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return *this; }

    //!
    //! Initialize and empty sparse compressed coordinate
    //! matrix of `nr` rows and `nc` columns
    //! with reserved room for `mnnz` nonzeros.
    //!
    void
    resize( integer nr, integer nc, integer mnnz = SPARSETOOL_DEFAULT_NNZ ) {
      SPARSE::setup(nr, nc);
      I.resize(mnnz);
      J.resize(mnnz);
      A.resize(mnnz+1);
      A(0) = T(0);
      C.resize( SPARSE::sp_ncols + 1 );
    }

    //!
    //! Initialize the sparse compressed coordinate matrix and insert
    //! the elements M(i,j) which satify `cmp(i,j)` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    void
    resize( SparseBase<MAT> const & M, Compare cmp )
    { convert(M,cmp); }

    //!
    //! Initialize the sparse compressed coordinate matrix and insert the element of the matrix `M`.
    //! to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    void
    resize( SparseBase<MAT> const & M )
    { convert(M,all_ok()); }

    //!
    //! Check if pattern is empty.
    //!
    bool empty(void) const { return SPARSE::sp_nnz == 0; }

    //!
    //! Order the internal structure to permit random access to nonzeros.
    //!
    void
    internal_order() {
      QuickSortIJ2<integer,T>( I.data(), J.data(), A.data(), SPARSE::sp_nnz );
      // eliminate duplicate elements
      integer i1{0}, i{0};
      SPARSE::sp_lower_nnz = 0;
      SPARSE::sp_diag_nnz  = 0;
      SPARSE::sp_upper_nnz = 0;
      while ( i1 < SPARSE::sp_nnz ) {
        // copy i1 => i
        I(i) = I(i1); J(i) = J(i1); A(i) = A(i1);
        // setup statistic
        SPARSE::ldu_count(I(i),J(i));
        // find i1 != i
        for ( ++i1; i1 < SPARSE::sp_nnz && I(i1) == I(i) && J(i1) == J(i); ++i1 )
          A(i) += A(i1);
        ++i;
      }
      SPARSE::sp_nnz        = i;
      SPARSE::sp_is_ordered = true;
      I.conservativeResize(SPARSE::sp_nnz);
      J.conservativeResize(SPARSE::sp_nnz);
      A.conservativeResize(SPARSE::sp_nnz+1);
      A(SPARSE::sp_nnz) = T(0);
      C.resize(SPARSE::sp_ncols+1);

      integer nc{0};
      C(0) = 0;
      for ( integer k = 1; k < SPARSE::sp_nnz; ++k ) {
        UTILS_ASSERT(
          J(k-1) <= J(k),
          "CCoorMatrix::internal_order() internal error\n"
          "J({}) is not <= J({})\n",
          k-1, k
        );
        integer kk = J(k) - J(k-1);
        while ( kk-- > 0 ) C(++nc) = k;
      }
      UTILS_ASSERT(
        nc < SPARSE::sp_ncols || SPARSE::sp_nnz == 0,
        "CCoorMatrix::internal_order()\n"
        "nc = {} is not less of sp_ncols = {}, sp_nnz = {}\n",
        nc, SPARSE::sp_ncols, SPARSE::sp_nnz
      );
      while ( ++nc <= SPARSE::sp_ncols ) C(nc) = SPARSE::sp_nnz;
    }

    //!
    //! Return the position of the element `(i,j)`
    //! in the vector storing elements.
    //!
    integer
    position( integer i, integer j ) const {
      SPARSE::test_index(i,j);
      integer lo  = C(j);
      integer hi  = C(j+1);
      integer len = hi - lo;

      while ( len > 0 ) {
        integer half = len / 2;
        integer mid  = lo + half;
        if ( I(mid) < i ) {
          lo = mid + 1;
          len -= half + 1;
        } else len = half;
      }
      if ( lo == hi || I(lo) != i ) lo = SPARSE::sp_nnz;
      return lo;
    }

    //!
    //! Insert `(i,j)` in the sparse pattern.
    //!
    real_type &
    insert( integer i, integer j ) {
      SPARSE::test_index(i,j);
      if ( I.size() <= SPARSE::sp_nnz ) {
        integer newdim = SPARSE::sp_nnz + (SPARSE::sp_nnz>>3) + 100;
        I.conservativeResize( newdim );
        J.conservativeResize( newdim );
        A.conservativeResize( newdim+1 );
      }
      I(SPARSE::sp_nnz) = i;
      J(SPARSE::sp_nnz) = j;
      SPARSE::sp_nnz++;
      SPARSE::sp_is_ordered = false;
      return A(SPARSE::sp_nnz-1);
    }

    real_type const &
    operator [] (integer idx) const {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    real_type &
    operator [] (integer idx) {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    void setZero() { A = real_type(0); }
    void scale_values( real_type const & s ) { A *= s; }

    real_type const &
    value( integer i, integer j ) const
    { return A(position(i,j)); }

    real_type const &
    operator ()( integer i, integer j ) const {
      integer pos = position(i,j);
      UTILS_ASSERT(
        pos != SPARSE::sp_nnz,
        "CCoorMatrix({},{}) referring to a non existent element", i, j
      );
      return A(pos);
    }

    real_type &
    operator ()( integer i, integer j ) {
      integer pos = position(i,j);
      UTILS_ASSERT(
        pos != SPARSE::sp_nnz,
        "CCoorMatrix({},{}) referring to a non existent element", i, j
      );
      return A(pos);
    }

    bool
    exists( integer i, integer j ) {
      return position(i,j) != SPARSE::sp_nnz;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<real_type> const & getA(void) const { return A; } //!< return the value vector
    Vector<real_type>       & getA(void)       { return A; } //!< return the value vector
    Vector<integer>   const & getI(void) const { return I; } //!< return the row index vector
    Vector<integer>   const & getJ(void) const { return J; } //!< return the column index vector
    Vector<integer>   const & getC(void) const { return C; } //!< return the pointer of the column index (valid after ordering)

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // ITERATOR
    inline void Begin (void) const { ipos = 0; }
    inline void Next  (void) const { ++ipos; }
    inline bool End   (void) const { return ipos < SPARSE::sp_nnz; }

    integer           row    (void) const { return I(ipos); }
    integer           column (void) const { return J(ipos); }
    real_type const & value  (void) const { return A(ipos); }
    real_type       & value  (void)       { return A(ipos); }

    //!
    //! Assign the pointed element to `rhs`.
    //!
    template <typename TS>
    void assign( TS & rhs ) const { rhs = A(ipos); }

    //!
    //! perform the operation `res += s * (A * x)`
    //!
    template <typename VEC> inline
    void
    add_S_mul_M_mul_V(
      VEC       & res,
      T   const & s,
      VEC const & x
    ) const {
      UTILS_ASSERT0(
        std::addressof(res) != std::addressof(x),
        "Sparse_tool: [res += s * (A * x)] add_S_mul_M_mul_V, `res` and `x` cant be the same"
      );
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_nrows,
        "Sparse_tool: [res += s * (A * x)] result vector too small"
      );
      integer   const * pI = I.data();
      integer   const * pJ = J.data();
      real_type const * pA = A.data();
      SPARSELIB_LOOP( SPARSE::sp_nnz, res(*pI++) += s * *pA++ * x(*pJ++) );
    }

    //!
    //! perform the operation `res += s * (A^T * x)`
    //!
    template <typename VEC> inline
    void
    add_S_mul_Mt_mul_V(
      VEC       & res,
      T   const & s,
      VEC const & x
    ) const {
      UTILS_ASSERT0(
        std::addressof(res) != std::addressof(x),
        "Sparse_tool: [res += s * (A^T * x)] add_S_mul_Mt_mul_V, `res` and `x` cant be the same"
      );
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_ncols,
        "Sparse_tool: [res += s * (A^T * x)] result vector too small"
      );
      integer   const * pI = I.data();
      integer   const * pJ = J.data();
      real_type const * pA = A.data();
      SPARSELIB_LOOP( SPARSE::sp_nnz, res(*pJ++) += s * *pA++ * x(*pI++) );
    }

    //!
    //! perform the operation `res += s * (A * x)`
    //!
    template <typename VRES, typename VB> inline
    void
    add_S_mul_M_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_nrows,
        "Sparse_tool: [res += s * (A * x)] result vector too small"
      );
      integer   const * pI = I.data();
      integer   const * pJ = J.data();
      real_type const * pA = A.data();
      SPARSELIB_LOOP( SPARSE::sp_nnz, res(*pI++) += s * *pA++ * x(*pJ++) );
    }

    //!
    //! perform the operation `res += s * (A^T * x)`
    //!
    template <typename VRES, typename VB> inline
    void
    add_S_mul_Mt_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_ncols,
        "Sparse_tool: [res += s * (A^T * x)] result vector too small"
      );
      integer   const * pI = I.data();
      integer   const * pJ = J.data();
      real_type const * pA = A.data();
      SPARSELIB_LOOP( SPARSE::sp_nnz, res(*pJ++) += s * *pA++ * x(*pI++) );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  SPARSELIB_MUL_STRUCTURES(CCoorMatrix)
  #endif

  /*
  //  #     # #     # #       ####### ### ######  #       #     #
  //  ##   ## #     # #          #     #  #     # #        #   #
  //  # # # # #     # #          #     #  #     # #         # #
  //  #  #  # #     # #          #     #  ######  #          #
  //  #     # #     # #          #     #  #       #          #
  //  #     # #     # #          #     #  #       #          #
  //  #     #  #####  #######    #    ### #       #######    #
  */
  //!
  //! Perform matrix multiplication of a matrix in Compressed Coordinate
  //! with a full vector `x`.
  //!
  //! \code res += s * (A * x) \endcode
  //! \param nrows numer of rows
  //! \param RR    vector of row pointers
  //! \param JJ    vector of column indexes
  //! \param AA    vector of values
  //! \param s     scalar
  //! \param res   vector or the result
  //! \param x     vector used in multiplication
  //!
  template <typename T, typename VR, typename VX> inline
  void
  S_mul_M_mul_V(
    integer                 nrows,
    Vector<integer> const & RR,
    Vector<integer> const & JJ,
    Vector<T>       const & AA,
    T               const & s,
    VR                    & res,
    VX              const & x
  ) {
    UTILS_ASSERT0(
      std::addressof(res) != std::addressof(x),
      "Sparse_tool: [res += s * (A * x)] S_mul_M_mul_V, `res` and `x` cant be the same"
    );
    UTILS_ASSERT0(
      res.size() >= nrows,
      "Sparse_tool: [res += s * (A * x)] result vector too small"
    );
    integer const * R = RR.data();
    integer const * J = JJ.data();
    T       const * A = AA.data();
    for ( integer ir = 0; ir < nrows; ++ir, ++R ) {
      T bf(0);
      SPARSELIB_LOOP( R[1] - R[0], bf += *A++ * x(*J++) );
      res(ir) += s * bf;
    }
  }

  //!
  //! Perform matrix multiplication of a transopose of a matrix in Compressed Coordinate
  //! with a full vector `x`.
  //!
  //! \code res += s * (A^T * x) \endcode
  //! \param nrows numer of rows
  //! \param RR    vector of row pointers
  //! \param JJ    vector of column indexes
  //! \param AA    vector of values
  //! \param s     scalar
  //! \param res   vector or the result
  //! \param x     vector used in multiplication
  //!
  template <typename T, typename VR, typename VX> inline
  void
  S_mul_Mt_mul_V(
    integer                 nrows,
    Vector<integer> const & RR,
    Vector<integer> const & JJ,
    Vector<T>       const & AA,
    T               const & s,
    VR                    & res,
    VX              const & x
  ) {
    UTILS_ASSERT0(
      std::addressof(res) != std::addressof(x),
      "Sparse_tool: [res += s * (A^T * x)] S_mul_Mt_mul_V, `res` and `x` cant be the same"
    );
    UTILS_ASSERT0(
      res.size() >= nrows,
      "Sparse_tool: [res += s * (A^T * x)] result vector too small"
    );
    integer const * R = RR.data();
    integer const * J = JJ.data();
    T       const * A = AA.data();
    for ( integer ir = 0; ir < nrows; ++ir, ++R ) {
      T bf = s*x(ir);
      SPARSELIB_LOOP( R[1] - R[0], res(*J++) += *A++ * bf );
    }
  }

  /*
  //
  //  #####  ######
  // #     # #     #   ####   #    #
  // #       #     #  #    #  #    #
  // #       ######   #    #  #    #
  // #       #   #    #    #  # ## #
  // #     # #    #   #    #  ##  ##
  //  #####  #     #   ####   #    #
  //
  // #     #
  // ##   ##    ##    #####  #####   #  #    #
  // # # # #   #  #     #    #    #  #   #  #
  // #  #  #  #    #    #    #    #  #    ##
  // #     #  ######    #    #####   #    ##
  // #     #  #    #    #    #   #   #   #  #
  // #     #  #    #    #    #    #  #  #    #
  */

  //!
  //!
  //! The class `CRowMatrix<T>` implement a `Compressed Rows`
  //! storage sparse scheme.  It consists of two big vector of
  //! integer which contain the coordinate of nonzero elements and a big one
  //! of real number which contain the values.  We call `R` the vector
  //! that store the start position of each row, `J` the vector that
  //! store the column coordinate and `A`  the vector that store the
  //! nonzero values.  For example the sparse matrix
  //!
  //! \f[
  //!  \mathbf{A} = \begin{pmatrix}
  //!   1 & 2    &    &   &   &   & 9 \\
  //!     & -1   & 0  &   &   &   &   \\
  //!     &      &    & 3 & 4 &   &   \\
  //!     & 2    &    &   & 5 &   &   \\
  //!   2 &      & -2 &   &   & 1 &   \\
  //!     & -1.5 &    &   &   &   & -1
  //!  \end{pmatrix}
  //! \f]
  //! \f[
  //!   \begin{array}{c|cccccccccccccc}
  //!      R & 0 & 3 & 5 & 7 & 9 & 12 \\
  //!      \hline
  //!      J & 0 & 1 & 6 & 1 & 2 & 3 & 4 & 1 & 4 & 0 & 2 & 5 & 1 & 6 \\
  //!      \hline
  //!      A & 1 & 2 & 9 &-1 & 0 & 3 & 4 & 2 & 5 & 2 &-2 & 1 &-1.5 & -1
  //!    \end{array}
  //! \f]
  //!
  //! To define a `CCoorMatrix<T>` class you can use one of the following
  //! scripture
  //!
  //! \code
  //! CRowMatrix<double> crow;         // instancean empty CRowMatrix<double> class.
  //! CRowMatrix<double> crow(sp);     // instance a CRowMatrix<double> class with the sparsity pattern
  //!                                  // defined in the sp SparsePattern class.
  //! CRowMatrix<double> crow(crow1);  // instance a CRowMatrix<double> class which is the copy
  //!                                  // of the CRowMatrix<double> object crow1.
  //! CRowMatrix<double> crow(sobj);   // instance a CRowMatrix<double> class which is the copy
  //!                                  // of the Sparse object crow1.
  //! \endcode
  //!
  //! It is possible in any moment to change the sizes and the maximum
  //! number of nonzero by the `resize` method:
  //!
  //! \code
  //! crow.resize(sp);
  //! crow.resize(crow1);
  //! crow.resize(sobj);
  //! \endcode
  //!
  template <typename T>
  class CRowMatrix : public Sparse<T,CRowMatrix<T> > {
    using MATRIX = CRowMatrix<T>;
    using SPARSE = Sparse<T,MATRIX>;
  public:
    using real_type = T; //!< the type of the elements of the matrix

  private:

    Vector<real_type> A;
    Vector<integer>   R;
    Vector<integer>   J;

    mutable integer iter_row;
    mutable integer iter_ptr;

    void
    internal_order() {
      UTILS_ASSERT0(
        R[SPARSE::sp_nrows] == SPARSE::sp_nnz,
        "Sparse_tool: CRowMatrix::internal_order() bad data for matrix\n"
      );
      integer ii, kk, rk, rk1;
      SPARSE::sp_lower_nnz = 0;
      SPARSE::sp_diag_nnz  = 0;
      SPARSE::sp_upper_nnz = 0;
      for ( ii = 0, rk = R(0); ii < SPARSE::sp_nrows; ++ii, rk = rk1 ) {
        rk1 = R(ii+1);
        if ( rk1 > rk ) { // skip empty rows
          QuickSortI<integer,T>( &J(rk), &A(rk), rk1 - rk );
          // setup statistic
          for ( kk = rk; kk < rk1; ++kk ) SPARSE::ldu_count(ii,J(kk));
  #ifdef SPARSETOOL_DEBUG
          for ( kk = rk+1; kk < rk1; ++kk )
            UTILS_ASSERT0(
              J(kk-1) < J(kk),
              "Sparse_tool: CRowMatrix::internal_order() failed\n"
            );
  #endif
        }
      }
      SPARSE::sp_nnz = R(SPARSE::sp_nrows);
    }

    template <typename MAT, typename Compare>
    void
    convert(SparseBase<MAT> const & M, Compare cmp) {

      SPARSE::setup( M.nrows(), M.ncols() );

      // step 0: Count nonzero
      SPARSE::sp_nnz = 0;
      for ( M.Begin(); M.End();  M.Next() )
        if ( cmp(M.row(), M.column() ) ) ++SPARSE::sp_nnz;

      A.resize( SPARSE::sp_nnz + 1 );
      J.resize( SPARSE::sp_nnz );
      R.resize( SPARSE::sp_nrows + 1 );
      A(SPARSE::sp_nnz) = 0;
      R.setZero();

      // step 1: Evaluate not zero pattern
      for ( M.Begin(); M.End();  M.Next() ) {
        integer i = M.row();
        integer j = M.column();
        if ( cmp(i, j) ) ++R(i);
      }
      for ( integer k = 0; k < SPARSE::sp_nrows; ++k ) R(k+1) += R(k);

      // step 2: Fill matrix
      for ( M.Begin(); M.End(); M.Next() ) {
        integer i = M.row();
        integer j = M.column();
        if ( cmp(i,j) ) {
          integer ii = --R(i);
          M.assign(A(ii));
          J(ii) = j;
        }
      }

      UTILS_ASSERT0(
        R(0) == 0 && R(SPARSE::sp_nrows) == SPARSE::sp_nnz,
        "CRowMatrix::resize(Sparse & A) failed\n"
      );
      // step 3: internal_order matrix
      internal_order();
    }

  public:

    //!
    //! Initialize and empty sparse compressed row matrix
    //! of `0` rows and `0` columns.
    //!
    CRowMatrix(void) {}

    //!
    //! Insert the element of the sparse compressed row matrix `M`.
    //! which satify `cmp` to `*this`.
    //!
    //! \param M   sparse compressed row matrix to be copied
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename Compare>
    CRowMatrix(CRowMatrix<T> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Copy the sparse compressed row matrix `M` to `*this`.
    //!
    //! \param M sparse compressed row matrix to be copied
    //!
    CRowMatrix(CRowMatrix<T> const & M)
    { convert(M,all_ok()); }

    //!
    //! Insert the element of the sparse compressed row matrix `M`.
    //! which satify `cmp` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    CRowMatrix(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Copy the sparse compressed row matrix `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    CRowMatrix(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    //!
    //! Copy the sparse compressed row matrix `M` to `*this`.
    //!
    //! \param M sparse compressed row matrix
    //!
    CRowMatrix<T> &
    operator = (CRowMatrix<T> const & M) {
      if ( &M == this ) return *this; // avoid copy to itself
      SPARSE::setup( M.nrows(), M.ncols() );
      SPARSE::sp_nnz        = M.nnz();
      SPARSE::sp_is_ordered = M.is_ordered();
      A.resize(SPARSE::sp_nnz);     A = M.A.head(SPARSE::sp_nnz);
      R.resize(SPARSE::sp_nrows+1); R = M.R.head(SPARSE::sp_nrows+1);
      J.resize(SPARSE::sp_nnz);     J = M.J.head(SPARSE::sp_nnz);
      if ( !SPARSE::sp_is_ordered ) internal_order();
      return *this;
    }

    // convert => CRowMatrix
    //!
    //! Copy the sparse matrix `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    CRowMatrix<T> & operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return *this; }

    //!
    //! Initialize the sparse compressed row matrix and insert
    //! the elements M(i,j) which satify `cmp(i,j)` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    void
    resize(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Initialize the sparse compressed row matrix and insert the element of the matrix `M`.
    //! to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    void
    resize(SparseBase<MAT> const & M)
    { resize(M,all_ok()); }

    void
    scaleRow( integer nr, real_type const & val ) {
      SPARSE::test_row(nr);
      real_type * pA = A.data() + R(nr);
      real_type * pB = A.data() + R(nr+1);
      while ( pA < pB ) *pA++ *= val;
    }

    void
    scaleColumn(integer nc, real_type const & val) {
      SPARSE::test_col(nc);
      for ( integer i = 0; i < SPARSE::sp_nrows; ++i ) {
        integer pos = position(i,nc,true);
        if ( pos != SPARSE::sp_nnz ) A(pos) *= val;
      }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<real_type> const & getA(void) const { return A; } //!< return the value vector
    Vector<real_type>       & getA(void)       { return A; } //!< return the value vector
    Vector<integer>   const & getR(void) const { return R; } //!< return the row pointer
    Vector<integer>   const & getJ(void) const { return J; } //!< return the column index vector

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    integer
    position( integer i, integer j ) const {
      SPARSE::test_index(i,j);
      integer lo  = R(i);
      integer hi  = R(i+1);
      integer len = hi - lo;
      while (len > 0) {
        integer half = len / 2;
        integer mid = lo + half;
        if ( J(mid) < j ) {
          lo = mid + 1;
          len -= half + 1;
        } else len = half;
      }
      if ( lo == hi || J(lo) != j ) lo = SPARSE::sp_nnz;
      return lo;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    void setZero() { A = real_type(0); }
    void scale_values( real_type const & s ) { A *= s; }

    real_type const &
    operator ()( integer i, integer j ) const {
      integer pos = position(i,j);
      UTILS_ASSERT(
        pos != SPARSE::sp_nnz,
        "CRowMatrix({},{}) referring to a non existent element\n", i, j
      );
      return A(pos);
    }

    real_type &
    operator ()( integer i, integer j ) {
      integer pos = position(i,j);
      UTILS_ASSERT(
        pos != SPARSE::sp_nnz,
        "CRowMatrix({},{}) referring to a non existent element\n", i, j
      );
      return A(pos);
    }

    real_type const &
    value( integer i, integer j) const
    { return A(position(i,j)); }

    bool
    exists( integer i, integer j ) {
      return position(i,j) != SPARSE::sp_nnz;
    }

    real_type const &
    operator [] (integer idx) const
    { SPARSE::test_nnz(idx); return A(idx); }

    real_type &
    operator [] (integer idx)
    { SPARSE::test_nnz(idx); return A(idx); }

    // ****************************************************************
    // ITERATOR
    void Begin(void) const {
      iter_row = iter_ptr = 0;
      while ( R(iter_row) == 0 ) ++iter_row;
    }
    void Next (void) const {
      ++iter_ptr;
      while ( iter_row <= SPARSE::sp_nrows && iter_ptr >= R(iter_row) ) ++iter_row;
    }
    bool End(void) const { return iter_ptr < SPARSE::sp_nnz; }

    integer           row    (void) const { return iter_row-1; }
    integer           column (void) const { return J(iter_ptr); }
    real_type const & value  (void) const { return A(iter_ptr); }
    real_type       & value  (void)       { return A(iter_ptr); }

    //!
    //! Assign the pointed element to `rhs`.
    //!
    template <typename TS>
    void assign( TS & rhs ) const { rhs = A(iter_ptr); }

    //!
    //! Perform the operation `res += s * (A * x)`
    //!
    template <typename VRES, typename VB> inline
    void
    add_S_mul_M_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      S_mul_M_mul_V( SPARSE::sp_nrows, R, J, A, s, res, x );
    }

    //!
    //! Perform the operation `res += s * (A ^ x)`
    //!
    template <typename VRES, typename VB> inline
    void
    add_S_mul_Mt_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      S_mul_Mt_mul_V( SPARSE::sp_nrows, R, J, A, s, res, x );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  SPARSELIB_MUL_STRUCTURES(CRowMatrix)
  #endif

  /*
  //  #####   #####
  // #     # #     #   ####   #
  // #       #        #    #  #
  // #       #        #    #  #
  // #       #        #    #  #
  // #     # #     #  #    #  #
  //  #####   #####    ####   ######
  //
  // #     #
  // ##   ##    ##    #####  #####   #  #    #
  // # # # #   #  #     #    #    #  #   #  #
  // #  #  #  #    #    #    #    #  #    ##
  // #     #  ######    #    #####   #    ##
  // #     #  #    #    #    #   #   #   #  #
  // #     #  #    #    #    #    #  #  #    #
  */

  //!
  //! Compressed Column Matrix Storage
  //! --------------------------------
  //!
  //! The class `CColMatrix<T>` implement a `Compressed Columns`
  //! storage sparse scheme.  It consists of two big vector of
  //! integer which contain the coordinate of nonzero elements and a big one
  //! of real number which contain the values.  We call `I` the vector
  //! that store the row index of each elements, `C` the vector that
  //! store the starting position of the columns and `A`  the vector
  //! that store the nonzero values.  For example the sparse matrix
  //!
  //! \f[
  //! \mathbf{A} = \begin{pmatrix}
  //!   1 & 2    &    &   &   &   & 9 \\
  //!     & -1   & 0  &   &   &   &   \\
  //!     &      &    & 3 & 4 &   &   \\
  //!     & 2    &    &   & 5 &   &   \\
  //!   2 &      & -2 &   &   & 1 &   \\
  //!     & -1.5 &    &   &   &   & -1
  //! \end{pmatrix}
  //! \f]
  //! \f[
  //! \begin{array}{c|cccccccccccccc}
  //!   I & 0 & 4 & 0 & 1 & 3 & 5 & 1 & 4 & 2 & 2 & 3 & 4 & 0 & 5 \\
  //!   \hline
  //!   C & 0 & 2 & 6 & 8 & 9 & 11 & 12 \\
  //!   \hline
  //!   A & 1 & 2 & 2 & -1 & 2 & -1.5 & 0 & -2 & 3 & 4 & 5 & 1 & 9 & -1
  //! \end{array}
  //! \f]
  //! Notice that the index follows the C convention, starting from **0**.
  //! The class `CColMatrix<T>` try to manage such a structure in a
  //! simple way for the user.  To define a `CColMatrix<T>` class you
  //! can use one of the following scripture
  //!
  //! \code
  //! CColMatrix<double> ccol;        // instance an empty CColMatrix<double> class.
  //! CColMatrix<double> ccol(sp);    // instance a CColMatrix<double> class with the sparsity
  //!                                 // pattern defined in the SparsePattern class sp.
  //! CColMatrix<double> ccol(ccol1); // instance a CColMatrix<double> class which is the copy
  //!                                 // of the CColMatrix<double> object ccol1.
  //! CColMatrix<double> ccol(sobj);  // instance a CColMatrix<double> class which is the copy
  //!                                 // of the Sparse object sobj.
  //! \endcode
  //!
  //! It is possible in any moment to change the sizes and the maximum
  //! number of nonzero by the `resize` method:
  //!
  //! \code
  //!   ccol.resize(sp);
  //!   ccol.resize(ccol1);
  //!   ccol.resize(sobj);
  //! \endcode
  //!
  //!
  template <typename T>
  class CColMatrix : public Sparse<T,CColMatrix<T> > {
    using MATRIX = CColMatrix<T>;
    using SPARSE = Sparse<T,MATRIX>;
  public:
    using real_type = T; //!< the type of the elements of the matrix

  private:
    Vector<real_type> A;
    Vector<integer>   C, I;

    mutable integer iter_col;
    mutable integer iter_ptr;

    void
    internal_order() {
      integer jj, kk, ck, ck1;
      SPARSE::sp_lower_nnz = 0;
      SPARSE::sp_diag_nnz  = 0;
      SPARSE::sp_upper_nnz = 0;
      for ( jj = 0, ck = C(0); jj < SPARSE::sp_ncols; ++jj, ck = ck1 ) {
        ck1 = C(jj+1);
        if ( ck1 > ck ) { // skip empty columns
          QuickSortI<integer,T>( &I(ck), &A(ck), ck1 - ck );
          // setup statistic
          for ( kk = ck; kk < ck1; ++kk ) SPARSE::ldu_count(I(kk),jj);
  #ifdef SPARSETOOL_DEBUG
          for ( kk = ck+1; kk < ck1; ++kk )
            UTILS_ASSERT0(
              I(kk-1) < I(kk),
              "Sparse_tool: CColMatrix::internal_order() failed\n"
            );
  #endif
        }
      }
      SPARSE::sp_nnz = C(SPARSE::sp_ncols);
    }

    template <typename MAT, typename Compare>
    void
    convert(SparseBase<MAT> const & M, Compare cmp) {
      SPARSE::setup( M.nrows(), M.ncols() );

      // step 0: Count nonzero
      SPARSE::sp_nnz = 0;
      for ( M.Begin(); M.End();  M.Next() )
        if ( cmp(M.row(), M.column() ) ) ++SPARSE::sp_nnz;

      A.resize( SPARSE::sp_nnz + 1 );
      I.resize( SPARSE::sp_nnz );
      C.resize( SPARSE::sp_ncols + 1 );
      C.setZero();
      A(SPARSE::sp_nnz) = 0;

      // step 1: Evaluate not zero
      for ( M.Begin(); M.End(); M.Next() ) {
        integer i = M.row();
        integer j = M.column();
        if ( cmp(i,j) ) ++C(j);
      }
      for ( integer k = 0; k < SPARSE::sp_ncols; ++k ) C(k+1) += C(k);

      // step 2: Fill matrix
      for ( M.Begin(); M.End(); M.Next() ) {
        integer i = M.row();
        integer j = M.column();
        if ( cmp(i,j) ) {
          integer jj = --C(j);
          M.assign(A(jj));
          I(jj) = i;
        }
      }
      UTILS_ASSERT0(
        C(0) == 0 && C(SPARSE::sp_ncols) == SPARSE::sp_nnz,
        "CColMatrix::resize(Sparse & A) failed\n"
      );
      // step 3: internal_order matrix
      internal_order();
    }

  public:

    //!
    //! Initialize and empty sparse compressed column matrix
    //! of `0` rows and `0` columns.
    //!
    CColMatrix(void) {}

    //!
    //! Insert the element of the sparse compressed column matrix `M`.
    //! which satify `cmp` to `*this`.
    //!
    //! \param M   sparse compressed row matrix to be copied
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename Compare>
    CColMatrix(CColMatrix<T> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Copy the sparse compressed row matrix `M` to `*this`.
    //!
    //! \param M sparse compressed row matrix to be copied
    //!
    CColMatrix(CColMatrix<T> const & M)
    { convert(M,all_ok()); }

    //!
    //! Insert the element of the sparse compressed column matrix `M`.
    //! which satify `cmp` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    CColMatrix(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Copy the sparse compressed column matrix `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    CColMatrix(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    //!
    //! Copy the sparse compressed column matrix `M` to `*this`.
    //!
    //! \param M sparse compressed column matrix
    //!
    CColMatrix<T> &
    operator = (CColMatrix<T> const & M) {
      if ( &M == this ) return *this; // avoid copy to itself
      resize( M.nrows(), M.ncols(), M.nnz() );
      SPARSE::sp_nnz        = M.nnz();
      SPARSE::sp_is_ordered = M.is_ordered();
      A.resize( SPARSE::sp_nnz );     A = M.A.head( SPARSE::sp_nnz );
      C.resize( SPARSE::sp_ncols+1 ); C = M.C.head( SPARSE::sp_ncols+1 );
      I.resize( SPARSE::sp_nnz );     I = M.I.head( SPARSE::sp_nnz );
      if ( !SPARSE::sp_is_ordered ) internal_order();
      return *this;
    }

    // convert => CColMatrix
    //!
    //! Copy the sparse matrix `M` to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    CColMatrix<T> &
    operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return * this; }

    //!
    //! Initialize the sparse compressed column matrix and insert
    //! the elements M(i,j) which satify `cmp(i,j)` to `*this`.
    //!
    //! \param M   object derived from `Sparse`.
    //! \param cmp comparator struct for the selection of the element in pattern `sp`
    //!
    template <typename MAT, typename Compare>
    void
    resize(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    //!
    //! Initialize the sparse compressed column matrix and insert the element of the matrix `M`.
    //! to `*this`.
    //!
    //! \param M object derived from `Sparse`.
    //!
    template <typename MAT>
    void
    resize( SparseBase<MAT> const & M )
    { resize(M,all_ok()); }

    void
    scaleRow( integer nr, real_type const & val ) {
      SPARSE::test_row(nr);
      for ( integer j = 0; j < SPARSE::sp_ncols; ++j ) {
        integer pos = position(nr,j,true);
        if ( pos != SPARSE::sp_nnz ) A(pos) *= val;
      }
    }

    void
    scaleColumn( integer nc, real_type const & val ) {
      SPARSE::test_col(nc);
      real_type * pA = A.data() + C(nc);
      real_type * pB = A.data() + C(nc+1);
      while ( pA < pB ) *pA++ *= val;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<real_type> const & getA(void) const { return A; } //!< return the value vector
    Vector<real_type>       & getA(void)       { return A; } //!< return the value vector
    Vector<integer>   const & getI(void) const { return I; } //!< return the row index vector
    Vector<integer>   const & getC(void) const { return C; } //!< return the column pointer vector

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    integer
    position( integer i, integer j ) const {
      SPARSE::test_index(i,j);
      integer lo  = C(j);
      integer hi  = C(j+1);
      integer len = hi - lo;
      while (len > 0) {
        integer half = len / 2;
        integer mid = lo + half;
        if ( I(mid) < i ) {
          lo = mid + 1;
          len -= half + 1;
        } else len = half;
      }
      if ( lo == hi || I(lo) != i ) lo = SPARSE::sp_nnz;
      return lo;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    void setZero() { A = real_type(0); }
    void scale_values( real_type const & s ) { A *= s; }

    real_type const &
    operator () (integer i, integer j) const {
      integer pos = position(i,j);
      UTILS_ASSERT(
        pos != SPARSE::sp_nnz,
        "CColMatrix({},{}) referring to a non existent element\n", i, j
      );
      return A(pos);
    }

    real_type &
    operator () (integer i, integer j) {
      integer pos = position(i,j);
      UTILS_ASSERT(
        pos != SPARSE::sp_nnz,
        "CColMatrix({},{}) referring to a non existent element\n", i, j
      );
      return A(pos);
    }

    real_type const &
    value(integer i, integer j) const
    { return A(position(i,j)); }

    bool
    exists( integer i, integer j )
    { return position(i,j) != SPARSE::sp_nnz; }

    real_type const &
    operator [] (integer idx) const
    { SPARSE::test_nnz(idx); return A(idx); }

    real_type &
    operator [] (integer idx)
    { SPARSE::test_nnz(idx); return A(idx); }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    void Begin(void) const {
      iter_col = iter_ptr = 0;
      while ( C(iter_col) == 0 ) ++iter_col;
    }
    void Next (void) const {
      ++iter_ptr;
      while ( iter_ptr >= C(iter_col) && iter_col <= SPARSE::sp_ncols ) ++iter_col;
    }
    bool End(void) const { return iter_ptr < SPARSE::sp_nnz; }

    integer           row    (void) const { return I(iter_ptr); }
    integer           column (void) const { return iter_col-1;  }
    real_type const & value  (void) const { return A(iter_ptr); }
    real_type       & value  (void)       { return A(iter_ptr); }

    //!
    //! Assign the pointed element to `rhs`.
    //!
    template <typename TS>
    void assign( TS & rhs ) const { rhs = A(iter_ptr); }

    //!
    //! Perform the operation `res += s * (A * x)`.
    //!
    template <typename VRES, typename VB> inline
    void
    add_S_mul_M_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      S_mul_Mt_mul_V( SPARSE::sp_ncols, C, I, A, s, res, x );
    }

    //!
    //! Perform the operation `res += s * (A ^ x)`.
    //!
    template <typename VRES, typename VB> inline
    void
    add_S_mul_Mt_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      S_mul_M_mul_V( SPARSE::sp_ncols, C, I, A, s, res, x );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  SPARSELIB_MUL_STRUCTURES(CColMatrix)
  #endif

  /*
  // #######
  //    #     #####   #  #####
  //    #     #    #  #  #    #
  //    #     #    #  #  #    #
  //    #     #####   #  #    #
  //    #     #   #   #  #    #
  //    #     #    #  #  #####
  //
  // #     #
  // ##   ##    ##    #####  #####   #  #    #
  // # # # #   #  #     #    #    #  #   #  #
  // #  #  #  #    #    #    #    #  #    ##
  // #     #  ######    #    #####   #    ##
  // #     #  #    #    #    #   #   #   #  #
  // #     #  #    #    #    #    #  #  #    #
  */

  //!
  //! Tridiagonal Matrix
  //! ------------------
  //!
  //! The class `TridMatrix<T>` implement a sparse band matrix.
  //! It consists of a big matrix of `real_type` which contain the values
  //! of nonzero elements.  We call `M` this big matrix which
  //! represents an \f$ n\times m \f$ matrix which stores the rows of nonzero.
  //! For example the following band matrix
  //!
  //! \f[
  //!   \mathbf{A} = \begin{pmatrix}
  //!     2 & -1 &    &    &    &    &    \\
  //!     1 & 2  & -2 &    &    &    &    \\
  //!       & 2  & 2  & -3 &    &    &    \\
  //!       &    & 3  & 2  & -4 &    &    \\
  //!       &    &    & 4  & 2  & -5 &    \\
  //!       &    &    &    & 5  & 2  & -6 \\
  //!       &    &    &    &    & 6  & 2
  //!    \end{pmatrix}
  //!    \qquad
  //!    \begin{array}{c|c|c}
  //!       L & D & U \\
  //!      \hline
  //!      * & 2 & -1 \\
  //!      1 & 2 & -2 \\
  //!      2 & 2 & -3 \\
  //!      3 & 2 & -4 \\
  //!      4 & 2 & -5 \\
  //!      5 & 2 & -6 \\
  //!      6 & 2 & *  \\
  //!      \hline
  //!    \end{array}
  //! \f]
  //!
  //! where `*` means unused elements.  Notice that the index follows
  //! the C convention, starting from `0`.
  //! The class `TridMatrix<T>` manage such a structure in a simple
  //! way for the user.  To define a `TridMatrix<T>` class you can
  //! use one of the following scripture
  //!
  //! \code
  //! TridMatrix<double> tm;      // instance an empty TridMatrix<double> class
  //! TridMatrix<double> tm(100); // instance a 100 x 100 TridMatrix<double> class
  //! TridMatrix<double> tm(tm1); // instance a TridMatrix<double> class which is the copy
  //!                             // of the TridMatrix<double> class tm1.
  //! \endcode
  //!
  //! It is possible in any moment to change the sizes and the maximum
  //! number of nonzero by the `resize` method:
  //! \code
  //! tm.resize(n);
  //! \endcode
  //!
  //!
  template <typename T>
  class TridMatrix : public Sparse<T,TridMatrix<T> > {
    using MATRIX = TridMatrix<T>;
    using SPARSE = Sparse<T,MATRIX>;
    using VBASE  = typename Vector<T>::V_base;
  public:
    using real_type = T; //!< type of the element of the matrix

  private:

    Vector<T>         A;
    Eigen::Map<VBASE> L, D, U;

    mutable integer iter_counter, iter_row, iter_col;
    mutable Vector<T> TMP;

    template <typename VECTOR>
    void
    solve( VECTOR const & b, VECTOR & x ) const {

      VECTOR ll(SPARSE::sp_nrows-1), dd(SPARSE::sp_nrows);

      // LU DEC
      integer k{0};
      dd(0) = D(0);
      for ( k = 1; k < SPARSE::sp_nrows; ++k ) {
        ll(k-1) = L(k-1) / dd(k-1);
        dd(k)   = D(k) - ll(k-1) * U(k-1);
      }

      // SOLVE
      k = 0;
      x(0) = b(0);
      while ( ++k < SPARSE::sp_nrows ) x(k) = b(k) - ll(k-1) * x(k-1);
      --k;
      x(k) /= dd(k);
      while ( k > 0 ) {
        --k;
        x(k) = (x(k) - U(k) * x(k+1)) / dd(k);
      };

    }

  public:

    //!
    //! Contruct and exmpty `0x0` tridiagonal matrix.
    //!
    TridMatrix()
    : L(NULL,0)
    , D(NULL,0)
    , U(NULL,0)
    {}

    //!
    //! Contruct a `n x n` tridiagonal matrix.
    //!
    TridMatrix( integer ns ) : TridMatrix()
    { resize(ns); }

    //!
    //! Contruct a copy of the tridiagonal matrix `TM`.
    //!
    TridMatrix( TridMatrix<T> const & TM ) : TridMatrix()
    { resize(TM); }

    //!
    //! Contruct a copy of the tridiagonal matrix `TM`.
    //!
    TridMatrix<T> const &
    operator = ( TridMatrix<T> const & TM ) {
      if ( &TM == this ) return *this; // avoid copy to itself
      resize( A.nrows() );
      A = TM.A;
      return *this;
    }

    //!
    //! Resize to a `n x n` tridiagonal matrix.
    //!
    void
    resize( integer ns ) {
      A.resize(3*ns-2);
      integer kk{0};
      new (&L) Eigen::Map<VBASE>(A.data()+kk, ns-1 ); kk += ns-1;
      new (&D) Eigen::Map<VBASE>(A.data()+kk, ns   ); kk += ns;
      new (&U) Eigen::Map<VBASE>(A.data()+kk, ns-1 );
      TMP.resize( ns );
      SPARSE::setup(ns,ns);
      SPARSE::sp_nnz = 3*ns-2;
    }

    real_type const &
    operator [] ( integer idx ) const {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    real_type &
    operator [] ( integer idx ) {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    void
    scaleRow( integer nr, real_type const & val ) {
      SPARSE::test_row(nr);
      D(nr) *= val;
      if ( nr < SPARSE::sp_now ) U(nr)   *= val;
      if ( nr > 0              ) L(nr-1) *= val;
    }

    void
    scaleColumn( integer nc, real_type const & val ) {
      SPARSE::test_col(nc);
      D(nc) *= val;
      if ( nc > 0               ) U(nc-1) *= val;
      if ( nc < SPARSE::sp_nrow ) L(nc)   *= val;
    }

    //!
    //! \name Diagonal internal operation.
    //!
    //@{
    //!
    //! Set `s` to all the components of the diagonal,
    //! and `0` the others elements.  For example
    //!
    //! \f[
    //!   \mathbf{A} = \begin{pmatrix}
    //!       2 & -1 &    &    \\
    //!       1 & 2  & -1 &    \\
    //!         & 1  & 2  & -1 \\
    //!         &    & 1  & 2  \\
    //!   \end{pmatrix},\qquad
    //!   (\mathbf{A}=2.1) = \begin{pmatrix}
    //!       2.1 & 0    &      &     \\
    //!       0   & 2.1  & 0    &     \\
    //!           & 0    & 2.1  & 0   \\
    //!           &      & 0    & 2.1 \\
    //!    \end{pmatrix}\right]
    //! \f]
    //!
    TridMatrix<T> &
    operator = ( real_type const & s )
    { A.setZero(); D = s; return *this; }

    //!
    //! Add `s` to all the components of the diagonal,
    //! and `0` the others elements.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! (\mathbf{A}+=2.1) = \begin{pmatrix}
    //!     4.1 & -1   &      &     \\
    //!     1   & 4.1  & -1   &     \\
    //!         & 1    & 4.1  & -1  \\
    //!         &      & 1    & 4.1 \\
    //!  \end{pmatrix}
    //! \f]
    //!
    TridMatrix<T> &
    operator += ( real_type const & s )
    { D.array() += s; return *this; }

    //!
    //! Subtract `s` to all the components of the diagonal,
    //! and `0` the others elements.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! (\mathbf{A}-=2.1) = \begin{pmatrix}
    //!     -0.1 & -1   &      &     \\
    //!     1   & -0.1  & -1   &     \\
    //!         & 1    & -0.1  & -1  \\
    //!         &      & 1    & -0.1 \\
    //! \end{pmatrix}
    //! \f]
    //!
    TridMatrix<T> &
    operator -= ( real_type const & s )
    { D.array() -= s; return *this; }

    //!
    //! Multiply by `s` to all the nonzeros.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! (\mathbf{A}*=2) = \begin{pmatrix}
    //!     4 & -2 &    &    \\
    //!     2 & 4  & -2 &    \\
    //!       & 2  & 4  & -2 \\
    //!       &    & 2  & 4  \\
    //! \end{pmatrix}
    //! \f]
    //!
    TridMatrix<T> &
    operator *= ( real_type const & s )
    { A *= s; return *this; }

    //!
    //! Divide by `s` to all the nonzeros.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! (\mathbf{A}/=2) = \begin{pmatrix}
    //!     1   & -0.5 &      &      \\
    //!     0.5 & 1    & -0.5 &      \\
    //!         & 0.5  & 1    & -0.5 \\
    //!         &      & 0.5  & 1    \\
    //! \end{pmatrix}
    //! \f]
    //!
    TridMatrix<T> &
    operator /= ( real_type const & s )
    { A /= s; return *this; }

    //!
    //! Set the diagonal values of the sparse
    //! matrix to `v`.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! \mathbf{v} = \begin{pmatrix}
    //! 1 \\ 2\\ 3
    //! \end{pmatrix},
    //! \qquad
    //! (\mathbf{A}/=2) = \begin{pmatrix}
    //!     1 & 0  &   &   \\
    //!     0 & 2  & 0 &   \\
    //!       & 0  & 3 & 0 \\
    //!       &    & 0 & 0 \\
    //! \end{pmatrix}
    //! \f]
    //!
    template <typename TV>
    inline
    TridMatrix<T> &
    operator = ( Vector<TV> const & v )
    { A.setZero(); D = v; return *this; }

    //!
    //! Add  `v` to the diagonal values of the sparse
    //! matrix.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! \mathbf{v} = \begin{pmatrix}
    //! 1 \\ 2\\ 3
    //! \end{pmatrix},
    //! \qquad
    //! (\mathbf{A}\verb|+=|\mathbf{v}) = \begin{pmatrix}
    //!     3 & -1 &    &    \\
    //!     1 & 4  & -1 &    \\
    //!       & 1  & 5  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix}
    //! \f]
    //! \htmlonly </TD></TR></TABLE> \endhtmlonly
    //!
    template <typename TV>
    inline
    TridMatrix<T> &
    operator += ( Vector<TV> const & v )
    { D += v; return *this; }

    //!
    //! Subtract  `v`  to the diagonal values of the sparse
    //! matrix.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! \mathbf{v} = \begin{pmatrix}
    //! 1 \\ 2\\ 3
    //! \end{pmatrix},
    //! \qquad
    //! (\mathbf{A}\verb|-=|\mathbf{v}) = \begin{pmatrix}
    //!     3 & -1 &    &    \\
    //!     1 & 4  & -1 &    \\
    //!       & 1  & 5  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix}
    //! \f]
    //!
    template <typename TV>
    inline
    TridMatrix<T> &
    operator -= ( Vector<TV> const & v )
    { D -= v; return *this; }

    //!
    //! Set the diagonal values of the sparse
    //! matrix to the vector expression `e`.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! \mathbf{a} = \begin{pmatrix}
    //! 1 \\ 2\\ 3
    //! \end{pmatrix},\qquad
    //! \mathbf{b} = \begin{pmatrix}
    //! 1 \\ -1 \\ 2 \\ 1
    //! \end{pmatrix},
    //! \qquad
    //! (\mathbf{A}=\mathbf{a}+2*\mathbf{b}) = \begin{pmatrix}
    //!     3 & 0  &   &   \\
    //!     0 & 0  & 0 &   \\
    //!       & 0  & 7 & 0 \\
    //!       &    & 0 & 0 \\
    //! \end{pmatrix}
    //! \f]
    //!
    template <typename VEC_EXPR>
    inline
    TridMatrix<T> &
    operator = ( VEC_EXPR const & e )
    { A.setZero(); D = e; return *this; }

    //!
    //! Add to the diagonal values of the sparse
    //! matrix the vector expression `e`.  For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! \mathbf{a} = \begin{pmatrix}
    //! 1 \\ 2\\ 3
    //! \end{pmatrix},\qquad
    //! \mathbf{b} = \begin{pmatrix}
    //! 1 \\ -1 \\ 2 \\ 1
    //! \end{pmatrix},
    //! \qquad
    //! (\mathbf{A}\verb|+=|\mathbf{a}+2*\mathbf{b}) = \begin{pmatrix}
    //!     5 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 9  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix}
    //! \f]
    //!
    template <typename VEC_EXPR>
    inline
    TridMatrix<T> &
    operator += ( VEC_EXPR const & e )
    { TMP = e; D += TMP; return *this; }

    //!
    //! Add to the diagonal values of the sparse
    //! matrix the vector expression `e`. For example
    //!
    //! \f[
    //! \mathbf{A} = \begin{pmatrix}
    //!     2 & -1 &    &    \\
    //!     1 & 2  & -1 &    \\
    //!       & 1  & 2  & -1 \\
    //!       &    & 1  & 2  \\
    //! \end{pmatrix},\qquad
    //! \mathbf{a} = \begin{pmatrix}
    //! 1 \\ 2\\ 3
    //! \end{pmatrix},\qquad
    //! \mathbf{b} = \begin{pmatrix}
    //! 1 \\ -1 \\ 2 \\ 1
    //! \end{pmatrix},
    //! \qquad
    //! (\mathbf{A}\verb|-=|\mathbf{a}+2*\mathbf{b}) = \begin{pmatrix}
    //!     -1 & -1 &    &    \\
    //!      1 & 2  & -1 &    \\
    //!        & 1  & -5 & -1 \\
    //!        &    & 1  & 2  \\
    //! \end{pmatrix}
    //! \f]
    //!
    template <typename VEC_EXPR>
    inline
    TridMatrix<T> &
    operator -= ( VEC_EXPR const & e )
    { TMP = e; D -= TMP; return *this; }

    //@}

    real_type const &
    operator () ( integer i, integer j ) const {
      SPARSE::test_index(i,j);
      integer kk = ((j+1)-i);
      UTILS_ASSERT(
        kk < 3,
        "TridMatrix({},{}) index out of range\n", i, j
      );
      return A( kk*SPARSE::sp_nrows+i-1);
    }

    real_type &
    operator () ( integer i, integer j ) {
      SPARSE::test_index(i,j);
      integer kk = ((j+1)-i);
      UTILS_ASSERT(
        kk < 3,
        "TridMatrix({},{}) index out of range\n", i, j
      );
      return A( kk*SPARSE::sp_nrows+i-1);
    }

    real_type const &
    value( integer i, integer j ) const {
      static const T zero(0);
      SPARSE::test_index(i,j);
      integer kk = ((j+1)-i);
      if ( kk < 3 ) return A( kk*SPARSE::sp_nrows+i-1);
      else          return zero;
    }

    bool
    exists( integer i, integer j ) {
      return i < SPARSE::sp_nrows &&
             j < SPARSE::sp_nrows &&
             (i <= j+1 || j <= i+1 );
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<real_type> const & getA(void) const { return A; } //!< return the value vector
    Vector<real_type>       & getA(void)       { return A; } //!< return the value vector
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    void
    Begin(void) const {
      iter_counter = 0;
      iter_row     = 0;
      iter_col     = 0;
    }

    void
    Next(void) const {
      ++iter_counter;
      iter_row = (iter_counter+1) % SPARSE::sp_nrows;
      iter_col = ( (iter_counter+1) / SPARSE::sp_nrows) +iter_row-1;
    }

    bool
    End(void) const
    { return iter_counter < SPARSE::sp_nnz; }

    integer         row    (void) const { return iter_row; }
    integer         column (void) const { return iter_col; }
    real_type const & value  (void) const { return A(iter_counter); }
    real_type       & value  (void)       { return A(iter_counter); }

    //!
    //! Assign the pointed element to `rhs`.
    //!
    template <typename TS>
    inline
    void
    assign( TS & rhs ) const { rhs = A(iter_counter); }

    //!
    //! Perform the operation `res += s * (A * x)`
    //!
    template <typename VEC>
    inline
    void
    add_S_mul_M_mul_V(
      VEC       & res,
      T   const & s,
      VEC const & x
    ) const {
      UTILS_ASSERT0(
        std::addressof(res) != std::addressof(x),
        "Sparse_tool: [res += s * (A * x)] add_S_mul_M_mul_V, `res` and `x` cant be the same\n"
      );
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_nrows,
        "Sparse_tool: [res += s * (A * x)] result vector too small\n"
      );
      TMP = s*x;
      res += D.array() * TMP.array();
      integer nn = SPARSE::sp_nrows-1;
      res.head(nn).array() += U.array() * TMP.tail(nn).array();
      res.tail(nn).array() += L.array() * TMP.head(nn).array();
    }

    //!
    //! Perform the operation `res += s * (A^T * x)`.
    //!
    template <typename VEC>
    inline
    void
    add_S_mul_Mt_mul_V(
      VEC       & res,
      T   const & s,
      VEC const & x
    ) const {
      UTILS_ASSERT0(
        std::addressof(res) != std::addressof(x),
        "Sparse_tool: [res += s * (A^T * x)] add_S_mul_Mt_mul_V, `res` and `x` cant be the same\n"
      );
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_nrows,
        "Sparse_tool: [res += s * (A^T * x)] result vector too small\n"
      );
      res += s * D * x;
      res += s * L * x;
      SPARSELIB_LOOP( SPARSE::sp_nrows - 1, res(i+1) += s * U(i) * x(i+1) );
    }

    //!
    //! Perform the operation `res += s * (A * x)`
    //!
    template <typename VRES, typename VB>
    inline
    void
    add_S_mul_M_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_nrows,
        "[res += s * (A * x)] result vector too small\n"
      );
      res += s * D * x;
      res += s * U * x;
      SPARSELIB_LOOP( SPARSE::sp_nrows - 1, res(i+1) += s * L(i) * x(i+1) );
    }

    //!
    //! Perform the operation `res += s * (A^T * x)`.
    //!
    template <typename VRES, typename VB>
    inline
    void
    add_S_mul_Mt_mul_V(
      VRES     & res,
      T  const & s,
      VB const & x
    ) const {
      UTILS_ASSERT0(
        res.size() >= SPARSE::sp_nrows,
        "[res += s * (A^T * x)] result vector too small\n"
      );
      res += s * D * x;
      res += s * L * x;
      SPARSELIB_LOOP( SPARSE::sp_nrows - 1, res(i+1) += s * U(i) * x(i+1) );
    }

    //!
    //! solve the tridiagonal system `A*x=b`.
    //!
    template <typename VECTOR>
    inline
    void
    ass_V_div_M(VECTOR & res, VECTOR const & b) const {
      solve(b, res);
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  SPARSELIB_MUL_STRUCTURES(TridMatrix)

  template <typename TM, typename T> inline
  Vector_V_div_M<Vector<T>,TridMatrix<TM> >
  operator / (Vector<T> const & a, TridMatrix<TM> const & M) {
    return Vector_V_div_M<Vector<T>,TridMatrix<TM> >(a,M);
  }
  #endif

  /*
   *  ###       # #######
   *   #       #  #     #
   *   #      #   #     #
   *   #     #    #     #
   *   #    #     #     #
   *   #   #      #     #
   *  ### #       #######
   */

  //!
  //! \name I/O
  //!
  //@{

  //!
  //! Print the contents of vector `v` to the stream `s`.
  //!
  template <typename T, typename VEC_EXPR>
  inline
  string
  to_string( VEC_EXPR const & v ) {
    integer sz1 = v.size() - 1;
    string res = "[ ";
    for ( integer i{0}; i < sz1; ++i ) res += fmt::format( "{}, ", v(i) );
    res += fmt::format( "{} ]", v(sz1) );
    return res;
  }

  //!
  //! Print the contents of vector `v` to the stream `s`.
  //!
  template <typename T, typename VEC_EXPR> inline
  ostream_type &
  operator << ( ostream_type & s, VEC_EXPR const & v ) {
    s << v.to_string();
    return s;
  }

  //!
  //! Read from the stream `s` to the vector `v`.
  //!
  template <typename T> inline
  istream_type &
  operator >> ( istream_type & s, Vector<T> & v ) {
    typename Vector<T>::pointer pv=v.begin();
    while ( pv != v.end() ) s >> *pv++;
    return s;
  }

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  //!
  //! Print the contents of the SparsePattern `S` to the stream `s`.
  //!
  inline
  ostream_type &
  operator << ( ostream_type & s, SparsePattern const & S ) {
    for ( S.Begin(); S.End(); S.Next() )
      s << S.row() << " " << S.column() << "\n";
    return s;
  }
  #endif

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  //!
  //! Print the contents of the `Sparse<T,MAT>`
  //! object `N` to the stream `s`.
  //!
  template <typename T, typename MAT> inline
  ostream_type &
  operator << ( ostream_type & s, Sparse<T,MAT> const & M ) {
    fmt::print( s,
      "size    = {} x {}\n"
      "nnz     = {}\n",
      M.nrows(), M.ncols(), M.nnz()
    );
    if ( M.is_ordered() ) {
      integer l,d,u;
      M.nnz( l, d, u );
      fmt::print( s,
        "  (diag={}, lower={}, upper={}) ordered = YES\n",
        d, l, u
       );
    } else {
      s << "ordered = NO\n";
    }
    for ( M.Begin(); M.End(); M.Next() )
      fmt::print( s, "({},{}) = {}\n", M.row(), M.column(), M.value() );
    return s;
  }
  #endif
  //@}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {

  using ::Sparse_tool::Vector;
  using ::Sparse_tool::Sparse;
  using ::Sparse_tool::SparsePattern;
  using ::Sparse_tool::CCoorMatrix;
  using ::Sparse_tool::CRowMatrix;
  using ::Sparse_tool::CColMatrix;
  using ::Sparse_tool::TridMatrix;

  // I/O
  using ::Sparse_tool::operator >>;
  using ::Sparse_tool::operator <<;
}
#endif

namespace fmt {
  template <> struct formatter<Sparse_tool::SparsePattern> : ostream_formatter {};
  template <typename TYPE, typename MAT> struct formatter<Sparse_tool::Sparse<TYPE,MAT>> : ostream_formatter {};
}

#endif

/*
// ####### ####### #######
// #       #     # #
// #       #     # #
// #####   #     # #####
// #       #     # #
// #       #     # #
// ####### ####### #
*/
