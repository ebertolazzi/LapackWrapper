/*!

  \file     SparseTool.hh
  \mainpage SparseTool: a Sparse Matrix Manager 
  \date     2011, July 21
  \version  1.1
  \note     based on SparseLib 11.0 (C 1999 -- 2008)

  \author Enrico Bertolazzi
 
  \par Affiliations:
       Dipartimento di Ingegneria Industriale<BR>
       Universita` degli Studi di Trento<BR>
       email: enrico.bertolazzi@unitn.it<BR>
 
  \section Preface
  \a SparseTool is a collection of simple, fast, and efficient classes for
  manipulating large vectors and large sparse matrices.
  
  \section The Basic Vector and Matrix Classes
  
  The library consist of the following templated classes:
  - \c Vector\<T\> which define a \b dense large column vector.
  - \c SparsePattern which define a \b pattern of the nonzero 
    elements which can be used to construct a sparse matrix.
  - \c TridMatrix\<T\> which implements a
    \b tridiagonal matrix.
  - \c CCoorMatrix\<T\> which implements a sparse
    <b> Compressed Coordinate Storage </b> matrix.
  - \c CRowMatrix\<T\> which implements a sparse
    <b> Compressed Rows Storage </b> matrix.
  - \c CColMatrix\<T\> which implements a sparse
    <b> Compressed Columns Storage </b> matrix.

  Those classes allow vectors and matrices to be formally treated in
  software implementations as mathematical objects in arithmetic
  expressions.  For example if \c A is a sparse matrix and \c b
  is a \c Vector\<T\> than \c A*b means the matrix-vector product.

  \section Loading the library
  To use the library you must include it by the following piece of code:

\code
#include "SparseTool.hh"
using namespace SparseToolLoad;
\endcode

  line \a 2 is recommended to avoid \c SparseTool:: prefix
  for the library call.

  \par Some Macros for customization
  To acivate some control when \c SparseTool is used in developing

\code
  #define SPARSETOOL_DEBUG
\endcode

*/

/*!
   \defgroup Promote Promotion Classes and Structures
 */

#ifndef SPARSETOOL_HH
#define SPARSETOOL_HH

#include <cmath>

// standard includes I/O
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utility>

// STL lib
#include <vector>
#include <algorithm>
#include <functional>

#include <complex>
#include <string>

#include <cstdint>

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

//! default number of pre-allocated nonzeros
#ifndef SPARSETOOL_DEFAULT_NNZ
  #define SPARSETOOL_DEFAULT_NNZ 100
#endif

//! issue an error message
#define SPARSETOOL_ERR(W)                             \
  { using namespace ::std;                            \
    cerr << "\n" << W  << "\nin file `" << __FILE__   \
         << "' on line " << __LINE__ << "\n";;        \
    exit(0); }

/*! \cond NODOC */

#define SPARSETOOL_ASSERT(X,W) if ( !(X) ) SPARSETOOL_ERR(W)

#ifdef SPARSETOOL_DEBUG
  #define SPARSETOOL_INLINE
  #define SPARSETOOL_TEST(X,W) SPARSETOOL_ASSERT(X,W)
#else
  #define SPARSETOOL_INLINE inline
  #define SPARSETOOL_TEST(X,W)
#endif

// loops
#define SPARSELIB_LOOP(N,DO)   { for ( indexType i = 0; i < N; ++i ) { DO; } }
#define SPARSELIB_V1LOOP(DO)   SPARSELIB_LOOP(this->size(),DO)
#define SPARSELIB_V2LOOP(V,DO) SPARSELIB_LOOP(minIndex(V.size(),this->size()),DO)

/*! \endcond */

/*! \cond NODOC */
namespace SparseToolFun {

  static inline float       conj(float       const & a) { return a; }
  static inline double      conj(double      const & a) { return a; }
  static inline long double conj(long double const & a) { return a; }

  static inline float       absval(float       const a) { return a > 0 ? a : -a; }
  static inline double      absval(double      const a) { return a > 0 ? a : -a; }
  static inline long double absval(long double const a) { return a > 0 ? a : -a; }

  static inline float       absval(std::complex<float>       const & a) { return std::abs(a); }
  static inline double      absval(std::complex<double>      const & a) { return std::abs(a); }
  static inline long double absval(std::complex<long double> const & a) { return std::abs(a); }

  static inline float       maxval(float       const a, float       const b) { return a > b ? a : b; }
  static inline double      maxval(double      const a, double      const b) { return a > b ? a : b; }
  static inline long double maxval(long double const a, long double const b) { return a > b ? a : b; }

  static inline float       minval(float       const a, float       const b) { return a < b ? a : b; }
  static inline double      minval(double      const a, double      const b) { return a < b ? a : b; }
  static inline long double minval(long double const a, long double const b) { return a < b ? a : b; }

  static inline float       maxabsval(float       const a, float       const b) { return std::max(absval(a),absval(b)); }
  static inline double      maxabsval(double      const a, double      const b) { return std::max(absval(a),absval(b)); }
  static inline long double maxabsval(long double const a, long double const b) { return std::max(absval(a),absval(b)); }

  static inline float       minabsval(float       const a, float       const b) { return std::min(absval(a),absval(b)); }
  static inline double      minabsval(double      const a, double      const b) { return std::min(absval(a),absval(b)); }
  static inline long double minabsval(long double const a, long double const b) { return std::min(absval(a),absval(b)); }

  template <typename T> inline T absval2(T const & a) { return a*a; }
  template <typename T> inline T absval2(std::complex<T> const & a) { T bf(absval(a)); return bf*bf; }

  using namespace ::std;

}
/*! \endcond */

//! The namespace with the SparseTool toolkit 
namespace SparseTool {

  using ::std::vector;
  using ::std::istream;
  using ::std::ostream;
  using ::std::cin;
  using ::std::cout;
  using ::std::cerr;
  using ::std::setw;
  using ::std::greater;
  using ::std::less;
  
  // FROM BLITZ++

  /*! \cond NODOC */

  /*! \brief
   * This class is copied from \c BLITZ++.
   * The C++ compiler supports partial specialization, so type promotion
   * can be done the elegant way.
   * This implementation is after ideas by Jean-Louis Leroy.
   */

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

  /*! \endcond */

  /*! \cond NODOC */

  template<typename T>
  struct return_trait {
    typedef T valueType; //!< type of the trait
  };

  // Partial specialization for \c complex type.
  template<typename T> struct return_trait<std::complex<T> > { typedef T valueType; };

  //! \internal Class for type promotion
  template<typename T>
  struct autopromote_trait { typedef T T_numtype; };

  template<> struct autopromote_trait<bool>
  { typedef int T_numtype; }; // promote \c bool to \c int
  
  template<> struct autopromote_trait<char>
  { typedef int T_numtype; }; // promote \c char to \c int

  template<> struct autopromote_trait<unsigned char>
  { typedef int T_numtype; }; // promote \c bool to \c int

  template<> struct autopromote_trait<short>
  { typedef int T_numtype; }; // promote \c short to \c int
  
  template<> struct autopromote_trait<unsigned short>
  { typedef unsigned T_numtype; }; // promote \c unsigned \c short to \c unsigned

  template<typename T1, typename T2, bool promoteToT1>
  struct promote2 { typedef T1 T_promote; };

  template<typename T1, typename T2>
  struct promote2<T1,T2,false> { typedef T2 T_promote; };

  /*! \endco; */

  /*! \cond NODOC */

  template<typename T1_orig, typename T2_orig>
  class promote_trait {
  private:
    //! Handle promotion of small integers to int/unsigned int
    typedef typename autopromote_trait<T1_orig>::T_numtype T1;
    //! Handle promotion of small integers to int/unsigned int
    typedef typename autopromote_trait<T2_orig>::T_numtype T2;

    //! True if T1 is higher ranked
    static bool const T1IsBetter = precision_trait<T1>::precisionRank >
                                   precision_trait<T2>::precisionRank;

    //! True if we know ranks for both T1 and T2
    static bool const knowBothRanks = precision_trait<T1>::knowPrecisionRank
                                   && precision_trait<T2>::knowPrecisionRank;

    //! True if we know T1 but not T2
    static bool const knowT1butNotT2 = precision_trait<T1>::knowPrecisionRank &&
                                      !precision_trait<T2>::knowPrecisionRank;

    //! True if we know T2 but not T1
    static bool const knowT2butNotT1 = precision_trait<T2>::knowPrecisionRank &&
                                      !precision_trait<T1>::knowPrecisionRank;

    //! True if T1 is bigger than T2
    static bool const T1IsLarger = sizeof(T1) >= sizeof(T2);

    /*!
     * - We know T1 but not T2: true
     * - We know T2 but not T1: false
     * - Otherwise, if T1 is bigger than T2: true
     */
    static bool const defaultPromotion =
      knowT1butNotT2 ? false : (knowT2butNotT1 ? true : T1IsLarger);

    /*!
     * - If we have both ranks, then use them.
     * - If we have only one rank, then use the unknown type.
     * - If we have neither rank, then promote to the larger type.
     */
    static bool const promoteToT1 = knowBothRanks ? T1IsBetter : defaultPromotion;

  public:
    //! promoted type
    typedef typename promote2<T1,T2,promoteToT1>::T_promote T_promote;
  };
  
  //! Class for promotion with three types
  template<typename T1_orig, typename T2_orig, typename T3_orig>
  class promote_trait3 {
    //! Handle promotion of small integers to int/unsigned int
    typedef typename promote_trait<T1_orig,T2_orig>::T_promote T12_promote;
  public:
    //! promoted type
    typedef typename promote_trait<T12_promote,T3_orig>::T_promote T_promote;
  };
    
  #ifndef SPARSELIB_FUNCTION_TYPE
  #define SPARSELIB_FUNCTION_TYPE float
  #endif

  #ifndef SPARSELIB_ACCUMULATOR_TYPE
  #define SPARSELIB_ACCUMULATOR_TYPE double
  #endif

  #define SPARSELIB_PROMOTE(Ta,Tb) promote_trait<Ta,Tb>::T_promote
  #define SPARSELIB_PROMOTE3(Ta,Tb,Tc) \
    promote_trait<Ta,typename promote_trait<Tb,Tc>::T_promote>::T_promote

  #define SPARSELIB_ACCUMULATOR_PROMOTE(T) \
    SPARSELIB_PROMOTE(SPARSELIB_ACCUMULATOR_TYPE,T)

  #define SPARSELIB_ACCUMULATOR_PROMOTE2(Ta,Tb) \
    SPARSELIB_PROMOTE3(SPARSELIB_ACCUMULATOR_TYPE,Ta,Tb)

  #define SPARSELIB_TYPES_FROM_TYPENAME(A) \
    typedef typename A::valueType valueType

  #define SPARSELIB_TYPES_FROM_TYPENAME2(A,B)                                \
    typedef typename                                                         \
    promote_trait<typename A::valueType,typename B::valueType>::T_promote    \
    valueType

  #define SPARSELIB_TYPES_FROM_TYPENAME3(S,A)                                \
    typedef typename promote_trait<S,typename A::valueType>::T_promote       \
    valueType

  #define SPARSELIB_TYPES_FROM_TYPENAME_FUN1(A) \
  SPARSELIB_TYPES_FROM_TYPENAME3(SPARSELIB_FUNCTION_TYPE,A)

  #define SPARSELIB_TYPES_FROM_TYPENAME_FUN2(A,B)                            \
    typedef typename                                                         \
      promote_trait<typename A::valueType,typename B::valueType>::T_promote  \
      promoted_valueType;                                                    \
    typedef typename promote_trait<SPARSELIB_FUNCTION_TYPE,                  \
                                   promoted_valueType>::T_promote valueType

  template <typename T> class Vector;

  /*! \endcond */

  //! define \c uint32_t as the type for indexing vector and matrices 
  typedef uint32_t indexType;

  //! return minimum value between \c a and \c b
  inline indexType
  minIndex(indexType a, indexType b)
  { return a < b ? a : b; }

  //! return maximum value between \c a and \c b
  inline indexType
  maxIndex(indexType a, indexType b)
  { return a > b ? a : b; }
  
  //! \defgroup Comparator Comparator for Sparse Matrix costructor and conversion 
  //@{

  //! comparator class for selecting \a all elements
  struct all_ok {
    //! return always true
    bool operator () ( indexType i, indexType j ) const { return true; }
  };

  //! comparator class for selecting;lements \a under the diagonal
  struct lower_ok {
    //! select lower diagonal indexes
    bool operator () ( indexType i, indexType j ) const { return i > j; }
  };

  //! comparator class for selecting elements \a under the diagonal
  struct lowereq_ok {
    //! select lower diagonal indexes
    bool operator () ( indexType i, indexType j ) const { return i >= j; }
  };

  //! comparator class for selec;ng elements \a over the diagonal
  struct upper_ok {
    //! select upper diagonal indexes
    bool operator () ( indexType i, indexType j ) const { return i < j; }
  };

  //! comparator class for selecting elements \a over the diagonal
  struct uppereq_ok {
    //! select upper diagonal indexes
    bool operator () ( indexType i, indexType j ) const { return i <= j; }
  };

  //! comparator class for selecting elements \a on the diagonal
  struct diag_ok {
    //! select diagonal indexes
    bool operator () ( indexType i, indexType j ) const { return i==j; }
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
  
  //! \cond NODOC
  
  //! \defgroup MVStructures Structures Storing Matrix/Vector operations
  //@{

  //! structure storing the operation \c a/M  (vector-matrix division i.e. solve \c A*x=a)
  template <typename VA, typename MATRIX>
  struct Vector_V_div_M {
    VA     const & a;
    MATRIX const & M;
    Vector_V_div_M( VA const & _a, MATRIX const & _M )
    : a(_a)
    , M(_M)
    {}
  };

  //! structure storing the operation \c M*a (matrix-vector multiply)
  template <typename MATRIX, typename VA>
  struct Vector_M_mul_V {
    MATRIX const & M;
    VA     const & a;
    Vector_M_mul_V( MATRIX const & _M, VA const & _a )
    : M(_M)
    , a(_a)
    {}
  };

  //! structure storing the operation \c s*M  (scalar-matrix multiply)
  template <typename S, typename MATRIX>
  struct Vector_S_mul_M {
    S      const & s;
    MATRIX const & M;
    Vector_S_mul_M( S const & _s, MATRIX const & _M )
    : s(_s)
    , M(_M)
    {}
  };

  //! structure storing the operation \c M^a  (transpose matrix-vector multiply)
  template <typename MATRIX, typename VA>
  struct Vector_Mt_mul_V {
    MATRIX const & M;
    VA     const & a;
    Vector_Mt_mul_V( MATRIX const & _M, VA const & _a )
    : M(_M)
    , a(_a)
    {}
  };

  //! structure storing the operation \c s*(M*a)  (scalar-matrix-vector multiply)
  template <typename SCALAR, typename MATRIX, typename VA>
  struct Vector_S_mul_M_mul_V {
    SCALAR const & s;
    MATRIX const & M;
    VA     const & a;
    Vector_S_mul_M_mul_V(
      SCALAR const & _s,
      MATRIX const & _M,
      VA     const & _a
    )
    : s(_s)
    , M(_M)
    , a(_a)
    {}
  };

  //! structure storing the operation \c s*(M^a)  (transpose scalar-matrix-vector multiply)
  template <typename SCALAR, typename MATRIX, typename VA>
  struct Vector_S_mul_Mt_mul_V {
    SCALAR const & s;
    MATRIX const & M;
    VA     const & a;
    Vector_S_mul_Mt_mul_V(
      SCALAR const & _s,
      MATRIX const & _M,
      VA     const & _a
    )
    : s(_s)
    , M(_M)
    , a(_a)
    {}
  };

  //! structure storing the operation \c a+M*b
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sum_M_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_M_mul_V(
      VA     const & _a,
      MATRIX const & _M,
      VB     const & _b
    )
    : a(_a)
    , M(_M)
    , b(_b) 
    {}
  };

  //! structure storing the operation \c a+M^b
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sum_Mt_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_Mt_mul_V(
      VA     const & _a,
      MATRIX const & _M,
      VB     const & _b
    )
    : a(_a)
    , M(_M)
    , b(_b) 
    {}
  };

  //! structure storing the operation \c a-M*b
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sub_M_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sub_M_mul_V(
      VA     const & _a,
      MATRIX const & _M,
      VB     const & _b
    )
    : a(_a)
    , M(_M)
    , b(_b)
    {}
  };

  //! structure storing the operation \c a+M^b
  template <typename VA, typename MATRIX, typename VB>
  struct Vector_V_sub_Mt_mul_V {
    VA     const & a;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sub_Mt_mul_V(
      VA     const & _a,
      MATRIX const & _M,
      VB     const & _b
    )
    : a(_a)
    , M(_M)
    , b(_b)
    {}
  };
  
  //! structure storing the operation \c a+s*M*b
  template <typename VA, typename SCALAR, typename MATRIX, typename VB>
  struct Vector_V_sum_S_mul_M_mul_V {
    VA     const & a;
    SCALAR const & s;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_S_mul_M_mul_V(
      VA     const & _a,
      SCALAR const & _s,
      MATRIX const & _M,
      VB     const & _b
    )
    : a(_a)
    , s(_s)
    , M(_M)
    , b(_b)
    {}
  };

  //! structure storing the operation \c a+s*M^b
  template <typename VA, typename SCALAR, typename MATRIX, typename VB>
  struct Vector_V_sum_S_mul_Mt_mul_V {
    VA     const & a;
    SCALAR const & s;
    MATRIX const & M;
    VB     const & b;
    Vector_V_sum_S_mul_Mt_mul_V(
      VA     const & _a,
      SCALAR const & _s,
      MATRIX const & _M,
      VB     const & _b
    )
    : a(_a)
    , s(_s)
    , M(_M)
    , b(_b)
    {}
  };

  //! structure storing the operation \c a/P, i.e. the application of Preconditioner
  template <typename VECTOR, typename PRECO>
  struct Vector_V_div_P {
    VECTOR const & a;
    PRECO  const & P;
    Vector_V_div_P( VECTOR const & _a, PRECO const & _P )
    : a(_a)
    , P(_P)
    {}
  };

  /*! \brief Template class storing (recursively) vector expression */
  template <typename T, typename A>
  class VectorE {
  public:
   typedef T valueType;

  private:
    A expr;

  public:
    VectorE(A const & a) : expr(a) { }
    valueType operator () (indexType i) const { return expr(i); }
    indexType size(void) const { return expr.size(); }
  };
  
  //@}
  
  //! \endcond


//! macro to define the standar matrix/vector operators
#define SPARSELIB_VECTOR_OPERATIONS(VECTOR) \
    /*! \brief Assign the vector \c v derived from \c VectorBase to \c *this */ \
    template <typename VEC> inline \
    VECTOR const & \
    operator = (VectorBase<T,VEC> const & v) { \
      SPARSELIB_V2LOOP(v, (*this)(i) = v(i)); \
      return *this; \
    } \
    /*! \brief  Assign the vector expression \c e to \c *this */ \
    template <typename R> inline \
    VECTOR const & \
    operator = (VectorE<T,R> const & e) { \
      SPARSELIB_V2LOOP(e, (*this)(i) = e(i)); \
      return *this; \
    } \
    /*! \brief  Fill the Vector with the constant \a s. */ \
    VECTOR const & \
    operator = (valueType const & s) { \
      SPARSELIB_V1LOOP( (*this)(i) = valueType(s) ); \
      return *this; \
    } \
    /*! \brief  Add the element of the \c v derived from \c VectorBase to \c *this. 
        If the vector are not of the same size then \c min(size(),v.size()) elements are copied. */ \
    template <typename VEC> inline \
    VECTOR const & \
    operator += (VectorBase<T,VEC> const & v) { \
      SPARSELIB_V2LOOP(v, (*this)(i) += v(i)); \
      return *this; \
    } \
    /*! \brief  Add the expression \c e to \c *this.
        If the vector are not of the same size then \c min(size(),e.size()) elements are evaluated. */ \
    template <typename R> inline \
    VECTOR const & \
    operator += (VectorE<T,R> const & e) { \
      SPARSELIB_V2LOOP(e, (*this)(i) += e(i)); \
      return *this; \
    } \
    /*! \brief  Add to each element of the vector \c *this the constant \c s. */ \
    VECTOR const & \
    operator += (valueType const & s) { \
      SPARSELIB_V1LOOP( (*this)(i) += valueType(s)); \
      return *this; \
    } \
    /*! \brief  Subtract the element of the Vector \c v to \c *this.
        If the vector are not of the same size then \c min(size(),v.size()) elements are copied. */ \
    template <typename VEC> inline \
    VECTOR const & \
    operator -= (VectorBase<T,VEC> const & v) { \
      SPARSELIB_V2LOOP(v, (*this)(i) -= v(i)); \
      return *this; \
    } \
    /*! \brief  Subtract the expression \c e to \c *this.
        If the vector are not of the same size then \c min(size(),e.size()) elements are evaluated. */ \
    template <typename R> inline \
    VECTOR const & \
    operator -= (VectorE<T,R> const & e) { \
      SPARSELIB_V2LOOP(e, (*this)(i) -= e(i)); \
      return *this; \
    } \
    /*! \brief  Subtract to each element of the vector \c *this the constant \c s. */ \
    VECTOR const & \
    operator -= (valueType const & s) { \
      SPARSELIB_V1LOOP( (*this)(i) -= valueType(s)); \
      return *this; \
    } \
    /*! \brief  Multiply the element of the Vector \c v to \c *this.
        If the vector are not of the same size then \c min(size(),v.size()) elements are copied. */ \
    template <typename VEC> inline \
    VECTOR const & \
    operator *= (VectorBase<T,VEC> const & v) { \
      SPARSELIB_V2LOOP(v, (*this)(i) *= v(i)); \
      return *this; \
    } \
    /*! \brief  Multiply the Vector expression \c e to \c *this.
        If the vector and expression are not of the same size then \c min(size(),e.size()) elements are copied. */ \
    template <typename R> inline \
    VECTOR const & \
    operator *= (VectorE<T,R> const & e) { \
      SPARSELIB_V2LOOP(e, (*this)(i) *= e(i)); \
      return *this; \
    } \
    /*! \brief  Multiply each element of the Vector \c *this by constant \a s. */ \
    VECTOR const & \
    operator *= (valueType const & s) { \
      SPARSELIB_V1LOOP( (*this)(i) *= valueType(s)); \
      return *this; \
    } \
    /*! \brief  Divide the element of the Vector \c *this to the components of vector \c v.
        If the vector are not of the same size then \c min(size(),v.size()) elements are divided. */ \
    template <typename VEC> inline \
    VECTOR const & \
    operator /= (VectorBase<T,VEC> const & v) { \
      SPARSELIB_V2LOOP(v, (*this)(i) /= v(i)); \
      return *this; \
    } \
    /*! \brief  Divide the element of the Vector \c *this to the components of expression \c e.
        If the vector and the expression are not of the same size then \c min(size(),e.size()) elements are divided. */ \
    template <typename R> inline \
    VECTOR const & \
    operator /= (VectorE<T,R> const & e) { \
      SPARSELIB_V2LOOP(e, (*this)(i) /= e(i)); \
      return *this; \
    } \
    /*! \brief  Divide each element of the Vector by the constant \c s. */ \
    VECTOR const & \
    operator /= (valueType const & s) { \
      SPARSELIB_V1LOOP( (*this)(i) /= valueType(s)); \
      return *this; \
    } \
    /*------------ MM_COMPLEX MATRIX OPERATIONS ---------------------*/ \
    /*! \brief  Assign to \c *this the evaluation of the expression \c a/M contained in \c op of type \c Vector_V_div_M */ \
    template <typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator = ( Vector_V_div_M<VA,MATRIX> const & op ) \
    { op.M.ass_V_div_M(*this, op.a); return *this; } \
    \
    /*! \brief  Assign to \c *this the evaluation of the expression \c M*a contained in \c op of type \c Vector_M_mul_V */ \
    template <typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator = ( Vector_M_mul_V<MATRIX,VA> const & op ) { \
      *this = valueType(0); \
      op.M.add_S_mul_M_mul_V(*this, valueType(+1), op.a); \
      return *this; \
    } \
    /*! \brief  Add to \c *this  the evaluation of the expression \c a/M contained in \c op of type \c Vector_V_div_M */ \
    template <typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator += ( Vector_M_mul_V<MATRIX,VA> const & op ) { \
      op.M.add_S_mul_M_mul_V(*this, valueType(+1), op.a); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c M*a contained in \c op of type \c Vector_M_mul_V */ \
    template <typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator -= ( Vector_M_mul_V<MATRIX,VA> const & op ) { \
      op.M.add_S_mul_M_mul_V(*this, valueType(-1), op.a); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c M^a contained in \c op of type \c Vector_Mt_mul_V */ \
    template <typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator = (Vector_Mt_mul_V<MATRIX,VA> const & op) { \
      *this = valueType(0); \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(+1), op.a); \
      return *this; \
    } \
    /*! \brief  Add to \c *this  the evaluation of the expression \c M^a contained in \c op of type \c Vector_Mt_mul_V */ \
    template <typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator += (Vector_Mt_mul_V<MATRIX,VA> const & op) { \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(+1), op.a); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c M^a contained in \c op of type \c Vector_Mt_mul_V */ \
    template <typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator -= (Vector_Mt_mul_V<MATRIX,VA> const & op) { \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(-1), op.a); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c s*(M*a) contained in \c op of type \c Vector_S_mul_M_mul_V */ \
    template <typename SCALAR, typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator = (Vector_S_mul_M_mul_V<SCALAR,MATRIX,VA> const & op) { \
      *this = valueType(0); \
      op.M.add_S_mul_M_mul_V(*this, op.s, op.a); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c s*(M*a) contained in \c op of type \c Vector_S_mul_M_mul_V */ \
    template <typename SCALAR, typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator += (Vector_S_mul_M_mul_V<SCALAR,MATRIX,VA> const & op) { \
      op.M.add_S_mul_M_mul_V(*this, op.s, op.a); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c s*(M*a) contained in \c op of type \c Vector_S_mul_M_mul_V */ \
    template <typename SCALAR, typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator -= (Vector_S_mul_M_mul_V<SCALAR,MATRIX,VA> const & op) { \
      op.M.add_S_mul_M_mul_V(*this, -op.s, op.a); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c s*(M^a) contained in \c op of type \c Vector_S_mul_Mt_mul_V */ \
    template <typename SCALAR, typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator = (Vector_S_mul_Mt_mul_V<SCALAR,MATRIX,VA> const & op) { \
      *this = valueType(0); \
      op.M.add_S_mul_Mt_mul_V(*this, op.s, op.a); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c s*(M^a) contained in \c op of type \c Vector_S_mul_Mt_mul_V */ \
    template <typename SCALAR, typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator += (Vector_S_mul_Mt_mul_V<SCALAR,MATRIX,VA> const & op) { \
      op.M.add_S_mul_Mt_mul_V(*this, op.s, op.a); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c s*(M^a) contained in \c op of type \c Vector_S_mul_Mt_mul_V */ \
    template <typename SCALAR, typename MATRIX, typename VA> inline \
    VECTOR const & \
    operator -= (Vector_S_mul_Mt_mul_V<SCALAR,MATRIX,VA> const & op) { \
      op.M.add_S_mul_Mt_mul_V(*this, -op.s, op.a); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c a+M*b  contained in \c op of type \c Vector_V_sum_M_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator = (Vector_V_sum_M_mul_V<VA,MATRIX,VB> const & op) { \
      *this = op.a; \
      op.M.add_S_mul_M_mul_V(*this, valueType(+1), op.b); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c a+M*b  contained in \c op of type \c Vector_V_sum_M_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator += (Vector_V_sum_M_mul_V<VA,MATRIX,VB> const & op) { \
      *this += op.a; \
      op.M.add_S_mul_M_mul_V(*this, valueType(+1), op.b); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c a+M*b  contained in \c op of type \c Vector_V_sum_M_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator -= (Vector_V_sum_M_mul_V<VA,MATRIX,VB> const & op) { \
      *this -= op.a; \
      op.M.add_S_mul_M_mul_V(*this, valueType(-1), op.b); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c a+M^b  contained in \c op of type \c Vector_V_sum_Mt_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator = (Vector_V_sum_Mt_mul_V<VA,MATRIX,VB> const & op) { \
      *this = op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(+1), op.b); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c a+M^b  contained in \c op of type \c Vector_V_sum_Mt_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator += (Vector_V_sum_Mt_mul_V<VA,MATRIX,VB> const & op) { \
      *this += op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(+1), op.b); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c a+M^b  contained in \c op of type \c Vector_V_sum_Mt_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator -= (Vector_V_sum_Mt_mul_V<VA,MATRIX,VB> const & op) { \
      *this -= op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(-1), op.b); \
      return *this; \
    } \
    /*-----------------------*/ \
    /*! \brief  Assign to \c *this the evaluation of the expression \c a-M*b  contained in \c op of type \c Vector_V_sub_M_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator = (Vector_V_sub_M_mul_V<VA,MATRIX,VB> const & op) { \
      *this = op.a; \
      op.M.add_S_mul_M_mul_V(*this, valueType(-1), op.b); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c a-M*b  contained in \c op of type \c Vector_V_sub_M_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator += (Vector_V_sub_M_mul_V<VA,MATRIX,VB> const & op) { \
      *this += op.a; \
      op.M.add_S_mul_M_mul_V(*this, valueType(-1), op.b); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c a-M*b  contained in \c op of type \c Vector_V_sub_M_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator -= (Vector_V_sub_M_mul_V<VA,MATRIX,VB> const & op) { \
      *this -= op.a; \
      op.M.add_S_mul_M_mul_V(*this, valueType(+1), op.b); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c a-M^b  contained in \c op of type \c Vector_V_sub_Mt_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator = (Vector_V_sub_Mt_mul_V<VA,MATRIX,VB> const & op) { \
      *this = op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(-1), op.b); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c a-M^b  contained in \c op of type \c Vector_V_sub_Mt_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator += (Vector_V_sub_Mt_mul_V<VA,MATRIX,VB> const & op) { \
      *this += op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(-1), op.b); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c a-M^b  contained in \c op of type \c Vector_V_sub_Mt_mul_V */ \
    template <typename VA, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator -= (Vector_V_sub_Mt_mul_V<VA,MATRIX,VB> const & op) { \
      *this -= op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, valueType(+1), op.b); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c  a+s*(M*b)  contained in \c op of type \c Vector_V_sum_S_mul_M_mul_V */ \
    template <typename VA, typename SCALAR, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator = (Vector_V_sum_S_mul_M_mul_V<VA,SCALAR,MATRIX,VB> const & op) { \
      *this = op.a; \
      op.M.add_S_mul_M_mul_V(*this, op.s, op.b); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c  a+s*(M*b)  contained in \c op of type \c Vector_V_sum_S_mul_M_mul_V */ \
    template <typename VA, typename SCALAR, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator += (Vector_V_sum_S_mul_M_mul_V<VA,SCALAR,MATRIX,VB> const & op) { \
      *this += op.a; \
      op.M.add_S_mul_M_mul_V(*this, op.s, op.b); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c  a+s*(M*b)  contained in \c op of type \c Vector_V_sum_S_mul_M_mul_V */ \
    template <typename VA, typename SCALAR, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator -= (Vector_V_sum_S_mul_M_mul_V<VA,SCALAR,MATRIX,VB> const & op) { \
      *this -= op.a; \
      op.M.add_S_mul_M_mul_V(*this, -op.s, op.b); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the expression \c  a+s*(M^b)  contained in \c op of type \c Vector_V_sum_S_mul_M_mul_V */ \
    template <typename VA, typename SCALAR, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator = (Vector_V_sum_S_mul_Mt_mul_V<VA,SCALAR,MATRIX,VB> const & op) { \
      *this = op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, op.s, op.b); \
      return *this; \
    } \
    /*! \brief  Add to \c *this the evaluation of the expression \c  a+s*(M^b)  contained in \c op of type \c Vector_V_sum_S_mul_M_mul_V */ \
    template <typename VA, typename SCALAR, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator += (Vector_V_sum_S_mul_Mt_mul_V<VA,SCALAR,MATRIX,VB> const & op) { \
      *this += op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, op.s, op.b); \
      return *this; \
    } \
    /*! \brief  Subtract to \c *this the evaluation of the expression \c  a+s*(M^b)  contained in \c op of type \c Vector_V_sum_S_mul_M_mul_V */ \
    template <typename VA, typename SCALAR, typename MATRIX, typename VB> inline \
    VECTOR const & \
    operator -= (Vector_V_sum_S_mul_Mt_mul_V<VA,SCALAR,MATRIX,VB> const & op) { \
      *this -= op.a; \
      op.M.add_S_mul_Mt_mul_V(*this, -op.s, op.b); \
      return *this; \
    } \
    /*! \brief  Assign to \c *this the evaluation of the preconditioner \c  a/P contained in \c op of type \c Vector_V_div_P */ \
    template <typename VEC, typename PRECO> inline \
    VECTOR const & \
    operator = (Vector_V_div_P<VEC,PRECO> const & op) { \
      op.P.assPreco(*this, op.a); \
      return *this; \
    }

  /*
  //  #     #                                   ######                       
  //  #     # ######  ####  #####  ####  #####  #     #   ##    ####  ###### 
  //  #     # #      #    #   #   #    # #    # #     #  #  #  #      #      
  //  #     # #####  #        #   #    # #    # ######  #    #  ####  #####  
  //   #   #  #      #        #   #    # #####  #     # ######      # #      
  //    # #   #      #    #   #   #    # #   #  #     # #    # #    # #      
  //     #    ######  ####    #    ####  #    # ######  #    #  ####  ###### 
  */

  //! Base vector class
  /*!
   * This class in the base class for all the
   * vector classes of \c SparseTool.
   * This class is incomplete and is used as a pivot for internal operations.
   */
  template <typename T, typename VECTOR>
  class VectorBase {
  public:
    typedef T                    valueType; //!< type of the element of the vector
    typedef VectorBase<T,VECTOR> VectorT;   //!< the vector type

  public:

    //! build and empty vector 
    VectorBase() {};

    //! forward the \c size() operator to the derived class 
    indexType size() const
    { return static_cast<VECTOR const *>(this) -> size(); }

    //! forward the access to the \a i-th element of the vector 
    valueType const &
    operator [] ( indexType i ) const
    { return static_cast<VECTOR const *>(this) -> operator [] (i); }

    //! forward the access to the \a i-th element of the vector 
    valueType &
    operator [] ( indexType i )
    { return static_cast<VECTOR *>(this) -> operator [] (i); }

    //! forward the access to the \a i-th element of the vector 
    valueType const &
    operator () ( indexType i ) const
    { return static_cast<VECTOR const *>(this) -> operator () (i); }

    //! forward the access to the \a i-th element of the vector 
    valueType &
    operator () ( indexType i )
    { return static_cast<VECTOR *>(this) -> operator () (i); }

    SPARSELIB_VECTOR_OPERATIONS(VectorT)

  };

  /*
  //  #     #                                   
  //  #     # ######  ####  #####  ####  #####  
  //  #     # #      #    #   #   #    # #    # 
  //  #     # #####  #        #   #    # #    # 
  //   #   #  #      #        #   #    # #####  
  //    # #   #      #    #   #   #    # #   #  
  //     #    ######  ####    #    ####  #    # 
  */

  //! Variable size full vector class
  /*!
    This class extend the STL vector class by adding some math operation
    and interaction with sparse matrix classes.
    
    \par Usage

    A \c Vector\<T\> is defined by specifying the type \c T  and
    optionally its size.  For example,

\code
  Vector<double> b, c(100);
\endcode

    defines \c b as a vector of \b double of size \b 0 while \c c is a
    vector of \b double of size \b 100. You can change the size of the
    vector by the methods \c resize as follows:

\code
  Vector<double> d;
  d.resize(200);
\endcode

    so that \c d is a \c Vector<double> of size \b 200.
    There are many methods associate to a \c Vector\<T\> in the following
    paragraph they are listed.

    \par Constructors
\code
1  Vector<T> v;
2  Vector<T> v(dim);
3  Vector<T> w(v);
\endcode

    On line \b 1 construct the \c Vector\<T\> \c v of size
    \b 0.  On line \b 2 construct the \c Vector\<T\> \c v of of size \c dim.
    On line \b 3 construct the \c Vector\<T\> \c w as a copy of
    \c Vector\<T\> \c v.

    \par Indexing
    Vector instances are indexed as one-dimensional C arrays, and the
    index numbering follows the standard C convention, starting from zero.

    Let us define \c v as \c Vector\<T\>

\code
  Vector<T> v;
\endcode

    Then, \c v[i] returns a reference to the \c T&-type
    \c i-th element of \c v;

    \par Changing Dimension
    It is possible to change the size of a \c Vector\<T\> object.
    For example

\code
  Vector<double> d;
  d.resize(200);
\endcode

   the \c Vector\<T\> \c d has size \b 200.
   The method \c size() return the actual size of the \c Vector\<T\>.
   For example defining

\code
  Vector<double> d(123);
\endcode

   the method \c d.size() return \b 123.

   \par Initialization

    It is possible to initialize all the components of a \c Vector\<T\>
    to a value, for example

\code
  Vector<float> v(100);
  v = 3.14;
\endcode

    is equivalent to

\code
  Vector<float> v(100);
  for ( int i = 0; i < v.size(); ++i ) v[i] = 3.14;
\endcode

    although is done more efficiently by the library.

    \par Assignment
    
    It is possible to copy the contents of a \c Vector\<T\> to another
    one as the following example show

\code
1:  Vector<double> a, b, c;

2:  a.resize(100);
3:  b.resize(200);
4:  c.resize(150);
  
5:  c = 3;
6:  b = c;
7:  a = c;
\endcode

    the meaning of line \b 5 should be clear.  Lines \b 6 and
    \b 7 are equivalent to

\code
  int i;
  for ( i = 0; i < min(b.size(), c.size()); ++i ) b[i] = c[i];
  for ( i = 0; i < min(a.size(), c.size()); ++i ) a[i] = c[i];
\endcode

    where you can notice that only the values that can be stored are
    assigned.

    It is possible to initialize many vectors same value as in the
    following expressions

\code
1:  Vector<double> a, b, c;   
2:  a.resize(100);
3:  b.resize(200);
4:  c.resize(150);
5:  a = b = c = 3;
\endcode

    but take attention because line \b 5 is not equivalent to

\code
  a = 3;
  b = 3;
  c = 3;
\endcode

in fact line \b 5 is equivalent to

\code
  c = 3;
  b = c;
  a = b;
\endcode

    so that we have

    - the \c Vector\<T\> \c c is initialized with \a all its
      \b 150 elements set to \b 5.
    - the \c Vector\<T\> \c b is initialized with \a only its
      first \b 150 elements set to \b 1 while the remaining
       are undefined
    - the \c Vector\<T\> \c a is initialized with \a all its
      first \b 100 elements set to \b 1.
 
    \par Arithmetic Operators on \c Vector\<T\>

    A set of usual arithmetic operators are explicitly defined on
    vector-type data.  If not otherwise specified, the operators extend
    the corresponding scalar operation in a \a component-wise fashion. 
    Hence, for vectors with size \c dim, the component index \c i
    in all the following expressions is supposed to run through \b 0
    to \c dim-1.

    Let us define the three double precision vectors \c a, \c b,
    and \c c, that we shall use in all the following examples

\code
  int const dim = 100;
  Vector<double> a(dim), b(dim), c(dim);
\endcode

    The arithmetic operators defined on vectors are given in the following
    sections.

    \par \c Scalar-Vector internal operations

\verbatim
Command    Equivalence
a += 2     for ( i=0; i < a.size(); ++i ) a[i] += 2;
a -= 2     for ( i=0; i < a.size(); ++i ) a[i] -= 2;
a *= 2     for ( i=0; i < a.size(); ++i ) a[i] *= 2;
a /= 2     for ( i=0; i < a.size(); ++i ) a[i] /= 2;
\endverbatim

    \par \c Scalar-Vector operations

\verbatim
Command      Equivalence
             sz = min( a.size(), b.size() );
a = b + 2    for ( i=0; i < sz; ++i ) a[i] = b[i] + 2;
a = 3 + b    for ( i=0; i < sz; ++i ) a[i] = 3 + b[i];
a = b - 2    for ( i=0; i < sz; ++i ) a[i] = b[i] - 2;
a = 3 - b    for ( i=0; i < sz; ++i ) a[i] = 3 - b[i];
a = b * 2    for ( i=0; i < sz; ++i ) a[i] = b[i] * 2;
a = 3 * b    for ( i=0; i < sz; ++i ) a[i] = 3 * b[i];
a = b / 2    for ( i=0; i < sz; ++i ) a[i] = b[i] / 2;
a = 3 / b    for ( i=0; i < sz; ++i ) a[i] = 3 / b[i];
\endverbatim

  \par \c Vector-Vector internal operations

\verbatim
Command      Equivalence
             sz = min( a.size(), b.size() );
a += b       for ( i=0; i < sz; ++i ) a[i] += b[i];
a -= b       for ( i=0; i < sz; ++i ) a[i] -= b[i];
a *= b       for ( i=0; i < sz; ++i ) a[i] *= b[i];
a /= b       for ( i=0; i < sz; ++i ) a[i] /= b[i];
\endverbatim

    \par \c Vector-Vector

\verbatim
Command      Equivalence
             sz = min( a.size(), b.size() );
b = +a       for ( i=0; i < sz; ++i ) b[i] = +a[i];
b = -a       for ( i=0; i < sz; ++i ) b[i] = -a[i];

             sz = min( sz, c.size() );
c = a + b    for ( i=0; i < sz; ++i ) c[i] = a[i] + b[i];
c = a - b    for ( i=0; i < sz; ++i ) c[i] = a[i] - b[i];
c = a * b    for ( i=0; i < sz; ++i ) c[i] = a[i] * b[i];
c = a / b    for ( i=0; i < sz; ++i ) c[i] = a[i] / b[i];
\endverbatim

    \par Function of \c Vector\<T\>
  
    Let be <c> n = min( a.size(), b.size() ) </c>,

    - \c dot(a,b)   \f$ = \displaystyle\sum_{i=0}^{n-1} \overline{a_i} b_i \f$
    - \c rdot(a,b)  \f$ = \displaystyle\sum_{i=0}^{n-1} a_i b_i \f$
    - \c dist(a,b)  \f$ = \sqrt{\sum_{i=0}^{n-1} \left|a_i - b_i\right|^2}\f$
    - \c normi(a)   \f$ = ||a||_{\infty} = \max\left\{|a_0|,|a_1|,\ldots,|a_{n-1}|\right\}\f$
    - \c norm1(a)   \f$ = ||a||_{1} = \sum_{i=0}^{n-1} |a_i|\f$
    - \c norm2(a)   \f$ = ||a||_{2} = \sqrt{ \sum_{i=0}^{n-1} a_i^2}\f$
    - \c normp(a,p) \f$ = ||a||_{p} = \left( \sum_{i=0}^{n-1} |a_i|^p \right)^{1/p} \f$
    - \c sum(a)     \f$ = \displaystyle\sum_{i=0}^{n-1} a_i \f$
    - \c prod(a)    \f$ = \displaystyle\prod_{i=0}^{n-1} a_i \f$
    - \c max(a)     \f$ = \max\left\{|a_0|,|a_1|,\ldots,|a_{n-1}|\right\} \f$
    - \c min(a)     \f$ = \min\left\{|a_0|,|a_1|,\ldots,|a_{n-1}|\right\} \f$
  */
  template <typename T>
  class Vector : public VectorBase<T,Vector<T> >,
                 public vector<T> {
  public:

    typedef T valueType;  //!< assign to \c valueType the type of the vector

    //! Build an empty Vector with 0 elements
    Vector() : vector<T>(0) {};

    //! Build a Vector with \a numElement elements
    Vector( indexType numElement ) : vector<T>(numElement) {};
    
    //! Return the total number of elements of the vector
    indexType size() const { return indexType(vector<T>::size()); }

    /*! \brief  Access to the \a i-th element of the Vector.
        If \c SPARSETOOL_DEBUG is defined an index bound check is performed.  
     */
    valueType const &
    operator [] ( indexType i ) const {
      SPARSETOOL_TEST( i < vector<T>::size(),"Vector[" << i << "] Vector size = " << vector<T>::size() << " bad index")
      return vector<T>::operator [] (i);
    }

    /*! \brief  Access to the \a i-th element of the Vector.
        If \c SPARSETOOL_DEBUG is defined an index bound check is performed.  
     */
    valueType & 
    operator [] ( indexType i ) {
      SPARSETOOL_TEST( i < indexType(vector<T>::size()), "Vector[" << i << "] Vector size = " << vector<T>::size() << " bad index")
      return vector<T>::operator [] (i);
    }

    /*!
     * \brief  Access to the \a i-th element of the Vector.
     * No index bound check is performed.
     */
    valueType const &
    operator () ( indexType i ) const
    { return vector<T>::operator [] (i); }

    /*!
     * \brief  Access to the \a i-th element of the Vector.
     * No index bound check is performed.
     */
    valueType &
    operator () ( indexType i )
    { return vector<T>::operator [] (i); }

    //! Get a copy of the Vector \a v.
    template <typename VECTOR> inline
    void
    load( VectorBase<T,VECTOR> const & v ) {
      this -> resize( v.size() );
      SPARSELIB_LOOP( indexType(vector<T>::size()), (*this)(i) = v(i) );
    }

    //! Set all elements to 0
    void
    setZero() {
      SPARSELIB_LOOP( indexType(vector<T>::size()), (*this)(i) = valueType(0) );
    }

    //! Set all elements to a
    void
    setTo( valueType const & a ) {
      SPARSELIB_LOOP( indexType(vector<T>::size()), (*this)(i) = a );
    }

    //! Fill the Vector with the value \a v from index \a b to \a e.  
    void
    fill( indexType b, indexType e, valueType const & v ) {
      SPARSETOOL_TEST(
        b <= e,
        "Vector::fill(" << b << "," << e << ") bad range"
      )
      SPARSETOOL_TEST(
        e < indexType(vector<T>::size()),
        "Vector::fill(" << b << "," << e << ") out of range"
      )
      valueType * p = vector<T>::begin() + b;
      SPARSELIB_LOOP( e-b, *p++ = v);
    }
    
    SPARSELIB_VECTOR_OPERATIONS(Vector<T>)

  };
  
  /*
  //  #     #                                       #####                         
  //  #     # ######  ####  #####  ####  #####     #     # #      #  ####  ###### 
  //  #     # #      #    #   #   #    # #    #    #       #      # #    # #      
  //  #     # #####  #        #   #    # #    #     #####  #      # #      #####  
  //   #   #  #      #        #   #    # #####           # #      # #      #      
  //    # #   #      #    #   #   #    # #   #     #     # #      # #    # #      
  //     #    ######  ####    #    ####  #    #     #####  ###### #  ####  ###### 
  */
  //! Remapping a piece of memory to a vector
  /*!
    This class extend the \c VectorBase class by adding the capacity
    of remapping a piece o \c Vector or a C-pointer
    
    \par Usage

    A \c VectorSlice\<T\> is defined by specifying the type \c T.  For example,

\code
  VectorSlice<double> b, c;
\endcode

    There are many methods associate to a \c Vector\<T\> in the following
    paragraph they are listed.

    \par Slicing
\code
1  Vector<double> v(100); double w[100];
2  b.slice( v, 10, 20 );
3  c.slice( w + 5, w + 45 );
\endcode

    On line \b 1 construct the \c Vector<double> \c v of size \b 100
    and the C array \c w of \b 100 elements.
    On line \b 2 remap the components of the vecotor \c v from \c 10 to \c 19 to the "vector" b.
    On line \b 3 remap the components of the Array \c w from \c 5 to \c 44 to the "vector" c.
    \c Vector\<T\> \c v.

    \par Indexing
    Vector instances are indexed as one-dimensional C arrays, and the
    index numbering follows the standard C convention, starting from zero.

    Then, \c b[i] returns a reference to the \c T&-type
    \c i-th element of \c v and due to slicing the \c b[i]==v[i+10]. 
    Analogously \c c[i]==w[i+5]. 
  */
  template <typename T>
  class VectorSlice : public VectorBase<T,VectorSlice<T> > {
  public:
    typedef T   valueType; //!< type of the element of the vector
    indexType   len;       //!< length of the slice
    valueType * values;    //!< poiter to the beginning of the slice

  public:

    //! Build an empty Vector with 0 elements
    VectorSlice() : len(0), values(nullptr) {};

    /*!
     *  \brief  Map the vector \c v from \c v(begin) to \c v(end-1) in
     *  (*this)(0) to (*this)(end-begin)
     */
    VectorSlice<T> &
    slice( Vector<T> & v, indexType begin, indexType end ) {
      values = &v(begin);
      len    = end - begin;
      return *this;
    };

    VectorSlice<T> const &
    slice( Vector<T> const & v, indexType begin, indexType end ) {
      values = &v(begin);
      len    = end - begin;
      return *this;
    };

    /*!
     *  \brief  Map the piece of memory from \c pBegin to \c pEnd-1 in
     *  (*this)(0) to (*this)(pEnd-pBegin)
     */
    VectorSlice<T> &
    slice( T * pBegin, T * pEnd )
    { values = pBegin; len = pEnd - pBegin; return *this; };

    //! return the total number of element of the vector slice
    indexType size() const { return len; }

    /*!
     *  \brief  Access to the \a i-th element of the Vector.
     *  If \c SPARSETOOL_DEBUG is defined an index bound check is performed.
     */
    valueType const &
    operator [] ( indexType i ) const {
      SPARSETOOL_TEST(
        i < len,
        "VectorSlice[" << i << "] size = " << len << " bad index"
      )
      return values[i];
    }

    /*!
     *  \brief  Access to the \a i-th element of the Vector.
     *  If \c SPARSETOOL_DEBUG is defined an index bound check is performed.
     */
    valueType & 
    operator [] ( indexType i ) {
      SPARSETOOL_TEST(
        i < len,
        "VectorSlice[" << i << "] size = " << len << " bad index"
      )
      return values[i];
    }

    /*!
     *  \brief  Access to the \a i-th element of the Vector.
     *  No index bound check is performed.
     */
    valueType const &
    operator () ( indexType i ) const
    { return values[i]; }

    /*!
     *  \brief  Access to the \a i-th element of the Vector.
     *  No index bound check is performed.
     */
    valueType & 
    operator () ( indexType i )
    { return values[i]; }

    SPARSELIB_VECTOR_OPERATIONS(VectorSlice<T>)

  };

  /*! \cond NODOC */

  #define SPARSELIB_V(OP,OP_V)                                                \
  template<typename A>                                                        \
  class OP_V {                                                                \
  public:                                                                     \
    SPARSELIB_TYPES_FROM_TYPENAME(A);                                        \
  private:                                                                    \
    A const & a;                                                             \
  public:                                                                     \
    OP_V(A const & aa) : a(aa) { }                                            \
    valueType operator () (indexType i) const { return OP a(i); }            \
    indexType size(void) const { return a.size(); }                          \
  };                                                                          \
                                                                              \
  template<typename T, typename VEC> inline                                   \
  VectorE<T,OP_V<VectorBase<T,VEC> > >                                        \
  operator OP (VectorBase<T,VEC> const & a) {                                 \
    typedef OP_V<VectorBase<T,VEC> > op;                                     \
    return VectorE<T,op>(op(a));                                              \
  }                                                                           \
                                                                              \
  template <typename T, typename A> inline                                    \
  VectorE<T,OP_V<VectorE<T,A> > >                                             \
  operator OP (VectorE<T,A> const & a) {                                      \
    typedef OP_V<VectorE<T,A> > op;                                          \
    return VectorE<T,op>(op(a));                                              \
  }

  #define SPARSELIB_VV(OP,V_OP_V)                                             \
  template<typename A, typename B>                                            \
  class V_OP_V {                                                              \
  public:                                                                     \
    SPARSELIB_TYPES_FROM_TYPENAME2(A,B);                                     \
  private:                                                                    \
    A const & a;                                                             \
    B const & b;                                                             \
  public:                                                                     \
    V_OP_V(A const & aa, B const & bb) : a(aa), b(bb) { }                     \
    valueType operator () (indexType i) const { return a(i) OP b(i); }       \
    indexType size(void) const { return minIndex(a.size(), b.size()); }      \
  };                                                                         \
                                                                              \
  template<typename T, typename VA, typename VB> inline                       \
  VectorE<T,V_OP_V<VectorBase<T,VA>,VectorBase<T,VB> > >                      \
  operator OP (VectorBase<T,VA> const & a, VectorBase<T,VB> const & b) {      \
    typedef V_OP_V<VectorBase<T,VA>,VectorBase<T,VB> > op;                   \
    return VectorE<T,op>(op(a,b));                                           \
  }                                                                           \
                                                                              \
  template <typename T, typename A, typename VB> inline                       \
  VectorE<T,V_OP_V<VectorE<T,A>,VectorBase<T,VB> > >                          \
  operator OP (VectorE<T,A> const & a, VectorBase<T,VB> const & b) {          \
    typedef V_OP_V<VectorE<T,A>,VectorBase<T,VB> > op;                       \
    return VectorE<T,op>(op(a,b));                                            \
  }                                                                           \
                                                                              \
  template <typename T, typename VA, typename B> inline                       \
  VectorE<T,V_OP_V<VectorBase<T,VA>,VectorE<T,B> > >                          \
  operator OP (VectorBase<T,VA> const & a, VectorE<T,B> const & b) {          \
    typedef V_OP_V<VectorBase<T,VA>,VectorE<T,B> > op;                       \
    return VectorE<T,op>(op(a,b));                                            \
  }                                                                           \
                                                                              \
  template <typename T, typename A, typename B> inline                        \
  VectorE<T,V_OP_V<VectorE<T,A>, VectorE<T,B> > >                             \
  operator OP (VectorE<T,A> const & a, VectorE<T,B> const & b) {              \
    typedef V_OP_V<VectorE<T,A>,VectorE<T,B> > op;                           \
    return VectorE<T,op>(op(a,b));                                            \
  }

  #define SPARSELIB_SV(OP,S_OP_V)                                             \
  template<typename S, typename B>                                            \
  class S_OP_V {                                                              \
  public:                                                                     \
    SPARSELIB_TYPES_FROM_TYPENAME3(S,B);                                     \
  private:                                                                    \
    S const & a;                                                             \
    B const & b;                                                             \
  public:                                                                     \
    S_OP_V(S const & aa, B const & bb) : a(aa), b(bb) { }                     \
    valueType operator () (indexType i) const { return a OP b(i); }          \
    indexType size(void) const { return b.size(); }                          \
  };                                                                         \
                                                                              \
  template <typename T, typename VEC> inline                                  \
  VectorE<T,S_OP_V<T,VectorBase<T,VEC> > >                                    \
  operator OP (T const & s, VectorBase<T,VEC> const & v) {                    \
    typedef S_OP_V<T,VectorBase<T,VEC> > op;                                 \
    return VectorE<T,op>(op(s,v));                                            \
  }                                                                           \
                                                                              \
  template <typename T, typename A> inline                                    \
  VectorE<T,S_OP_V<T,VectorE<T,A> > >                                         \
  operator OP (T const & s, VectorE<T,A> const & v) {                         \
    typedef S_OP_V<T,VectorE<T,A> > op;                                      \
    return VectorE<T,op>(op(s,v));                                            \
  }

  #define SPARSELIB_VS(OP,V_OP_S)                                             \
  template<typename S, typename A>                                            \
  class V_OP_S {                                                              \
  public:                                                                     \
    SPARSELIB_TYPES_FROM_TYPENAME3(S,A);                                     \
  private:                                                                    \
    A const & a;                                                             \
    S const & b;                                                             \
  public:                                                                     \
    V_OP_S(A const & aa, S const & bb) : a(aa), b(bb) { }                     \
    valueType operator () (indexType i) const { return a(i) OP b; }          \
    indexType size(void) const { return a.size(); }                          \
  };                                                                         \
                                                                              \
  template <typename T, typename VEC> inline                                  \
  VectorE<T,V_OP_S<T,VectorBase<T,VEC> > >                                    \
  operator OP (VectorBase<T,VEC> const & v, T const & s) {                    \
    typedef V_OP_S<T,VectorBase<T,VEC> > op;                                 \
    return VectorE<T,op>(op(v,s));                                           \
  }                                                                           \
                                                                              \
  template <typename T, typename A> inline                                    \
  VectorE<T, V_OP_S<T,VectorE<T,A> > >                                        \
  operator OP (VectorE<T,A> const & v, T const & s) {                         \
    typedef V_OP_S<T,VectorE<T,A> > op;                                      \
    return VectorE<T,op>(op(v,s));                                           \
  }
  
  /*! \endcond */

  /*! \cond NODOC */

  SPARSELIB_V(-,Vector_neg_V)

  SPARSELIB_VV(+,Vector_V_sum_V)
  SPARSELIB_VV(-,Vector_V_sub_V)
  SPARSELIB_VV(*,Vector_V_mul_V)
  SPARSELIB_VV(/,Vector_V_div_V)

  SPARSELIB_SV(+,Vector_S_sum_V)
  SPARSELIB_SV(-,Vector_S_sub_V)
  SPARSELIB_SV(*,Vector_S_mul_V)
  SPARSELIB_SV(/,Vector_S_div_V)

  SPARSELIB_VS(+,Vector_V_sum_S)
  SPARSELIB_VS(-,Vector_V_sub_S)
  SPARSELIB_VS(*,Vector_V_mul_S)
  SPARSELIB_VS(/,Vector_V_div_S)

  /*! \endcond */
  
  /*! \cond NODOC */  

  #define SPARSELIB_F_V(FUN)                                                  \
  template<typename A>                                                        \
  class Vector_##FUN {                                                        \
  public:                                                                     \
    SPARSELIB_TYPES_FROM_TYPENAME_FUN1(A);                                   \
  private:                                                                    \
    A const & a;                                                             \
  public:                                                                     \
    Vector_##FUN(A const & aa) : a(aa) { }                                    \
    valueType operator () (indexType i) const                                 \
      { return ::SparseToolFun::FUN(a(i)); }                                 \
    indexType size(void) const { return a.size(); }                          \
  };

  #define SPARSELIB_F_V_TREE(FUN)                                             \
  template<typename T> inline                                                 \
  VectorE<T,Vector_##FUN<Vector<T> > >                                        \
  FUN(Vector<T> const & a) {                                                  \
    typedef Vector_##FUN<Vector<T> > op;                                     \
    return VectorE<T,op>(op(a));                                              \
  }                                                                           \
                                                                              \
  template <typename T, typename A> inline                                    \
  VectorE<T,Vector_##FUN<VectorE<T,A> > >                                     \
  FUN(VectorE<T,A> const & a) {                                               \
    typedef Vector_##FUN<VectorE<T,A> > op;                                  \
    return VectorE<T,op>(op(a));                                              \
  }

  #define SPARSELIB_F_VV(FUN)                                                 \
  template<typename A, typename B>                                            \
  class Vector_##FUN {                                                        \
  public:                                                                     \
    SPARSELIB_TYPES_FROM_TYPENAME_FUN2(A,B);                                 \
  private:                                                                    \
    A const & a;                                                             \
    B const & b;                                                             \
  public:                                                                     \
    Vector_##FUN(A const & aa, B const & bb) : a(aa), b(bb) { }               \
    valueType operator () (indexType i) const                                 \
      { return ::SparseToolFun::FUN(a(i), b(i)); }                           \
    indexType size(void) const { return minIndex(a.size(), b.size()); }      \
  };

  #define SPARSELIB_F_VV_TREE(FUN)                                            \
  template<typename T> inline                                                 \
  VectorE<T,Vector_##FUN<Vector<T>,Vector<T> > >                              \
  FUN(Vector<T> const & a, Vector<T> const & b) {                             \
    typedef Vector_##FUN<Vector<T>,Vector<T> > op;                           \
    return VectorE<T,op>(op(a,b));                                           \
  }                                                                           \
                                                                              \
  template <typename T, typename A> inline                                    \
  VectorE<T,Vector_##FUN<VectorE<T,A>,Vector<T> > >                           \
  FUN(VectorE<T,A> const & a, Vector<T> const & b) {                          \
    typedef Vector_##FUN<VectorE<T,A>,Vector<T> > op;                        \
    return VectorE<T,op>(op(a,b));                                            \
  }                                                                           \
                                                                              \
  template <typename T, typename B> inline                                    \
  VectorE<T,Vector_##FUN<Vector<T>,VectorE<T,B> > >                           \
  FUN(Vector<T> const & a, VectorE<T,B> const & b) {                          \
    typedef Vector_##FUN<Vector<T>,VectorE<T,B> > op;                        \
    return VectorE<T,op>(op(a,b));                                            \
  }                                                                           \
                                                                              \
  template <typename T, typename A, typename B> inline                        \
  VectorE<T,Vector_##FUN<VectorE<T,A>, VectorE<T,B> > >                       \
  FUN(VectorE<T,A> const & a, VectorE<T,B> const & b) {                       \
    typedef Vector_##FUN<VectorE<T,A>,VectorE<T,B> > op;                     \
    return VectorE<T,op>(op(a,b));                                            \
  }
  /*! \endcond */

  /*! \cond NODOC */

  SPARSELIB_F_V(absval)
  SPARSELIB_F_V(sin)
  SPARSELIB_F_V(cos)
  SPARSELIB_F_V(tan)
  SPARSELIB_F_V(asin)
  SPARSELIB_F_V(acos)
  SPARSELIB_F_V(atan)
  SPARSELIB_F_V(cosh)
  SPARSELIB_F_V(sinh)
  SPARSELIB_F_V(tanh)

  SPARSELIB_F_V(sqrt)
  SPARSELIB_F_V(ceil)
  SPARSELIB_F_V(floor)
  SPARSELIB_F_V(log)
  SPARSELIB_F_V(log10)

  SPARSELIB_F_VV(pow)
  SPARSELIB_F_VV(atan2)

  SPARSELIB_F_VV(minval)
  SPARSELIB_F_VV(maxval)

  SPARSELIB_F_VV(minabsval)
  SPARSELIB_F_VV(maxabsval)

  /*! \endcond */

  /*! \cond NODOC */

  template<typename A>
  class F1 {
  public:
    typedef typename A::valueType valueType;

    static inline
    typename return_trait<valueType>::valueType
    norm1(A const & a) {
      using ::SparseToolFun::absval;
      valueType res = 0;
      SPARSELIB_LOOP( a.size(), res += absval(a(i)) );
      return res;
    }

    static inline
    typename return_trait<valueType>::valueType
    normi(A const & a) {
      using ::SparseToolFun::absval;
      typename return_trait<valueType>::valueType bf, res = 0;
      SPARSELIB_LOOP( a.size(), bf = absval(a(i)); if ( bf > res ) res = bf );
      return res;
    }

    static inline
    typename return_trait<valueType>::valueType
    norm2(A const & a) {
      using ::SparseToolFun::absval2;
      using ::SparseToolFun::sqrt;
      typename return_trait<valueType>::valueType res = 0;
      SPARSELIB_LOOP( a.size(), res += absval2(a(i)) );
      return sqrt(res);
    }

    static inline
    typename return_trait<valueType>::valueType
    normp(A const & a, valueType const & p) {
      using ::SparseToolFun::absval;
      using ::SparseToolFun::pow;
      typename return_trait<valueType>::valueType res = 0;
      SPARSELIB_LOOP( a.size(), res += pow(absval(a(i)),p) );
      return pow(res,1/p);
    }

    static inline
    valueType
    sum(A const & a) {
      valueType res = 0;
      SPARSELIB_LOOP( a.size(), res += a(i) );
      return res;
    }

    static inline
    valueType
    prod(A const & a) {
      valueType res = 1;
      SPARSELIB_LOOP( a.size(), res *= a(i) );
      return res;
    }

    static inline
    valueType
    maxval(A const & a) {
      valueType res = a(0);
      SPARSELIB_LOOP( a.size(), if ( a(i) > res ) res = a(i) );
      return res;
    }

    static inline
    valueType
    minval(A const & a) {
      valueType res = a(0);
      SPARSELIB_LOOP( a.size(), if ( a(i) < res ) res = a(i) );
      return res;
    }

    static inline
    typename return_trait<valueType>::valueType
    maxabsval(A const & a) {
      using ::SparseToolFun::absval;
      typename return_trait<valueType>::valueType tmp, res = absval(a(0));
      SPARSELIB_LOOP( a.size(), tmp = absval(a(i)); if ( tmp > res ) res = tmp );
      return res;
    }

    static inline
    typename return_trait<valueType>::valueType
    minabsval(A const & a) {
      using ::SparseToolFun::absval;
      typename return_trait<valueType>::valueType tmp, res = absval(a(0));
      SPARSELIB_LOOP( a.size(), tmp = absval(a(i)); if ( tmp < res ) res = tmp );
      return res;
    }

  };

  template<typename A, typename B>
  class F2 {
  public:

    SPARSELIB_TYPES_FROM_TYPENAME2(A,B);

    static inline
    valueType
    dot(A const & a, B const & b) {
      using ::SparseToolFun::conj;
      valueType res = 0;
      SPARSELIB_LOOP( minIndex(a.size(),b.size()), res += conj(a(i)) * b(i) );
      return res;
    }

    static inline
    valueType
    rdot(A const & a, B const & b) {
      using ::SparseToolFun::conj;
      valueType res = 0;
      SPARSELIB_LOOP( minIndex(a.size(),b.size()), res += a(i) * b(i) );
      return res;
    }

    static inline
    typename return_trait<valueType>::valueType
    dist2(A const & a, B const & b) {
      using ::SparseToolFun::absval2;
      typename return_trait<valueType>::valueType res = 0;
      SPARSELIB_LOOP( minIndex( a.size(), b.size() ), res += absval2(a(i) - b(i)) );
      return res;
    }

    static inline
    typename return_trait<valueType>::valueType
    dist(A const & a, B const & b) {
      using ::SparseToolFun::absval2;
      typename return_trait<valueType>::valueType res = 0;
      SPARSELIB_LOOP( minIndex( a.size(), b.size() ), res += absval2(a(i) - b(i)) );
      return ::SparseToolFun::sqrt(res);
    }

  };
  
  /*! \endcond */

  /*! \cond NODOC */

  SPARSELIB_F_V_TREE(absval)
  SPARSELIB_F_V_TREE(sin)
  SPARSELIB_F_V_TREE(cos)
  SPARSELIB_F_V_TREE(tan)
  SPARSELIB_F_V_TREE(asin)
  SPARSELIB_F_V_TREE(acos)
  SPARSELIB_F_V_TREE(atan)
  SPARSELIB_F_V_TREE(cosh)
  SPARSELIB_F_V_TREE(sinh)
  SPARSELIB_F_V_TREE(tanh)

  SPARSELIB_F_V_TREE(sqrt)
  SPARSELIB_F_V_TREE(ceil)
  SPARSELIB_F_V_TREE(floor)
  SPARSELIB_F_V_TREE(log)
  SPARSELIB_F_V_TREE(log10)

  SPARSELIB_F_VV_TREE(pow)
  SPARSELIB_F_VV_TREE(atan2)
  SPARSELIB_F_VV_TREE(minval)
  SPARSELIB_F_VV_TREE(maxval)
  SPARSELIB_F_VV_TREE(minabsval)
  SPARSELIB_F_VV_TREE(maxabsval)

  /*! \endcond */

  #undef SPARSELIB_V
  #undef SPARSELIB_VV
  #undef SPARSELIB_SV
  #undef SPARSELIB_VS
  #undef SPARSELIB_F_V
  #undef SPARSELIB_F_V_TREE
  #undef SPARSELIB_F_VV
  #undef SPARSELIB_F_VV_TREE
  
  //! \name Vector norm
  //@{

  // N O R M 1
  //! Evaluate the 1-norm of the vector \c v
  template <typename T, typename V> inline
  typename return_trait<T>::valueType norm1(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::norm1(v); }

  //! Evaluate the 1-norm of the vector expression \c e
  template <typename T, typename A> inline
  typename return_trait<T>::valueType norm1(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::norm1(e); }

  // N O R M 2
  //! Evaluate the 2-norm of the vector \c v 
  template <typename T, typename V> inline
  typename return_trait<T>::valueType norm2(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::norm2(v); }

  //! Evaluate the 2-norm of the vector expression \c e
  template <typename T, typename A> inline
  typename return_trait<T>::valueType norm2(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::norm2(e); }

  // N O R M I
  //! Evaluate the infinity-norm of the vector \c v
  template <typename T, typename V> inline
  typename return_trait<T>::valueType normi(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::normi(v); }

  //! Evaluate the infinity-norm of the vector expression \c e
  template <typename T, typename A> inline
  typename return_trait<T>::valueType normi(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::normi(e); }

  // N O R M P
  //! Evaluate the p-norm of the vector \c v
  template <typename T, typename S, typename V> inline
  typename return_trait<T>::valueType normp(VectorBase<T,V> const & v, S const & p)
  { return F1<VectorBase<T,V> >::normp(v,T(p)); }

  //! Evaluate the p-norm of the vector \c e
  template <typename T, typename A, typename S> inline
  typename return_trait<T>::valueType normp(VectorE<T,A> const & e, S const & p)
  { return F1<VectorE<T,A> >::normp(e,T(p)); }
  //@}

  //! \name sum and prod
  //@{

  // S U M
  //! Evaluate sum of the entries the vector \c v
  template <typename T, typename V> inline
  T sum(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::sum(v); }

  //! Evaluate sum of the entries the vector expression \c v
  template <typename T, typename A> inline
  T sum(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::sum(e); }
  // S U M

  //! Evaluate product of the entries the vector \c v
  template <typename T, typename V> inline
  T prod(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::prod(v); }

  //! Evaluate product of the entries the vector expression \c v
  template <typename T, typename A> inline
  T prod(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::prod(e); }

  //@}

  //! \name Vector minimum maximum values
  //@{

  // M A X
  //! Evaluate the maximum value of the vector \c v
  template <typename T, typename V> inline
  T maxval(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::maxval(v); }

  //! Evaluate the maximum value of the vector expression \c e
  template <typename T, typename A> inline
  T maxval(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::maxval(e); }

  // M I N
  //! Evaluate the minimum value of the vector \c v
  template <typename T, typename V> inline
  T minval(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::minval(v); }

  //! Evaluate the maximum value of the vector \c e
  template <typename T, typename A> inline
  T minval(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::minval(e); }

  // M A X
  //! Evaluate the maximum absolute value of the vector \c v
  template <typename T, typename V> inline
  T maxabsval(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::maxabsval(v); }

  //! Evaluate the maximum value of the vector expression \c e
  template <typename T, typename A> inline
  T maxabsval(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::maxabsval(e); }

  // M I N
  //! Evaluate the minimum value of the vector \c v
  template <typename T, typename V> inline
  T minabsval(VectorBase<T,V> const & v)
  { return F1<VectorBase<T,V> >::minabsval(v); }

  //! Evaluate the maximum value of the vector \c e
  template <typename T, typename A> inline
  T minabsval(VectorE<T,A> const & e)
  { return F1<VectorE<T,A> >::minabsval(e); }

  //@}

  //! \name Vector Dot Product
  //@{

  // D O T
  //! Evaluate dot product between vector \c a and vector \c b 
  template <typename T, typename V1, typename V2> inline
  T dot(VectorBase<T,V1> const & a, VectorBase<T,V2> const & b)
  { return F2<VectorBase<T,V1>,VectorBase<T,V2> >::dot(a,b); }

  //! Evaluate dot product between vector expression \c a and vector \c b 
  template <typename T, typename A, typename V> inline
  T dot(VectorE<T,A> const & a, VectorBase<T,V> const & b)
  { return F2<VectorE<T,A>,VectorBase<T,V> >::dot(a,b); }

  //! Evaluate dot product between vector \c a and vector expression \c b 
  template <typename T, typename A, typename V> inline
  T dot(VectorBase<T,V> const & a, VectorE<T,A> const & b)
  { return F2<VectorBase<T,V>,VectorE<T,A> >::dot(a,b); }

  //! Evaluate dot product between vector expression  \c a and vector expression \c b 
  template <typename T, typename A, typename B> inline
  T dot(VectorE<T,A> const & a, VectorE<T,B> const & b)
  { return F2<VectorE<T,A>,VectorE<T,B> >::dot(a,b); }

  // R D O T
  //! Evaluate dot product between vector \c a and vector \c b
  template <typename T, typename V1, typename V2> inline
  T rdot(VectorBase<T,V1> const & a, VectorBase<T,V2> const & b)
  { return F2<VectorBase<T,V1>,VectorBase<T,V2> >::rdot(a,b); }

  //! Evaluate dot product between vector expression \c a and vector \c b
  template <typename T, typename A, typename V> inline
  T rdot(VectorE<T,A> const & a, VectorBase<T,V> const & b)
  { return F2<VectorE<T,A>,VectorBase<T,V> >::rdot(a,b); }

  //! Evaluate dot product between vector \c a and vector expression \c b
  template <typename T, typename A, typename V> inline
  T rdot(VectorBase<T,V> const & a, VectorE<T,A> const & b)
  { return F2<VectorBase<T,V>,VectorE<T,A> >::rdot(a,b); }

  //! Evaluate dot product between vector expression  \c a and vector expression \c b
  template <typename T, typename A, typename B> inline
  T rdot(VectorE<T,A> const & a, VectorE<T,B> const & b)
  { return F2<VectorE<T,A>,VectorE<T,B> >::rdot(a,b); }

  //@}

  //! \name Vector Distance
  //@{

  // D I S T
  //! Evaluate euclidean distance between vector \c a and vector \c b   
  template <typename T, typename V1, typename V2> inline
  typename return_trait<T>::valueType dist(VectorBase<T,V1> const & a, VectorBase<T,V2> const & b)
  { return F2<VectorBase<T,V1>,VectorBase<T,V2> >::dist(a,b); }

  //! Evaluate euclidean distance between vector expression \c a and vector \c b   
  template <typename T, typename A, typename V> inline
  typename return_trait<T>::valueType dist(VectorE<T,A> const & a, VectorBase<T,V> const & b)
  { return F2<VectorE<T,A>,VectorBase<T,V> >::dist(a,b); }

  //! Evaluate euclidean distance between vector \c a and vector expression \c b   
  template <typename T, typename A, typename V> inline
  typename return_trait<T>::valueType dist(VectorBase<T,V> const & a, VectorE<T,A> const & b)
  { return F2<VectorBase<T,V>,VectorE<T,A> >::dist(a,b); }

  //! Evaluate euclidean distance between vector expression \c a and vector expression \c b   
  template <typename T, typename A, typename B> inline
  typename return_trait<T>::valueType dist(VectorE<T,A> const & a, VectorE<T,B> const & b)
  { return F2<VectorE<T,A>,VectorE<T,B> >::dist(a,b); }

  // D I S T 2

  //! Evaluate square of euclidean distance between vector \c a and vector \c b   
  template <typename T, typename V1, typename V2> inline
  typename return_trait<T>::valueType dist2(VectorBase<T,V1> const & a, VectorBase<T,V2> const & b)
  { return F2<VectorBase<T,V1>,VectorBase<T,V2> >::dist2(a,b); }

  //! Evaluate square of euclidean distance between vector expression \c a and vector \c b   
  template <typename T, typename A, typename V> inline
  typename return_trait<T>::valueType dist2(VectorE<T,A> const & a, VectorBase<T,V> const & b)
  { return F2<VectorE<T,A>,Vector<T> >::dist2(a,b); }

  //! Evaluate square of euclidean distance between vector \c a and vector expression \c b   
  template <typename T, typename A, typename V> inline
  typename return_trait<T>::valueType dist2(VectorBase<T,V> const & a, VectorE<T,A> const & b)
  { return F2<VectorBase<T,V>,VectorE<T,A> >::dist2(a,b); }

  //! Evaluate square of euclidean distance between vector expression \c a and vector expression \c b   
  template <typename T, typename A, typename B> inline
  typename return_trait<T>::valueType dist2(VectorE<T,A> const & a, VectorE<T,B> const & b)
  { return F2<VectorE<T,A>,VectorE<T,B> >::dist2(a,b); }

  //@}

  /*
  //  #####  #     #  ###   #####  #    #        #####  ####### ######  #######
  // #     # #     #   #   #     # #   #        #     # #     # #     #    #
  // #     # #     #   #   #       #  #         #       #     # #     #    #
  // #     # #     #   #   #       ###           #####  #     # ######     #
  // #   # # #     #   #   #       #  #               # #     # #   #      #
  // #    #  #     #   #   #     # #   #        #     # #     # #    #     #
  //  #### #  #####   ###   #####  #    #        #####  ####### #     #    #
  */
  
  //! \name Sorting Structures
  //@{

  /*! \cond NODOC */

  typedef struct { indexType lo, hi; } stack_node;

  static inline
  bool GT( indexType i, indexType i1 )
  { return i>i1; }
    
  static inline
  bool GT( indexType i,  indexType j,
           indexType i1, indexType j1 )
  { return j==j1 ? i>i1 : j>j1; }

  /*! \endcond */

  /*! 
   * Function to sort two vector of indexes and values
   * ascending respect to index vector \c I.
   * \param I           C-array of row index vector
   * \param A           C-array of values
   * \param total_elems total number of element to sort
   * \param MAX_THRESH  threshold value to swich to insert sort
   */
  template <typename I_type, typename T_type>
  static
  void
  QuickSortI(
    I_type    I[],
    T_type    A[],
    indexType total_elems,
    indexType MAX_THRESH = 4
  ) {

    using ::std::swap;

    if ( total_elems <= 1 ) return;
    if ( total_elems <= MAX_THRESH ) goto insert_sort;
    {
      indexType lo = 0;
      indexType hi = total_elems - 1;

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

        indexType Pivot     = I_mid;
        indexType left_ptr  = lo + 1;
        indexType right_ptr = hi - 1;

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

    for ( indexType i = 1; i < total_elems; ++i ) {
      for ( indexType j = i; j > 0; --j ) {
        if ( GT(I[j], I[j-1]) ) break;
        swap( I[j], I[j-1] );
        swap( A[j], A[j-1] );
      }
    }
  }

  /*! 
   * Function to sort two vector of indexes
   * ascending respect to index vector \c I and \c J.
   * \param I           C-array of row index vector
   * \param J           C-array of column index vector
   * \param total_elems total number of element to sort
   * \param MAX_THRESH  threshold value to swich to insert sort
   */
  template <typename I_type>
  static
  void
  QuickSortIJ(
    I_type    I[],
    I_type    J[],
    indexType total_elems,
    indexType MAX_THRESH = 4
  ) {

    using ::std::swap;

    if ( total_elems <= 1 ) return;
    if ( total_elems <= MAX_THRESH ) goto insert_sort;

    {
      indexType lo = 0;
      indexType hi = total_elems - 1;

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

        indexType IPivot    = I_mid;
        indexType JPivot    = J_mid;
        indexType left_ptr  = lo + 1;
        indexType right_ptr = hi - 1;

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

    for ( indexType i = 1; i < total_elems; ++i ) {
      for ( indexType j = i; j > 0; --j ) {
        if ( GT(I[j], J[j], I[j-1], J[j-1]) ) break;
        swap( I[j], I[j-1] );
        swap( J[j], J[j-1] );
      }
    }
  }

  /*! 
   * Function to sort three vector of indexes and values
   * ascending respect to index vector \c I and \c J.
   * \param I           C-array of row index vector
   * \param J           C-array of column index vector
   * \param A           C-array of values
   * \param total_elems total number of element to sort
   * \param MAX_THRESH  threshold value to swich to insert sort
   */
  template <typename I_type, typename T_type>
  static
  void
  QuickSortIJ(
    I_type    I[],
    I_type    J[],
    T_type    A[],
    indexType total_elems,
    indexType MAX_THRESH = 4
  ) {
 
    using ::std::swap;

    if ( total_elems <= 1 ) return;
    if ( total_elems <= MAX_THRESH ) goto insert_sort;

    {
      indexType lo = 0;
      indexType hi = total_elems - 1;

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

        indexType IPivot    = I_mid;
        indexType JPivot    = J_mid;
        indexType left_ptr  = lo + 1;
        indexType right_ptr = hi - 1;

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

    for ( indexType i = 1; i < total_elems; ++i ) {
      for ( indexType j = i; j > 0; --j ) {
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

  //! macro to define the operator overloading of standard matrix-vector expressions
  #define SPARSELIB_MUL_STRUCTURES(MATRIX)                                    \
  template <typename TM, typename T, typename VEC> inline                     \
  Vector_M_mul_V<MATRIX<TM>, VectorBase<T,VEC> >                              \
  operator * (MATRIX<TM> const & M, VectorBase<T,VEC> const & a) {            \
    return Vector_M_mul_V<MATRIX<TM>, VectorBase<T,VEC> >(M,a);              \
  }                                                                           \
                                                                              \
  template <typename T, typename TM> inline                                   \
  Vector_S_mul_M<T,MATRIX<TM> >                                               \
  operator * (T const & s, MATRIX<TM> const & M) {                            \
    return Vector_S_mul_M<T, MATRIX<TM> >(s,M);                              \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VEC> inline                     \
  Vector_S_mul_M_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >                     \
  operator * (Vector_S_mul_M<T,MATRIX<TM> > const & sM,                       \
              VectorBase<T,VEC>             const & v) {                      \
    return Vector_S_mul_M_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >            \
           (sM.s,sM.M,v);                                                    \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VEC> inline                     \
  Vector_S_mul_M_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >                     \
  operator * (T const & s,                                                    \
              Vector_M_mul_V<MATRIX<TM>, VectorBase<T,VEC> > const & Mv) {    \
    return Vector_S_mul_M_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >            \
           (s,Mv.M,Mv.a);                                                    \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VA, typename VB> inline         \
  Vector_V_sum_M_mul_V<VectorBase<T,VA>, MATRIX<TM>, VectorBase<T,VB> >       \
  operator + (VectorBase<T,VA> const & a,                                     \
              Vector_M_mul_V<MATRIX<TM>, VectorBase<T,VB> > const & Mv) {     \
    return Vector_V_sum_M_mul_V<VectorBase<T,VA>, MATRIX<TM>,                 \
                                VectorBase<T,VB> >(a,Mv.M,Mv.a);             \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VA, typename VB> inline         \
  Vector_V_sub_M_mul_V<VectorBase<T,VA>,MATRIX<TM>,VectorBase<T,VB> >         \
  operator - (VectorBase<T,VA> const & a,                                     \
              Vector_M_mul_V<MATRIX<TM>, VectorBase<T,VB> > const & Mv) {     \
    return Vector_V_sub_M_mul_V<VectorBase<T,VA>,MATRIX<TM>,                  \
                                VectorBase<T,VB> >(a,Mv.M,Mv.a);             \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VA, typename VB> inline         \
  Vector_V_sum_S_mul_M_mul_V<VectorBase<T,VA>,T,MATRIX<TM>,VectorBase<T,VB> > \
  operator + (VectorBase<T,VA> const & a,                                     \
              Vector_S_mul_M_mul_V<T,MATRIX<TM>,VectorBase<T,VB> > const & sMv) { \
    return Vector_V_sum_S_mul_M_mul_V<VectorBase<T,VA>,T,MATRIX<TM>,          \
                                      VectorBase<T,VB> >                      \
                                      (a,sMv.s,sMv.M,sMv.a);                 \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VEC> inline                     \
  Vector_Mt_mul_V<MATRIX<TM>, VectorBase<T,VEC> >                             \
  operator ^ (MATRIX<TM> const & M, VectorBase<T,VEC> const & a) {            \
    return Vector_Mt_mul_V<MATRIX<TM>, VectorBase<T,VEC> >(M,a);             \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VEC> inline                     \
  Vector_S_mul_Mt_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >                    \
  operator ^ (Vector_S_mul_M<T,MATRIX<TM> > const & sM,                       \
              VectorBase<T,VEC>             const & v) {                      \
    return Vector_S_mul_Mt_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >           \
           (sM.s,sM.M,v);                                                    \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VEC> inline                     \
  Vector_S_mul_Mt_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >                    \
  operator * (T const & s,                                                    \
              Vector_Mt_mul_V<MATRIX<TM>, VectorBase<T,VEC> > const & Mv) {   \
    return Vector_S_mul_Mt_mul_V<T, MATRIX<TM>, VectorBase<T,VEC> >           \
           (s,Mv.M,Mv.a);                                                    \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VA, typename VB> inline         \
  Vector_V_sum_Mt_mul_V<VectorBase<T,VA>,MATRIX<TM>,VectorBase<T,VB> >        \
  operator + (VectorBase<T,VA> const & a,                                     \
              Vector_Mt_mul_V<MATRIX<TM>, VectorBase<T,VB> > const & Mv) {    \
    return Vector_V_sum_Mt_mul_V<VectorBase<T,VA>,MATRIX<TM>,                 \
                                 VectorBase<T,VB> >(a,Mv.M,Mv.a);            \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VA, typename VB> inline         \
  Vector_V_sub_Mt_mul_V<VectorBase<T,VA>,MATRIX<TM>,VectorBase<T,VB> >        \
  operator - (VectorBase<T,VA> const & a,                                     \
              Vector_Mt_mul_V<MATRIX<TM>,VectorBase<T,VB> > const & Mv) {     \
    return Vector_V_sub_Mt_mul_V<VectorBase<T,VA>,MATRIX<TM>,                 \
                                 VectorBase<T,VB> >(a,Mv.M,Mv.a);            \
  }                                                                           \
                                                                              \
  template <typename T, typename TM, typename VA, typename VB> inline         \
  Vector_V_sum_S_mul_Mt_mul_V<VectorBase<T,VA>,T,MATRIX<TM>,VectorBase<T,VB> >\
  operator + (VectorBase<T,VA> const & a,                                     \
              Vector_S_mul_Mt_mul_V<T,MATRIX<TM>,VectorBase<T,VB> > const & sMv) { \
    return Vector_V_sum_S_mul_Mt_mul_V<VectorBase<T,VA>,T,MATRIX<TM>,         \
                                       VectorBase<T,VB> >                     \
                                       (a,sMv.s,sMv.M,sMv.a);                \
  }

  /*!
   * \class SparseBase
   * \brief
   * This class in the base class for all the
   * sparse matrix classes of \c SparseTool.
   * This class is incomplete and is used as a pivot for internal operations.
   */

  template <typename Matrix>
  class SparseBase {
  protected:

    indexType sp_nrows;     //!< Number of rows of the derived class
    indexType sp_ncols;     //!< Number of columns of the derived class
    indexType sp_min_size;  //!< Minimum between \c sp_nrows and \c sp_ncols
    indexType sp_max_size;  //!< Minimum between \c sp_nrows and \c sp_ncols
    indexType sp_nnz;       //!< Total number of nonzeros of the derived sparse matrix  
    indexType sp_lower_nnz; //!< Total number of nonzeros under the main diagonal
    indexType sp_upper_nnz; //!< Total number of nonzeros over the main diagonal
    indexType sp_diag_nnz;  //!< Total number of nonzeros on the main diagonal
    bool      sp_isOrdered; 
    /*!< \brief Some sparse matrix can be internally in a state not ordered.
         When in this state the random access to element is unpredictable
         and can result in runtime error. 
         This method return \c true if the matrix is ordered. */

    //! check the index \c idx. If out of the range \c 0..nnz-1 and error is issued.
    void
    test_nnz(indexType idx) const {
      SPARSETOOL_TEST(
        idx < sp_nnz,
        "Sparse::operator [" << idx << "] index out of range"
      )
    }
    
    //! check the indices \c i and \c j. If out of the range of the matrix and error is issued.
    void
    test_index(indexType i, indexType j) const {
      SPARSETOOL_TEST(
        i < sp_nrows && j < sp_ncols,
        "Sparse::test_index(" << i << "," << j << ") index out of range"
      )
    }

    //! check the indices \c i. If out of the row range of the matrix and error is issued.
    void
    test_row(indexType i) const {
      SPARSETOOL_TEST(
        i < sp_nrows,
        "Sparse::test_row(" << i << ") index out of range"
      )
    }

    //! check the indices \c j. If out of the column range of the matrix and error is issued.
    void
    test_col(indexType j) const {
      SPARSETOOL_TEST(
        j < sp_ncols,
        "Sparse::test_col(" << j << ") index out of range"
      )
    }

    //! Initialize the class \c Sparse with \c nr rows and \c nc columns
    void
    setup( indexType nr, indexType nc ) {
      sp_nrows     = nr;
      sp_ncols     = nc;
      sp_min_size  = std::min(nr,nc);
      sp_max_size  = std::max(nr,nc);
      sp_nnz       = 0;
      sp_lower_nnz = 0;
      sp_upper_nnz = 0;
      sp_diag_nnz  = 0;
      sp_isOrdered = true;
    }

    //! update the counter for the lower, upper and diagonal elements
    void
    ldu_count(indexType i, indexType j) {
      if      ( j < i ) ++sp_lower_nnz;
      else if ( j > i ) ++sp_upper_nnz;
      else              ++sp_diag_nnz;
    }

  public:

    SparseBase(void) { setup(0, 0); }
    ~SparseBase(void) {};

    // common data
    indexType numRows (void) const { return sp_nrows; }    //!< return the number of rows of the \c Sparse object
    indexType numCols (void) const { return sp_ncols; }    //!< return the number of columns of the \c Sparse object
    indexType minSize (void) const { return sp_min_size; } //!< return the minimum between rows and columns
    indexType maxSize (void) const { return sp_max_size; } //!< return the maximum between rows and columns
    indexType nnz     (void) const { return sp_nnz; }      //!< return the total number of nonzeros

    /*!
     *  \brief
     *  Return the number of nonzeros of the derived sparse
     *  matrix under the main diagonal (\c lower)
     *  over the main diagonal (\c upper) and on the diagonal (\c diag)
     */
    void
    nnz ( indexType & lower, indexType & diag, indexType & upper) const {
      lower = sp_lower_nnz;
      diag  = sp_diag_nnz;
      upper = sp_upper_nnz;
    }
    //! return \c true if the internal data is ordered
    bool isOrdered (void) const { return sp_isOrdered; }

    /*! \brief
        \name Iterator
        These methods are useful for accessing all the nonzero elements
        of the derived sparse matrix.  The methods 

        - \c void \c Begin()
        - \c void \c Next()
        - \c bool \c End()
        
        permits to loops on all the elements, while the methods

        - \c indexType \c row()
        - \c indexType \c column()
        - \c valueType \c value()

        permits to access values of the actual elements 
        pointed by the iterator. For example to print all the stored
        values of the \c Sparse object \c S we can do:

\code
  for ( S.Begin(); S.End(); S.Next() ) {
    cout << " row   = " << S.row()
         << " col   = " << S.column()
         << " value = " << S.value()
         << '\n';  
  }
\endcode
    */

    //! \name Iterator for accessing matrix elements
    //@{
    //! set the iterator at the begin of the loop
    void Begin (void) const { static_cast<Matrix const *>(this) -> Begin(); }
    //! go the next item
    void Next  (void) const { static_cast<Matrix const *>(this) -> Next(); }
    //! return \c true unless we are at the end of the loop
    bool End   (void) const { return static_cast<Matrix const *>(this) -> End(); }

    //! The \a row of the pointed element
    indexType row   (void) const { return static_cast<Matrix const *>(this) -> row(); }

    //! The \a column of the pointed element
    indexType column(void) const { return static_cast<Matrix const *>(this) -> column(); }

    //! Assign the pointed element to \c rhs
    template <typename T>
    void assign( T & rhs ) const {
      Matrix const * This = static_cast<Matrix const *>(this);
      This -> assign(rhs);
    }

    //@}

    //! return \c true if the element \c (i,j) is a nonzeros of the derived matrix
    bool
    exists( indexType i, indexType j )
    { return static_cast<Matrix const *>(this) -> exists( i, j ); }

    /*! \brief
        Return the position of the \c (i,j) elements
        in the internal structure. If the element \c (i,j) do not
        exist it return \c nnz() */
    indexType
    position( indexType i, indexType j ) const
    { return static_cast<Matrix const *>(this) -> position( i, j ); }

  };

  /*!
   * \class Sparse
   * \brief
   * This class in the base class for all the
   * sparse matrix classes of \c SparseTool.
   * This class is incomplete and is used as a pivot for internal operations.
   */

  template <typename T, typename Matrix>
  class Sparse : public SparseBase<Matrix> {
    typedef SparseBase<Matrix> SBASE;
  public:
    typedef T valueType; //!< type of the element of the matrix

    Sparse(void) { SBASE::setup(0, 0); }
    ~Sparse(void) {};

    //! The \a value of the pointed element
    valueType value (void) const { return static_cast<Matrix const *>(this) -> value(); }

    //! Assign the pointed element to \c rhs
    template <typename TS>
    void assign( TS & rhs ) const { rhs = this -> value(); }

    //@}

    /*!
     * \brief
     * \name Random access to elements
     * The following methods permits random access to
     * nonzeros elements of the derived class.
     */
    //@{
    /*!
     * \brief
     * Access to the elements of the sparse matrix in
     * a random way.  For example given a sparse matrix \c A the code
     * \c A(i,j) return the \a value of the elements of the
     * matrix at the \c i-th row and \c j-th column.
     * If at the \c (i,j) coordinate there are no elements an error is issued.
     */
    valueType const &
    operator () ( indexType i, indexType j ) const
    { return static_cast<Matrix const *>(this) -> operator () (i,j); }

    /*!
     * \brief
     * Access to the elements of the sparse matrix in
     * a random way.  For example given a sparse matrix \c A the code
     * \c A(i,j) return the \a value of the elements of the
     * matrix at the \c i-th row and \c j-th column.
     * If at the \c (i,j) coordinate there are no elements an error is issued.
     */
    valueType &
    operator () ( indexType i, indexType j )
    { return static_cast<Matrix *>(this) -> operator () (i,j); }

    /*!
     * \brief
     * This method permits to access the elements of the sparse matrix in
     * a random way.  For example given a sparse matrix \c A the code
     * \c A(i,j) return the \a reference of the elements of the
     * matrix at the \c i-th row and \c j-th column.  If at the \c (i,j)
     * coordinate there are no elements a reference to a value \b 0 is returned
     */
    valueType const &
    value( indexType i, indexType j ) const
    { return static_cast<Matrix const *>(this) -> value(i,j); }

    //! return \c true if the element \c (i,j) is a nonzeros of the derived matrix
    bool
    exists( indexType i, indexType j )
    { return static_cast<Matrix const *>(this) -> exists( i, j ); }

    /*!
     * \brief
     * Return the position of the \c (i,j) elements
     * in the internal structure. If the element \c (i,j) do not
     * exist it return \c nnz()
     */
    indexType
    position( indexType i, indexType j ) const
    { return static_cast<Matrix const *>(this) -> position( i, j ); }

    //! Return the reference of the \c idx-th element of the vector of stored values
    valueType const &
    operator [] (indexType idx) const
    { return static_cast<Matrix const *>(this) -> operator[] (idx); }

    //! Return the reference of the \c idx-th element of the vector of stored values
    valueType &
    operator [] (indexType idx)
    { return static_cast<Matrix *>(this) -> operator[] (idx); }

    //@}

    //! set all the nonzeros of the sparse matrix to the value \c 0
    void
    setZero()
    { return static_cast<Matrix const *>(this) -> setZero(); }

    //! multiply all the nonzeros of the sparse matrix by \c s
    void
    scaleValues( valueType const & s )
    { return static_cast<Matrix const *>(this) -> scaleValues( s ); }

    //! multiply all the nonzeros of the row \c nr matrix by \c val
    void
    scaleRow( indexType nr, valueType const & val )
    { return static_cast<Matrix const *>(this) -> scaleRow( nr, val ); }

    //! multiply all the nonzeros of the column \c nc matrix by \c val
    void
    scaleColumn( indexType nc, valueType const & val )
    { return static_cast<Matrix const *>(this) -> scaleColumn( nc, val ); }

    /*! \name Diagonal internal operation */
    //@{
    /*! \brief
        Set \c s to all the components of the diagonal,
        and \c 0 the others elements.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
            2 & -1 &    & 1  &    \\ 
            1 & 2  &    &    & 5  \\ 
              &    & 2  & -3 &    \\ 
              &    & 3  & 2  & -4 \\ 
              &    &    & 4  & 2
        \end{BMAT}\right],\qquad
        (\bm{A}=2.1) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
           2.1 & 0   &     & 0   &    \\ 
           0   & 2.1 &     &     & 0  \\ 
               &     & 2.1 & 0   &    \\ 
               &     & 0   & 2.1 & 0  \\ 
               &     &     & 0   & 2.1
        \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    Matrix &
    operator = ( valueType const & s ) {
      setZero();
      if ( s != 0 )
        for ( indexType k = 0; k < SBASE::sp_min_size; ++k )
          (*this)(k,k) = s;
      return *this;
    }
    /*! \brief Add \c s to all the components of the diagonal.
        For example
      \htmlonly <TABLE><TR><TD> \endhtmlonly
      \f[
      \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      2 & -1 &    & 1  &    \\ 
      1 & 2  &    &    & 5  \\ 
        &    & 2  & -3 &    \\ 
        &    & 3  & 2  & -4 \\ 
        &    &    & 4  & 2
      \end{BMAT}\right],\qquad
      (\bm{A}\verb|+=| 2) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      4 & -1 &    & 1  &    \\ 
      1 & 4  &    &    & 5  \\ 
        &    & 4  & -3 &    \\ 
        &    & 3  & 4  & -4 \\ 
        &    &    & 4  & 4
      \end{BMAT}\right]
      \f]
      \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    Matrix &
    operator += ( valueType const & s ) {
      for ( indexType k = 0; k < SBASE::sp_min_size; ++k )
        (*this)(k,k) += s;
      return *this;
    }
    /*! \brief Subtract \c s to all the components of the diagonal.
        For example
      \htmlonly <TABLE><TR><TD> \endhtmlonly
      \f[
      \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      2 & -1 &    & 1  &    \\ 
      1 & 2  &    &    & 5  \\ 
        &    & 2  & -3 &    \\ 
        &    & 3  & 2  & -4 \\ 
        &    &    & 4  & 2
      \end{BMAT}\right],\qquad
      (\bm{A} \verb|-=| 2) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      0 & -1 &    & 1  &    \\ 
      1 & 0  &    &    & 5  \\ 
        &    & 0  & -3 &    \\ 
        &    & 3  & 0  & -4 \\ 
        &    &    & 4  & 0
      \end{BMAT}\right]
      \f]
      \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    Matrix &
    operator -= ( valueType const & s ) {
      for ( indexType k = 0; k < SBASE::sp_min_size; ++k )
        (*this)(k,k) -= s;
      return *this;
    }
    /*! \brief Multiply by \c s to all the components of the diagonal.
        For example
      \htmlonly <TABLE><TR><TD> \endhtmlonly
      \f[
      \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      2 & -1 &    & 1  &    \\ 
      1 & 2  &    &    & 5  \\ 
        &    & 2  & -3 &    \\ 
        &    & 3  & 2  & -4 \\ 
        &    &    & 4  & 2
      \end{BMAT}\right],\qquad
      (\bm{A} \verb|*=| 2) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      4 & -2 &    & 2  &    \\ 
      2 &  4 &    &    & 10 \\ 
        &    & 4  & -6 &    \\ 
        &    & 6  & 4  & -8 \\ 
        &    &    & 8  & 4
      \end{BMAT}\right]
      \f]
      \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    Matrix &
    operator *= ( valueType const & s )
    { scaleValues(s); return * this; }
    /*! \brief Divide all the nonzeros of the derived matrix by \c s.
        For example
      \htmlonly <TABLE><TR><TD> \endhtmlonly
      \f[
      \bm{A}= \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      2 & -1 &    & 1  &    \\ 
      1 & 2  &    &    & 5  \\ 
        &    & 2  & -3 &    \\ 
        &    & 3  & 2  & -4 \\ 
        &    &    & 4  & 2
      \end{BMAT}\right],\qquad
      (\bm{A} \verb|/=| 10) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
      0.2 & -0.1 &      & 0.1  &      \\ 
      0.1 & 0.2  &      &      & 0.5  \\ 
          &      & 0.2  & -0.3 &      \\ 
          &      & 0.3  & 0.2  & -0.4 \\ 
          &      &      & 0.4  & 0.2
      \end{BMAT}\right]
      \f]
      \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    Matrix &
    operator /= ( valueType const & s)
    { valueType ss = 1/s; scaleValues(ss); return * this; }

    // ------------------------
    /*! \brief
        Set the diagonal values of the sparse
        matrix to \c v.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
           \bm{A}= \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
           2 & -1 &    & 1  &    \\ 
           1 & 2  &    &    & 5  \\ 
             &    & 2  & -3 &    \\ 
             &    & 3  & 2  & -4 \\ 
             &    &    & 4  & 2
           \end{BMAT}\right],\qquad
           \bm{v} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
           1 \\ 2\\ 3
           \end{BMAT}\right],
           \qquad
           (\bm{A}=\bm{v}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
           1 & 0 &   & 0 &   \\ 
           0 & 2 &   &   & 0 \\ 
             &   & 3 & 0 &   \\ 
             &   & 0 & 0 & 0 \\ 
             &   &   & 0 & 0
          \end{BMAT}\right]
        \f]
       \htmlonly </TD></TR></TABLE> \endhtmlonly
       notice that only the first \c v.size() diagonal elements are
       set while the rest of the matrix is set to \c 0.
    */
    template <typename V> inline
    Matrix &
    operator = ( VectorBase<T,V> const & v ) {
      indexType n = minIndex(v.size(),  SBASE::sp_min_size );
      setZero();
      for ( indexType k = 0; k < n; ++k )
        (*this)(k,k) = v(k);
      return *this;
    }
    /*! \brief 
        Add the components of \c v to the diagonal.
        For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
          2 & -1 &    & 1  &    \\ 
          1 & 2  &    &    & 5  \\ 
            &    & 2  & -3 &    \\ 
            &    & 3  & 2  & -4 \\ 
            &    &    & 4  & 2
        \end{BMAT}\right],\qquad
        \bm{v} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
          1 \\ 2\\ 3
        \end{BMAT}\right],
        \qquad
        (\bm{A}\verb|-=| \bm{v}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
          3 & -1 &    & 1  &    \\ 
          1 & 4  &    &    & 5  \\ 
            &    & 5  & -3 &    \\ 
            &    & 3  & 2  & -4 \\ 
            &    &    & 4  & 2
        \end{BMAT}\right]
      \f]
      \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename V> inline
    Matrix &
    operator += ( VectorBase<T,V> const & v ) {
      indexType n = minIndex(v.size(),  SBASE::sp_min_size );
      for ( indexType k = 0; k < n; ++k )
        (*this)(k,k) += v(k);
      return *this;
    }
    /*! \brief 
        Subtract the components of \c v to the diagonal.
        For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
          2 & -1 &    & 1  &    \\ 
          1 & 2  &    &    & 5  \\ 
            &    & 2  & -3 &    \\ 
            &    & 3  & 2  & -4 \\ 
            &    &    & 4  & 2
        \end{BMAT}\right],\qquad
        \bm{v} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
          1 \\ 2\\ 3
        \end{BMAT}\right],
        \qquad
        (\bm{A}\verb|-=| \bm{v}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
          1 & -1 &    & 1  &    \\ 
          1 & 0  &    &    & 5  \\ 
            &    & -1  & -3 &    \\ 
            &    & 3  & 2  & -4 \\ 
            &    &    & 4  & 2
        \end{BMAT}\right]
      \f]
      \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename V> inline
    Matrix &
    operator -= ( VectorBase<T,V> const & v ) {
      indexType n = minIndex(v.size(),  SBASE::sp_min_size);
      for ( indexType k = 0; k < n; ++k )
        (*this)(k,k) -= v(k);
      return *this;
    }
    /*! \brief
        Assign the vector expression \c e to the components of the diagonal,
        and \c 0 the others elements.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
            2 & -1 &    & 1  &    \\ 
            1 & 2  &    &    & 5  \\ 
              &    & 2  & -3 &    \\ 
              &    & 3  & 2  & -4 \\ 
              &    &    & 4  & 2
        \end{BMAT}\right],\qquad
        \bm{a} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c.c.c}
            1 \\ 2 \\ 3 \\ 4 \\ 5
        \end{BMAT}\right],\qquad
        \bm{b} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
            1 \\ -1 \\ 2
        \end{BMAT}\right],
        \f]
        \htmlonly </TD><TD>&nbsp;&nbsp;&nbsp;</TD><TD> \endhtmlonly
        \f[
        (\bm{A}=\bm{a}+2\bm{b}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
           3 & 0 &   & 0  &    \\ 
           0 & 0 &   &    & 0  \\ 
             &   & 7 & 0  &    \\ 
             &   & 0 & 0  & 0  \\ 
             &   &   & 0  & 0
         \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename R> inline
    Matrix &
    operator = ( VectorE<T,R> const & e ) {
      indexType n = minIndex(e.size(),  SBASE::sp_min_size );
      setZero();
      for ( indexType k = 0; k < n; ++k )
        (*this)(k,k) = e(k);
      return *this;
    }
    /*! \brief
        Add the vector expression \c e to the components of the diagonal,
        and \c 0 the others elements.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
            2 & -1 &    & 1  &    \\ 
            1 & 2  &    &    & 5  \\ 
              &    & 2  & -3 &    \\ 
              &    & 3  & 2  & -4 \\ 
              &    &    & 4  & 2
        \end{BMAT}\right],\qquad
        \bm{a} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c.c.c}
            1 \\ 2 \\ 3 \\ 4 \\ 5
        \end{BMAT}\right],\qquad
        \bm{b} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
            1 \\ -1 \\ 2
        \end{BMAT}\right],
        \f]
        \htmlonly </TD><TD>&nbsp;&nbsp;&nbsp;</TD><TD> \endhtmlonly
        \f[
        (\bm{A}\verb|+=| \bm{a}+2\bm{b}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
            5 & -1 &    & 1  &    \\ 
            1 & 2  &    &    & 5  \\ 
              &    & 9  & -3 &    \\ 
              &    & 3  & 2  & -4 \\ 
              &    &    & 4  & 2
        \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename R> inline
    Matrix &
    operator += ( VectorE<T,R> const & e ) {
      indexType n = minIndex(e.size(),  SBASE::sp_min_size );
      for ( indexType k = 0; k < n; ++k )
        (*this)(k,k) += e(k);
      return *this;
    }
    /*! \brief
        Add the vector expression \c e to the components of the diagonal,
        and \c 0 the others elements.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
            2 & -1 &    & 1  &    \\ 
            1 & 2  &    &    & 5  \\ 
              &    & 2  & -3 &    \\ 
              &    & 3  & 2  & -4 \\ 
              &    &    & 4  & 2
        \end{BMAT}\right],\qquad
        \bm{a} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c.c.c}
            1 \\ 2 \\ 3 \\ 4 \\ 5
        \end{BMAT}\right],\qquad
        \bm{b} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
            1 \\ -1 \\ 2
        \end{BMAT}\right],
        \f]
        \htmlonly </TD><TD>&nbsp;&nbsp;&nbsp;</TD><TD> \endhtmlonly
        \f[
        (\bm{A} \verb|-=| \bm{a}+2\bm{b}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c}{c.c.c.c.c}
           -1 & -1 &    & 1  &    \\ 
            1 & 2  &    &    & 5  \\ 
              &    & -5 & -3 &    \\ 
              &    & 3  & 2  & -4 \\ 
              &    &    & 4  & 2
        \end{BMAT}\right]
        \f]    
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename R> inline
    Matrix &
    operator -= ( VectorE<T,R> const & e ) {
      indexType n = minIndex(e.size(),  SBASE::sp_min_size );
      for ( indexType k = 0; k < n; ++k )
        (*this)(k,k) -= e(k);
      return *this;
    }
    //@}

    /*! \brief 
        \name Matrix-Vector multiplication
    */

    //@{

    //! perform the operation res += s * (A * x)
    template <typename VA, typename VB>
    void
    add_S_mul_M_mul_V(
      VectorBase<T,VA>       & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      SPARSETOOL_TEST(
        (void*)&res != (void*)&x,
        "add_S_mul_M_mul_V equal pointer"
      )
      SPARSETOOL_TEST(
        res.size() >= SBASE::sp_nrows,
        "result vector too small"
      )
      for ( SBASE::Begin(); SBASE::End(); SBASE::Next() )
        res(SBASE::row()) += s * value() * x(SBASE::column());
    }
    //! perform the operation res += s * (A ^ x)
    template <typename VA, typename VB>
    void
    add_S_mul_Mt_mul_V(
      VectorBase<T,VA>       & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      SPARSETOOL_TEST(
        (void*)&res != (void*)&x,
        "add_S_mul_M_mul_V equal pointer"
      )
      SPARSETOOL_TEST(
        res.size() >= SBASE::sp_ncols,
        "result vector too small"
      )
      for ( SBASE::Begin(); SBASE::End(); SBASE::Next() )
        res(SBASE::column()) += s * value() * x(SBASE::row());
    }
    //@}

  };

  /*!
   * \class Preco
   *
   * This class in the base class for all the preconditioner
   * of sparse matrix classes of \c SparseTool.
   * This class is incomplete and is used as a pivot for internal operations.
   */
  //! Base preconditioner class
  template <typename PRECO>
  class Preco {
  protected:
    indexType pr_size; //!< the number of row/column of the preconditioner matrix
  public:
    //! Create an empty preconditioner
    Preco(void) : pr_size(0) {}
    ~Preco(void) {};
    //! return the number of row/column of the preconditioner matrix 
    indexType size(void) const { return pr_size; }
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

  //! Sparse Pattern in compressed coordinate.
  /*!
    The sparse pattern internal structure is essentially a sparse
    compressed coordinate ones.  It consists of two big vector of unsigned
    integer which contain the coordinate of nonzero elements.  We call
    \c I the vector that store the first coordinate, while we call
    \c J the vector that store the second coordinate.  For example the
    following \b 6 x \b 7 sparse matrix pattern

    \htmlonly <TABLE><TR><TD> \endhtmlonly
    \f[ 
      \bm{S} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c.c.c}{c.c.c.c.c.c}
     * & * &   &   &   &   & * \\ 
       & * & * &   &   &   &   \\ 
       &   &   & * & * &   &   \\ 
       & * &   &   & * &   &   \\ 
     * &   & * &   &   & * &   \\ 
       & * &   &   &   &   & * 
       \end{BMAT}\right]
    \f]
    \htmlonly </TD><TD>&nbsp;&nbsp;&nbsp;</TD><TD> \endhtmlonly
    \f[
       \begin{tabular}{c|cccccccccccccc}
         I & 0 & 0 & 0 & 1 & 1 & 2 & 2 & 3 & 3 & 4 & 4 & 4 & 5 & 5 \\
         \hline
         J & 0 & 1 & 6 & 1 & 2 & 3 & 4 & 1 & 4 & 0 & 2 & 5 & 1 & 6
       \end{tabular}
    \f]
    \htmlonly </TD></TR></TABLE> \endhtmlonly
    Notice that the index follows the C convention, starting from \b 0. 
    The class \c SparsePattern try to manage such a structure in a
    simple way for the user.  To define a \c SparsePattern class you
    can use uno of the following scripture

\code
SparsePatter sp; // define an empty SparsePattern class
SparsePatter sp(nr, nc, nnz); // define a SparsePattern class on a pattern of nr rows and nc columns.
                              // The number nnz is an unsigned integer which determine the preallocated
                              // number of pattern elements stored in the class, in practice is the
                              // initial dimension of vectors I and J.
SparsePatter sp(sp1);         // define a SparsePattern object which is the copy of SparsePattern object sp1.
SparsePatter sp(sobj);        // define a SparsePattern object which is the pattern of nonzero of the
                              // Sparse object sobj. In this way it is possible to obtain the sparse
                              // pattern of any object derived from Sparse.
\endcode

    It is possible in any moment to change the sizes and the maximum 
    number of nonzero by the \c resize methods:

\code
sp.resize(nr, nc, nnz);
sp.resize(sp1);
sp.resize(sobj);
\endcode

    to complete the basic description of the \c SparsePattern a sintetic
    description of the remaining methods of the class are presented

    - \c insert(i,j) this method permit to insert an item of nonzero.  For example
      the first pattern can be constructed as
\code
  SparsePattern sp(6,7);
  sp.insert(0, 0); sp.insert(0, 1); sp.insert(0, 6); sp.insert(1, 1);
  sp.insert(1, 2); sp.insert(2, 3); sp.insert(2, 4); sp.insert(3, 1);
  sp.insert(3, 4); sp.insert(4, 0); sp.insert(4, 2); sp.insert(4, 5);
  sp.insert(5, 1); sp.insert(5, 6); sp.internalOrder();
\endcode

    - \c internalOrder()
      This methods reorder internally the nonzero elements of the sparse pattern
      and remove all duplicated entries.
      If <c> k1 <= k2 </c> we have one of the two following cases

      -# <c> J(k1) < J(k2) </c>
      -# <c> J(k1) == J(k2) </c> and <c> I(k1) <= I(k2) </c>

    - \c bool \c isOrdered()
      this methods return \c true if the elements inside the sparse
      pattern are ordered, \c false otherwise.

  */
  class SparsePattern : public SparseBase<SparsePattern> {
    typedef SparsePattern             MATRIX;
    typedef SparseBase<SparsePattern> SPARSE;

    Vector<indexType> I;    //!< Vector of row index 
    Vector<indexType> J;    //!< Vector of column index
    mutable indexType ipos; //!< Actual position of iterator
    
    //! \cond NODOC
    template <typename MAT, typename Compare>
    void
    convert( SparseBase<MAT> const & M, Compare cmp ) {
      // count nonzero
      indexType nz = 0;
      for ( M.Begin(); M.End(); M.Next() )
        if ( cmp( M.row(), M.column() ) ) ++nz;

      resize( M.numRows(), M.numCols(), nz );

      for ( M.Begin(); M.End(); M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) insert(i,j);
      }
      internalOrder();
    }
    //! \endcond

  public:

    /*!
     *  \brief
     *  Initialize and empty sparse pattern of \c 0 rows and \c 0 columns.
     */
    SparsePattern(void)
    { }

    /*!
     * \brief
     * Initialize and empty sparse pattern of \c nr rows and \c nc columns
     * with reserved room for \c mnnz nonzeros.
     */
    SparsePattern(
      indexType nr,
      indexType nc,
      indexType mnnz = SPARSETOOL_DEFAULT_NNZ
    )
    { resize(nr, nc, mnnz); }

    //! Assign the pointed element to \c rhs
    template <typename T>
    void assign( T & rhs ) const { /* do nothing! */ }

    /*!
     *  \brief
     *  Insert the element of the sparse pattern \c sp
     *  which satify \c cmp to \c *this.
     *  \param sp  sparse pattern to be copied
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename Compare>
    SparsePattern( SparsePattern const & sp, Compare cmp )
    { convert(sp,cmp); }

    /*!
     * \brief
     * Copy the sparse pattern \c sp to \c *this.
     * \param sp sparse pattern to be copied
     */
    SparsePattern(SparsePattern const & sp) : SPARSE() 
    { convert(sp,all_ok()); }

    // convert => SparsePattern
    /*!
     * \brief
     * Insert the element of the pattern of \c M
     * which satify \c cmp to \c *this.
     * \param M   object derived from \c Sparse
     * \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    SparsePattern(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    /*!
     * \brief
     * Copy the sparse pattern of \c M to \c *this.
     * \param M object derived from \c Sparse
     */
    template <typename MAT>
    SparsePattern(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    /*!
     * \brief
     * Copy the sparse pattern of \c SP to \c *this.
     * \param SP sparse pattern object
     */
    SparsePattern &
    operator = ( SparsePattern const & SP ) {
      if ( &SP != this ) { // avoid copy to itself
        SPARSE::setup( SP.numRows(), SP.numCols() );
        SPARSE::sp_nnz       = SP.nnz();
        SPARSE::sp_isOrdered = SP.isOrdered();
        I.load( SP.I );
        J.load( SP.J );
        if ( !SPARSE::sp_isOrdered ) internalOrder();
      }
      return *this;
    }

    // convert => SparsePattern
    /*!
     * \brief
     * Copy the sparse pattern of \c M to \c *this.
     * \param M object derived from \c Sparse
     */
    template <typename MAT>
    SparsePattern &
    operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return *this; }

    /*!
     * \brief
     * Initialize and empty sparse pattern of \c nr rows and \c nc columns
     * with reserved room for \c mnnz nonzeros.
     */
    void
    resize( indexType nr, indexType nc, indexType mnnz = SPARSETOOL_DEFAULT_NNZ ) {
      SPARSE::setup( nr, nc );
      I.reserve( mnnz ); I.resize( 0 );
      J.reserve( mnnz ); J.resize( 0 );
    }

    /*! \brief
     *  Initialize the pattern and insert the element of the pattern of \c M
     *  to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    void
    resize(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    /*! \brief
     *  Initialize the pattern and insert the elements M(i,j)
     *  which satify \c cmp(i,j) to \c *this.
     *  \param M   object derived from \c Sparse
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    void
    resize(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Order the internal structure to permit random
     *  access to nonzeros.
     */
    void
    internalOrder() {
      QuickSortIJ<indexType>( &I.front(), &J.front(), SPARSE::sp_nnz );
      // eliminate duplicate elements
      indexType i1 = 0, i = 0;
      SPARSE::sp_lower_nnz = 0;
      SPARSE::sp_diag_nnz  = 0;
      SPARSE::sp_upper_nnz = 0;
      while ( i1 < SPARSE::sp_nnz ) {
        // copy i1 => i
        I(i) = I(i1); J(i) = J(i1);
        // setup statistic
        SPARSE::ldu_count(I(i),J(i));
        // find i1 != i
        for ( ++i1; i1 < SPARSE::sp_nnz && I(i1) == I(i) && J(i1) == J(i); ++i1 ) {};
        ++i;
      }
      SPARSE::sp_nnz       = i;
      SPARSE::sp_isOrdered = true;
      I.resize( SPARSE::sp_nnz );
      J.resize( SPARSE::sp_nnz );
    }

    //! check if pattern is empty
    void empty(void) { SPARSE::sp_nnz = 0; }

    //! Insert \c (i,j) in the sparse pattern
    SparsePattern &
    insert( indexType i, indexType j ) {
      SPARSE::test_index(i,j);
      I.push_back(i);
      J.push_back(j);
      ++SPARSE::sp_nnz;
      SPARSE::sp_isOrdered = false;
      return *this;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<indexType> const & getI(void) const { return I; } //!< return the row index vector
    Vector<indexType> const & getJ(void) const { return J; } //!< return the column index vector

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // ITERATORS
    void Begin (void) const { ipos = 0; } //!< initialize iterator
    void Next  (void) const { ++ipos; }   //!< go to the next element
    bool End   (void) const { return ipos < SPARSE::sp_nnz; } //! check if iterator terminate the loop

    indexType row    (void) const { return I(ipos); } //! the row of the element pointed by iterator
    indexType column (void) const { return J(ipos); } //! the column of the element pointed by iterator

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
  
  /*!

    \class CCoorMatrix

    The class \c CCoorMatrix\<T\> implement a 
    <B> Compressed Coordinate </B> storage sparse scheme.
    It consists of two big vector of unsigned integer which contain
    the coordinate of nonzero elements and a big one of real number
    which contain the values.  We call \c I the vector that store the
    rows coordinate, \c J the vector that store the columns coordinate 
    and \c A the vector that store the nonzero values.
    For example the following \b 6 x \b 7 sparse matrix

    \htmlonly <TABLE><TR><TD> \endhtmlonly
    \f[
    \bm{A} = \left[\begin{BMAT}(e){c.c.c.c.c.c.c}{c.c.c.c.c.c}
     1 & 2    &    &   &   &   & 9 \\ 
       & -1   & 0  &   &   &   &   \\ 
       &      &    & 3 & 4 &   &   \\ 
       & 2    &    &   & 5 &   &   \\ 
     2 &      & -2 &   &   & 1 &   \\ 
       & -1.5 &    &   &   &   & -1 
    \end{BMAT}\right]
    \f]
    \htmlonly </TD><TD>&nbsp;&nbsp;&nbsp;</TD><TD> \endhtmlonly
    \f[
    \begin{tabular}{c|cccccccccccccc}
       I & 0 & 0 & 0 & 1 & 1 & 2 & 2 & 3 & 3 & 4 & 4 & 4 &  5  & 5 \\
       \hline
       J & 0 & 1 & 6 & 1 & 2 & 3 & 4 & 1 & 4 & 0 & 2 & 5 &  1  & 6 \\
       \hline
       A & 1 & 2 & 9 &-1 & 0 & 3 & 4 & 2 & 5 & 2 &-2 & 1 &-1.5 &-1
    \end{tabular}
    \f]
    \htmlonly </TD></TR></TABLE> \endhtmlonly
    
    Notice that the index follows the C convention, starting from
    \b 0.  The class \c CCoorMatrix\<T\> try to manage such a
    structure in a simple way for the user.  To define a
    \c CCoorMatrix\<T\> class you can use one of the following scripture

\code
CCoorMatrix<double> ccoor; // define an empty CCoorMatrix<double> class

// instance a CCoorMatrix<double> class of nr rows and nc columns.
// The number nnz is an unsigned integer which determine the pre-allocated number of elements
// stored in the class, in practice it is the initial dimension of vectors I, J and \c A.
//  If needed the vectors are enlarged.
CCoorMatrix<double> ccoor(nr, nc, nnz);

// instance a CCoorMatrix<double> class with the sparsity pattern defined of the SparsePattern class sp.
// nr will be sp.numRows(), ncol will be sp.numCols().
CCoorMatrix<double> ccoor(sp);

// instance a CCoorMatrix<double> class with which is the copy of the CCoorMatrix<double> object ccoor1.
CCoorMatrix<double> ccoor(ccoor1);

// instance a CCoorMatrix<double> class with which is the copy of the Sparse object sobj.
CCoorMatrix<double> ccoor(sobj);
\endcode

    It is possible in any moment to change the sizes and the maximum
    number of nonzero by the \c resize methods:

\code
  ccoor.resize(nrow, ncol, nnz);
  ccoor.resize(sp);
  ccoor.resize(ccoor1);
  ccoor.resize(sobj);
\endcode

    - <c>valueType & insert(i,j)</c>
    this method permit to insert an item in the matrix.  For example
    matrix \c A previously defined can be constructed with
\code
  CCoorMatrix<double> A(6,7);
  A.insert(0, 0) = 1;    A.insert(0, 1) = 2; A.insert(0, 6) = 9;  A.insert(1, 1) = -1;
  A.insert(1, 2) = 0;    A.insert(2, 3) = 3; A.insert(2, 4) = 4;  A.insert(3, 1) = 2;
  A.insert(3, 4) = 5;    A.insert(4, 0) = 2; A.insert(4, 2) = -2; A.insert(4, 5) = 1;
  A.insert(5, 1) = -1.5; A.insert(5, 6) = 1; A.internalOrder();
\endcode

    - \c internalOrder()
      This methods reorder internally the nonzero elements of
      \c CCoorMatrix\<double\> in such a way if 
      <c> k1 < k2 </c>we have one of the two following cases
      -# <c> I(k1) <  I(k2) </c>
      -# <c> I(k1) == I(k2) </c> and <c> J(k1) <= J(k2)</c>
      \n
      Moreover all duplicated entries are added togheter.
    
    - <c> bool isOrdered() </c>
      this methods return \c true if the elements inside the class
      <c> CCoorMatrix<double> </c>are ordered, \c false otherwise. 
  
  */
  //! Compressed Coordinate Matrix Storage
  template <typename T>
  class CCoorMatrix : public Sparse<T,CCoorMatrix<T> > {
    typedef CCoorMatrix<T>   MATRIX;
    typedef Sparse<T,MATRIX> SPARSE;
  public:
    typedef T valueType; //!< the type of the elements of the matrix

  private:

    Vector<indexType> I, J, C;
    Vector<valueType> A;

    mutable indexType ipos;

    template <typename MAT, typename Compare> inline
    void
    convert( SparseBase<MAT> const & M, Compare cmp ) {
      // count nonzero
      indexType nz = 0;
      for ( M.Begin(); M.End(); M.Next() )
        if ( cmp( M.row(), M.column() ) ) ++nz;

      resize( M.numRows(), M.numCols(), nz );

      for ( M.Begin(); M.End(); M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) M.assign ( insert(i,j) );
      }
      internalOrder();
    }

  public:

    /*! \brief 
     *  Initialize and empty sparse compressed coordinate matrix
     *  of \c 0 rows and \c 0 columns.
     */
    CCoorMatrix(void) {};

    /*! \brief
     *  Initialize and empty empty sparse compressed coordinate matrix
     *  of \c nr rows and \c nc columns with reserved room for \c mnnz nonzeros.
     */
    CCoorMatrix( indexType nr, indexType nc, indexType mnnz = SPARSETOOL_DEFAULT_NNZ )
    { resize(nr, nc, mnnz); }

    /*! \brief
     *  Insert the element of the sparse compressed coordinate matrix \c M
     *  which satify \c cmp to \c *this.
     *  \param M   sparse compressed coordinate matrix to be copied
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename Compare>
    CCoorMatrix( CCoorMatrix<T> const & M, Compare cmp )
    { convert(M,cmp); }

    /*! \brief
     *  Copy the sparse compressed coordinate matrix \c M to \c *this.
     *  \param M sparse compressed coordinate matrix to be copied
     */
    CCoorMatrix(CCoorMatrix<T> const & M)
    { convert(M,all_ok()); }

    /*! \brief
     *  Insert the element of the sparse compressed coordinate matrix \c M
     *  which satify \c cmp to \c *this.
     *  \param M   object derived from \c Sparse
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    CCoorMatrix(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Copy the sparse compressed coordinate matrix \c M to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    CCoorMatrix(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    /*! \brief
     *  Copy the sparse compressed coordinate matrix \c M to \c *this.
     *  \param M sparse compressed coordinate matrix
     */
    CCoorMatrix<T> &
    operator = (CCoorMatrix<T> const & M) {
      if ( &M == this ) return *this; // avoid copy to itself
      resize( M.numRows(), M.numCols(), M.nnz() );
      SPARSE::sp_nnz       = M.nnz();
      SPARSE::sp_isOrdered = M.isOrdered();
      A.load( M.A );
      I.load( M.I );
      J.load( M.J );
      C.load( M.C );
      if ( !SPARSE::sp_isOrdered ) internalOrder();
      return *this;
    }

    // convert => CCoorMatrix
    /*! \brief
     *  Copy the sparse matrix \c M to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    CCoorMatrix<T> & operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return *this; }

    /*! \brief
     *  Initialize and empty sparse compressed coordinate matrix of \c nr rows and \c nc columns
     *  with reserved room for \c mnnz nonzeros.
     */
    void
    resize( indexType nr, indexType nc, indexType mnnz = SPARSETOOL_DEFAULT_NNZ ) {
      SPARSE::setup(nr, nc);
      I.resize(0); I.reserve(mnnz);
      J.resize(0); J.reserve(mnnz);
      A.resize(1); A.reserve(mnnz+1);
      A.back() = T(0);
      C.resize( SPARSE::sp_ncols + 1 );
    }

    /*! \brief
     *  Initialize the sparse compressed coordinate matrix and insert
     *  the elements M(i,j) which satify \c cmp(i,j) to \c *this.
     *  \param M   object derived from \c Sparse
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    void
    resize( SparseBase<MAT> const & M, Compare cmp )
    { convert(M,cmp); }

    /*! \brief
     *  Initialize the sparse compressed coordinate matrix and insert the element of the matrix \c M
     *  to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    void
    resize( SparseBase<MAT> const & M )
    { convert(M,all_ok()); }

    //! check if pattern is empty
    void empty(void) { SPARSE::sp_nnz = 0; }

    /*! \brief
     *  Order the internal structure to permit random access to nonzeros.
     */
    void
    internalOrder() {
      QuickSortIJ<indexType,T>( &I.front(), &J.front(), &A.front(), SPARSE::sp_nnz);
      // eliminate duplicate elements
      indexType i1 = 0, i = 0;
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
      SPARSE::sp_nnz       = i;
      SPARSE::sp_isOrdered = true;
      I.resize(SPARSE::sp_nnz);
      J.resize(SPARSE::sp_nnz);
      A.resize(SPARSE::sp_nnz+1);
      A.back() = T(0);

      indexType nc = 0;
      C(0) = 0;
      for ( indexType k = 1; k < SPARSE::sp_nnz; ++k ) {
        SPARSETOOL_TEST(
          J(k-1) <= J(k),
          "CCoorMatrix::internalOrderRowMajor() internal error J(" << k-1 <<
          ") is not <= J(" << k << ")"
        )
        indexType kk = J(k) - J(k-1);
        while ( kk-- > 0 ) C(++nc) = k;
      }
      SPARSETOOL_TEST(
        nc < SPARSE::sp_ncols || SPARSE::sp_nnz == 0,
        "CCoorMatrix::internalOrder() nc = " << nc <<
        " is not less of sp_ncols = " << SPARSE::sp_ncols <<
        " sp_nnz = " << SPARSE::sp_nnz
      )
      while ( ++nc <= SPARSE::sp_ncols ) C(nc) = SPARSE::sp_nnz;
    }

    //! Return the position of the element \c (i,j) in the vector storing elements
    indexType
    position( indexType i, indexType j ) const {
      SPARSE::test_index(i,j);
      indexType lo  = C(j);
      indexType hi  = C(j+1);
      indexType len = hi - lo;

      while ( len > 0 ) {
        indexType half = len / 2;
        indexType mid  = lo + half;
        if ( I(mid) < i ) {
          lo = mid + 1;
          len -= half + 1;
        } else len = half;
      }
      if ( lo == hi || I(lo) != i ) lo = SPARSE::sp_nnz;
      return lo;
    }

    //! Insert \c (i,j) in the sparse pattern
    valueType &
    insert( indexType i, indexType j ) {
      SPARSE::test_index(i,j);
      I.push_back(i);
      J.push_back(j);
      A.push_back(T(0));
      SPARSE::sp_nnz++;
      SPARSE::sp_isOrdered = false;
      return A(SPARSE::sp_nnz-1);
    }

    valueType const &
    operator [] (indexType idx) const {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    valueType &
    operator [] (indexType idx) {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    void setZero() { A = valueType(0); }
    void scaleValues( valueType const & s ) { A *= s; }

    valueType const &
    value( indexType i, indexType j ) const
    { return A(position(i,j)); }

    valueType const & 
    operator ()( indexType i, indexType j ) const {
      indexType pos = position(i,j);
      SPARSETOOL_TEST(
        pos != SPARSE::sp_nnz,
        "CCoorMatrix(" << i << "," << j <<
        ") referring to a non existent element"
      )
      return A(pos);
    }

    valueType & 
    operator ()( indexType i, indexType j ) {
      indexType pos = position(i,j);
      SPARSETOOL_TEST(
        pos != SPARSE::sp_nnz,
        "CCoorMatrix(" << i << "," << j <<
        ") referring to a non existent element"
      )
      return A(pos);
    }


    bool
    exists( indexType i, indexType j ) {
      return position(i,j) != SPARSE::sp_nnz;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<valueType> const & getA(void) const { return A; } //!< return the value vector
    Vector<valueType>       & getA(void)       { return A; } //!< return the value vector
    Vector<indexType> const & getI(void) const { return I; } //!< return the row index vector
    Vector<indexType> const & getJ(void) const { return J; } //!< return the column index vector
    Vector<indexType> const & getC(void) const { return C; } //!< return the pointer of the column index (valid after ordering)

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // ITERATOR
    inline void Begin (void) const { ipos = 0; }
    inline void Next  (void) const { ++ipos; }
    inline bool End   (void) const { return ipos < SPARSE::sp_nnz; }

    indexType         row    (void) const { return I(ipos); }
    indexType         column (void) const { return J(ipos); }
    valueType const & value  (void) const { return A(ipos); }
    valueType       & value  (void)       { return A(ipos); }

    //! Assign the pointed element to \c rhs
    template <typename TS>
    void assign( TS & rhs ) const { rhs = A(ipos); }

    //! perform the operation res += s * (A * x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_M_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      SPARSETOOL_TEST(
        (void*)&res != (void*)&x,
        "add_S_mul_M_mul_V equal pointer"
      )
      SPARSETOOL_TEST(
        res.size() >= SPARSE::sp_nrows,
        "result vector too small"
      )
      indexType const * pI = & I.front();
      indexType const * pJ = & J.front();
      valueType const * pA = & A.front();
      SPARSELIB_LOOP( SPARSE::sp_nnz, res(*pI++) += s * *pA++ * x(*pJ++) );
    }

    //! perform the operation res += s * (A ^ x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_Mt_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      SPARSETOOL_TEST(
        (void*)&res != (void*)&x,
        "add_S_mul_Mt_mul_V equal pointer"
      )
      SPARSETOOL_TEST(
        res.size() >= SPARSE::sp_ncols,
        "result vector too small"
      )
      indexType const * pI = & I.front();
      indexType const * pJ = & J.front();
      valueType const * pA = & A.front();
      SPARSELIB_LOOP( SPARSE::sp_nnz, res(*pJ++) += s * *pA++ * x(*pI++) );
    }
  };

  /*! \cond NODOC */
  SPARSELIB_MUL_STRUCTURES(CCoorMatrix)
  /*! \endcond */
  
  /*
  //  #     # #     # #       ####### ### ######  #       #     # 
  //  ##   ## #     # #          #     #  #     # #        #   #  
  //  # # # # #     # #          #     #  #     # #         # #   
  //  #  #  # #     # #          #     #  ######  #          #    
  //  #     # #     # #          #     #  #       #          #    
  //  #     # #     # #          #     #  #       #          #    
  //  #     #  #####  #######    #    ### #       #######    #    
  */
  /*!
   * Perform matrix multiplication of a matrix in Compressed Coordinate
   * with a full vector \c x
   * \code res += s * (A * x) \endcode
   * \param numRows numer of rows
   * \param RR      vector of row pointers
   * \param JJ      vector of column indexes
   * \param AA      vector of values
   * \param s       scalar
   * \param res     vector or the result
   * \param x       vector used in multiplication
   */
  template <typename T, typename VR, typename VX> inline
  void
  S_mul_M_mul_V(
    indexType                 numRows,
    Vector<indexType> const & RR,
    Vector<indexType> const & JJ,
    Vector<T>         const & AA,
    T const &                 s,
    VectorBase<T,VR>        & res,
    VectorBase<T,VX>  const & x
  ) {
    SPARSETOOL_TEST(
      (void*)&res != (void*)&x,
      "M_mul_V equal pointer"
    )
    SPARSETOOL_TEST(
      res.size() >= numRows,
      "result vector too small"
    )
    indexType const * R = & RR.front();
    indexType const * J = & JJ.front();
    T const *         A = & AA.front();
    for ( indexType ir = 0; ir < numRows; ++ir, ++R ) {
      T bf(0);
      SPARSELIB_LOOP( R[1] - R[0], bf += *A++ * x(*J++) );
      res(ir) += s * bf;
    }
  }

  /*!
   * Perform matrix multiplication of a transopose of a matrix in Compressed Coordinate
   * with a full vector \c x
   * \code res += s * (A ^ x) \endcode
   * \param numRows numer of rows
   * \param RR      vector of row pointers
   * \param JJ      vector of column indexes
   * \param AA      vector of values
   * \param s       scalar
   * \param res     vector or the result
   * \param x       vector used in multiplication
   */
  template <typename T, typename VR, typename VX> inline
  void
  S_mul_Mt_mul_V(
    indexType                 numRows,
    Vector<indexType> const & RR,
    Vector<indexType> const & JJ,
    Vector<T>         const & AA,
    T const &                 s,
    VectorBase<T,VR>        & res,
    VectorBase<T,VX>  const & x
  ) {
    SPARSETOOL_TEST(
      (void*)&res != (void*)&x,
      "M_mul_V equal pointer"
    )
    SPARSETOOL_TEST(
      res.size() >= numRows,
      "result vector too small"
    )
    indexType const * R = & RR.front();
    indexType const * J = & JJ.front();
    T const *         A = & AA.front();
    for ( indexType ir = 0; ir < numRows; ++ir, ++R ) {
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
  
  /*!
  
    The class \c CRowMatrix\<T\> implement a <B> Compressed Rows </B>
    storage sparse scheme.  It consists of two big vector of unsigned
    integer which contain the coordinate of nonzero elements and a big one
    of real number which contain the values.  We call \c R the vector
    that store the start position of each row, \c J the vector that
    store the column coordinate and \c A the vector that store the
    nonzero values.  For example the sparse matrix

    \htmlonly <TABLE><TR><TD> \endhtmlonly
    \f[
    \bm{A} = \left[\begin{BMAT}(e){c.c.c.c.c.c.c}{c.c.c.c.c.c}
     1 & 2    &    &   &   &   & 9 \\ 
       & -1   & 0  &   &   &   &   \\ 
       &      &    & 3 & 4 &   &   \\ 
       & 2    &    &   & 5 &   &   \\ 
     2 &      & -2 &   &   & 1 &   \\ 
       & -1.5 &    &   &   &   & -1 
   \end{BMAT}\right]
   \f]
   \htmlonly </TD><TD>&nbsp;&nbsp;&nbsp;</TD><TD> \endhtmlonly
   \f[
   \begin{tabular}{c|cccccccccccccc}
      R & 0 & 3 & 5 & 7 & 9 & 12 \\
      \hline
      J & 0 & 1 & 6 & 1 & 2 & 3 & 4 & 1 & 4 & 0 & 2 & 5 & 1 & 6 \\
      \hline
      A & 1 & 2 & 9 &-1 & 0 & 3 & 4 & 2 & 5 & 2 &-2 & 1 &-1.5 & -1
    \end{tabular}
    \f]
    \htmlonly </TD></TR></TABLE> \endhtmlonly

    To define a \c CCoorMatrix\<T\> class you can use one of the following
    scripture

\code
CRowMatrix<double> crow;         // instancean empty CRowMatrix<double> class.
CRowMatrix<double> crow(sp);     // instance a CRowMatrix<double> class with the sparsity pattern
                                  // defined in the sp SparsePattern class.
CRowMatrix<double> crow(crow1);  // instance a CRowMatrix<double> class which is the copy
                                  // of the CRowMatrix<double> object crow1.
CRowMatrix<double> crow(sobj);   // instance a CRowMatrix<double> class which is the copy
                                  // of the Sparse object crow1.
\endcode

    It is possible in any moment to change the sizes and the maximum
    number of nonzero by the \c resize method:

\code
  crow.resize(sp);
  crow.resize(crow1);
  crow.resize(sobj);
\endcode

  */
  //! Compressed Row Matrix Storage

  template <typename T>
  class CRowMatrix : public Sparse<T,CRowMatrix<T> > {
    typedef CRowMatrix<T>    MATRIX;
    typedef Sparse<T,MATRIX> SPARSE;
  public:
    typedef T valueType; //!< the type of the elements of the matrix

  private:

    Vector<valueType> A;
    Vector<indexType> R;
    Vector<indexType> J;

    mutable indexType iter_row;
    mutable indexType iter_ptr;

    void
    internalOrder() {
      SPARSETOOL_ASSERT(
        R[SPARSE::sp_nrows] == SPARSE::sp_nnz,
        "CRowMatrix::internalOrder() bad data for matrix"
      )
      indexType ii, kk, rk, rk1;
      SPARSE::sp_lower_nnz = 0;
      SPARSE::sp_diag_nnz  = 0;
      SPARSE::sp_upper_nnz = 0;
      for ( ii = 0, rk = R(0); ii < SPARSE::sp_nrows; ++ii, rk = rk1 ) {
        rk1 = R(ii+1);
        if ( rk1 > rk ) { // skip empty rows
          QuickSortI<indexType,T>( &J(rk), &A(rk), rk1 - rk);
          // setup statistic
          for ( kk = rk; kk < rk1; ++kk ) SPARSE::ldu_count(ii,J(kk));
  #ifdef SPARSETOOL_DEBUG
          for ( kk = rk+1; kk < rk1; ++kk )
            SPARSETOOL_ASSERT(
              J(kk-1) < J(kk),
              "CRowMatrix::internalOrder() failed"
            )
  #endif
        }
      }
      SPARSE::sp_nnz = R(SPARSE::sp_nrows);
    }

    template <typename MAT, typename Compare>
    void
    convert(SparseBase<MAT> const & M, Compare cmp) {

      SPARSE::setup( M.numRows(), M.numCols() );

      // step 0: Count nonzero
      SPARSE::sp_nnz = 0;
      for ( M.Begin(); M.End();  M.Next() )
        if ( cmp(M.row(), M.column() ) ) ++SPARSE::sp_nnz;

      A.resize( SPARSE::sp_nnz + 1 );
      J.resize( SPARSE::sp_nnz );
      R.resize( SPARSE::sp_nrows + 1 );
      A(SPARSE::sp_nnz) = 0;
      R = 0;

      // step 1: Evaluate not zero pattern
      for ( M.Begin(); M.End();  M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i, j) ) ++R(i);
      }
      for ( indexType k = 0; k < SPARSE::sp_nrows; ++k ) R(k+1) += R(k);

      // step 2: Fill matrix
      for ( M.Begin(); M.End(); M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) {
          indexType ii = --R(i);
          M.assign(A(ii));
          J(ii) = j;
        }
      }

      SPARSETOOL_TEST(
        R(0) == 0 && R(SPARSE::sp_nrows) == SPARSE::sp_nnz,
        "CRowMatrix::resize(Sparse & A) failed"
      )
      // step 3: internalOrder matrix
      internalOrder();
    }

  public:

    /*! \brief 
     *  Initialize and empty sparse compressed row matrix
     *  of \c 0 rows and \c 0 columns.
     */
    CRowMatrix(void) {};

    /*! \brief
     *  Insert the element of the sparse compressed row matrix \c M
     *  which satify \c cmp to \c *this.
     *  \param M   sparse compressed row matrix to be copied
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename Compare>
    CRowMatrix(CRowMatrix<T> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Copy the sparse compressed row matrix \c M to \c *this.
     *  \param M sparse compressed row matrix to be copied
     */
    CRowMatrix(CRowMatrix<T> const & M)
    { convert(M,all_ok()); }

    /*! \brief
     *  Insert the element of the sparse compressed row matrix \c M
     *  which satify \c cmp to \c *this.
     *  \param M   object derived from \c Sparse
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    CRowMatrix(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Copy the sparse compressed row matrix \c M to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    CRowMatrix(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    /*! \brief
     *  Copy the sparse compressed row matrix \c M to \c *this.
     *  \param M sparse compressed row matrix
     */
    CRowMatrix<T> &
    operator = (CRowMatrix<T> const & M) {
      if ( &M == this ) return *this; // avoid copy to itself
      SPARSE::setup( M.numRows(), M.numCols() );
      SPARSE::sp_nnz       = M.nnz();
      SPARSE::sp_isOrdered = M.isOrdered();
      A.load( M.A );
      R.load( M.R );
      J.load( M.J );
      if ( !SPARSE::sp_isOrdered ) internalOrder();
      return *this;
    }

    // convert => CRowMatrix
    /*! \brief
     *  Copy the sparse matrix \c M to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    CRowMatrix<T> & operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return *this; }

    /*! \brief
     *  Initialize the sparse compressed row matrix and insert
     *  the elements M(i,j) which satify \c cmp(i,j) to \c *this.
     *  \param M   object derived from \c Sparse
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    void
    resize(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Initialize the sparse compressed row matrix and insert the element of the matrix \c M
     *  to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    void
    resize(SparseBase<MAT> const & M)
    { resize(M,all_ok()); }

    void
    scaleRow( indexType nr, valueType const & val ) {
      SPARSE::test_row(nr);
      valueType * pA = &A.front() + R(nr);
      valueType * pB = &A.front() + R(nr+1);
      while ( pA < pB ) *pA++ *= val;
    }

    void
    scaleColumn(indexType nc, valueType const & val) {
      SPARSE::test_col(nc);
      for ( indexType i = 0; i < SPARSE::sp_nrows; ++i ) {
        indexType pos = position(i,nc,true);
        if ( pos != SPARSE::sp_nnz ) A(pos) *= val;
      }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<valueType> const & getA(void) const { return A; } //!< return the value vector
    Vector<valueType>       & getA(void)       { return A; } //!< return the value vector
    Vector<indexType> const & getR(void) const { return R; } //!< return the row pointer
    Vector<indexType> const & getJ(void) const { return J; } //!< return the column index vector

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    indexType
    position( indexType i, indexType j ) const {
      SPARSE::test_index(i,j);
      indexType lo  = R(i);
      indexType hi  = R(i+1);
      indexType len = hi - lo;
      while (len > 0) {
        indexType half = len / 2;
        indexType mid = lo + half;
        if ( J(mid) < j ) {
          lo = mid + 1;
          len -= half + 1;
        } else len = half;
      }
      if ( lo == hi || J(lo) != j ) lo = SPARSE::sp_nnz;
      return lo;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    void setZero() { A = valueType(0); }
    void scaleValues( valueType const & s ) { A *= s; }

    valueType const & 
    operator ()( indexType i, indexType j ) const {
      indexType pos = position(i,j);
      SPARSETOOL_TEST(
        pos != SPARSE::sp_nnz,
        "CRowMatrix(" << i << "," << j <<
        ") referring to a non existent element"
      )
      return A(pos);
    }

    valueType & 
    operator ()( indexType i, indexType j ) {
      indexType pos = position(i,j);
      SPARSETOOL_TEST(
        pos != SPARSE::sp_nnz,
        "CRowMatrix(" << i << "," << j <<
        ") referring to a non existent element"
      )
      return A(pos);
    }

    valueType const &
    value( indexType i, indexType j) const
    { return A(position(i,j)); }

    bool
    exists( indexType i, indexType j ) {
      return position(i,j) != SPARSE::sp_nnz;
    }

    valueType const &
    operator [] (indexType idx) const
    { SPARSE::test_nnz(idx); return A(idx); }

    valueType &
    operator [] (indexType idx)
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

    indexType         row    (void) const { return iter_row-1; }
    indexType         column (void) const { return J(iter_ptr); }
    valueType const & value  (void) const { return A(iter_ptr); }
    valueType       & value  (void)       { return A(iter_ptr); }

    //! Assign the pointed element to \c rhs
    template <typename TS>
    void assign( TS & rhs ) const { rhs = A(iter_ptr); }

    //! perform the operation res += s * (A * x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_M_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      S_mul_M_mul_V( SPARSE::sp_nrows, R, J, A, s, res, x );
    }

    //! perform the operation res += s * (A ^ x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_Mt_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      S_mul_Mt_mul_V( SPARSE::sp_nrows, R, J, A, s, res, x );
    }

  };

  /*! \cond NODOC */
  SPARSELIB_MUL_STRUCTURES(CRowMatrix)
  /*! \endcond */

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
  
  /*!

    The class \c CColMatrix\<T\> implement a <B> Compressed Columns </B>
    storage sparse scheme.  It consists of two big vector of unsigned
    integer which contain the coordinate of nonzero elements and a big one
    of real number which contain the values.  We call \c I the vector
    that store the row index of each elements, \c C the vector that
    store the starting position of the columns and \c A the vector
    that store the nonzero values.  For example the sparse matrix

    \htmlonly <TABLE><TR><TD> \endhtmlonly
    \f[
    \bm{A} = \left[\begin{BMAT}(e){c.c.c.c.c.c.c}{c.c.c.c.c.c}
      1 & 2    &    &   &   &   & 9 \\ 
        & -1   & 0  &   &   &   &   \\ 
        &      &    & 3 & 4 &   &   \\ 
        & 2    &    &   & 5 &   &   \\ 
      2 &      & -2 &   &   & 1 &   \\ 
        & -1.5 &    &   &   &   & -1 
    \end{BMAT}\right]
    \f]
    \htmlonly </TD><TD>&nbsp;&nbsp;&nbsp;</TD><TD> \endhtmlonly
    \f[
    \begin{tabular}{c|cccccccccccccc}
      I & 0 & 4 & 0 & 1 & 3 & 5 & 1 & 4 & 2 & 2 & 3 & 4 & 0 & 5 \\
      \hline
      C & 0 & 2 & 6 & 8 & 9 & 11 & 12 \\
      \hline
      A & 1 & 2 & 2 & -1 & 2 & -1.5 & 0 & -2 & 3 & 4 & 5 & 1 & 9 & -1
    \end{tabular}
    \f]
    \htmlonly </TD></TR></TABLE> \endhtmlonly

    Notice that the index follows the C convention, starting from \b 0. 
    The class \c CColMatrix\<T\> try to manage such a structure in a
    simple way for the user.  To define a \c CColMatrix\<T\> class you
    can use one of the following scripture

\code
CColMatrix<double> ccol;        // instance an empty CColMatrix<double> class.
CColMatrix<double> ccol(sp);    // instance a CColMatrix<double> class with the sparsity
                                 // pattern defined in the SparsePattern class sp.  
CColMatrix<double> ccol(ccol1); // instance a CColMatrix<double> class which is the copy
                                 // of the CColMatrix<double> object ccol1.
CColMatrix<double> ccol(sobj);  // instance a CColMatrix<double> class which is the copy
                                 // of the Sparse object sobj.
\endncode

    It is possible in any moment to change the sizes and the maximum
    number of nonzero by the \c resize method:

\code
  ccol.resize(sp);
  ccol.resize(ccol1);
  ccol.resize(sobj);
\endcode

  */
  //! Compressed Column Matrix Storage
  template <typename T>
  class CColMatrix : public Sparse<T,CColMatrix<T> > {
    typedef CColMatrix<T>    MATRIX;
    typedef Sparse<T,MATRIX> SPARSE;
  public:
    typedef T valueType; //!< the type of the elements of the matrix

  private:
    Vector<valueType> A;
    Vector<indexType> C, I;

    mutable indexType iter_col;
    mutable indexType iter_ptr;

    void
    internalOrder() {
      indexType jj, kk, ck, ck1;
      SPARSE::sp_lower_nnz = 0;
      SPARSE::sp_diag_nnz  = 0;
      SPARSE::sp_upper_nnz = 0;
      for ( jj = 0, ck = C(0); jj < SPARSE::sp_ncols; ++jj, ck = ck1 ) {
        ck1 = C(jj+1);
        if ( ck1 > ck ) { // skip empty columns
          QuickSortI<indexType,T>( &I(ck), &A(ck), ck1 - ck );
          // setup statistic
          for ( kk = ck; kk < ck1; ++kk ) SPARSE::ldu_count(I(kk),jj);
  #ifdef SPARSETOOL_DEBUG
          for ( kk = ck+1; kk < ck1; ++kk )
            SPARSETOOL_ASSERT(
              I(kk-1) < I(kk),
              "CColMatrix::internalOrder() failed"
            )
  #endif
        }
      }
      SPARSE::sp_nnz = C(SPARSE::sp_ncols);
    }

    template <typename MAT, typename Compare>
    void
    convert(SparseBase<MAT> const & M, Compare cmp) {
      SPARSE::setup( M.numRows(), M.numCols() );

      // step 0: Count nonzero
      SPARSE::sp_nnz = 0;
      for ( M.Begin(); M.End();  M.Next() )
        if ( cmp(M.row(), M.column() ) ) ++SPARSE::sp_nnz;

      A.resize( SPARSE::sp_nnz + 1 );
      I.resize( SPARSE::sp_nnz );
      C.resize( SPARSE::sp_ncols + 1 );
      C = 0;
      A(SPARSE::sp_nnz) = 0;

      // step 1: Evaluate not zero
      for ( M.Begin(); M.End(); M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) ++C(j);
      }
      for ( indexType k = 0; k < SPARSE::sp_ncols; ++k ) C(k+1) += C(k);

      // step 2: Fill matrix
      for ( M.Begin(); M.End(); M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) {
          indexType jj = --C(j);
          M.assign(A(jj));
          I(jj) = i;
        }
      }
      SPARSETOOL_TEST(
        C(0) == 0 && C(SPARSE::sp_ncols) == SPARSE::sp_nnz,
        "CColMatrix::resize(Sparse & A) failed"
      )
      // step 3: internalOrder matrix
      internalOrder();
    }

  public:

    /*! \brief 
     *  Initialize and empty sparse compressed column matrix
     *  of \c 0 rows and \c 0 columns.
     */
    CColMatrix(void) {};

    /*! \brief
     *  Insert the element of the sparse compressed column matrix \c M
     *  which satify \c cmp to \c *this.
     *  \param M   sparse compressed row matrix to be copied
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename Compare>
    CColMatrix(CColMatrix<T> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Copy the sparse compressed row matrix \c M to \c *this.
     *  \param M sparse compressed row matrix to be copied
     */
    CColMatrix(CColMatrix<T> const & M)
    { convert(M,all_ok()); }

    /*! \brief
     *  Insert the element of the sparse compressed column matrix \c M
     *  which satify \c cmp to \c *this.
     *  \param M   object derived from \c Sparse
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    CColMatrix(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Copy the sparse compressed column matrix \c M to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    CColMatrix(SparseBase<MAT> const & M)
    { convert(M,all_ok()); }

    /*! \brief
     *  Copy the sparse compressed column matrix \c M to \c *this.
     *  \param M sparse compressed column matrix
     */
    CColMatrix<T> &
    operator = (CColMatrix<T> const & M) {
      if ( &M == this ) return *this; // avoid copy to itself
      resize( M.numRows(), M.numCols(), M.nnz() );
      SPARSE::sp_nnz       = M.nnz();
      SPARSE::sp_isOrdered = M.isOrdered();
      A.load( M.A );
      C.load( M.C );
      I.load( M.I );
      if ( !SPARSE::sp_isOrdered ) internalOrder();
      return *this;
    }

    // convert => CColMatrix
    /*! \brief
     *  Copy the sparse matrix \c M to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    CColMatrix<T> &
    operator = (SparseBase<MAT> const & M)
    { convert(M,all_ok()); return * this; }

    /*! \brief
     *  Initialize the sparse compressed column matrix and insert
     *  the elements M(i,j) which satify \c cmp(i,j) to \c *this.
     *  \param M   object derived from \c Sparse
     *  \param cmp comparator struct for the selection of the element in pattern \c sp
     */
    template <typename MAT, typename Compare>
    void
    resize(SparseBase<MAT> const & M, Compare cmp)
    { convert(M,cmp); }

    /*! \brief
     *  Initialize the sparse compressed column matrix and insert the element of the matrix \c M
     *  to \c *this.
     *  \param M object derived from \c Sparse
     */
    template <typename MAT>
    void
    resize( SparseBase<MAT> const & M )
    { resize(M,all_ok()); }

    void
    scaleRow( indexType nr, valueType const & val ) {
      SPARSE::test_row(nr);
      for ( indexType j = 0; j < SPARSE::sp_ncols; ++j ) {
        indexType pos = position(nr,j,true);
        if ( pos != SPARSE::sp_nnz ) A(pos) *= val;
      }
    }

    void
    scaleColumn( indexType nc, valueType const & val ) {
      SPARSE::test_col(nc);
      valueType * pA = &A.front() + C(nc);
      valueType * pB = &A.front() + C(nc+1);
      while ( pA < pB ) *pA++ *= val;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<valueType> const & getA(void) const { return A; } //!< return the value vector
    Vector<valueType>       & getA(void)       { return A; } //!< return the value vector
    Vector<indexType> const & getI(void) const { return I; } //!< return the row index vector
    Vector<indexType> const & getC(void) const { return C; } //!< return the column pointer vector

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    indexType
    position( indexType i, indexType j ) const {
      SPARSE::test_index(i,j);
      indexType lo  = C(j);
      indexType hi  = C(j+1);
      indexType len = hi - lo;
      while (len > 0) {
        indexType half = len / 2;
        indexType mid = lo + half;
        if ( I(mid) < i ) {
          lo = mid + 1;
          len -= half + 1;
        } else len = half;
      }
      if ( lo == hi || I(lo) != i ) lo = SPARSE::sp_nnz;
      return lo;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    void setZero() { A = valueType(0); }
    void scaleValues( valueType const & s ) { A *= s; }

    valueType const &
    operator () (indexType i, indexType j) const {
      indexType pos = position(i,j);
      SPARSETOOL_TEST(
        pos != SPARSE::sp_nnz,
        "CColMatrix(" << i << "," << j <<
        ") referring to a non existent element"
      )
      return A(pos);
    }

    valueType &
    operator () (indexType i, indexType j) {
      indexType pos = position(i,j);
      SPARSETOOL_TEST(
        pos != SPARSE::sp_nnz,
        "CColMatrix(" << i << "," << j <<
        ") referring to a non existent element"
      )
      return A(pos);
    }

    valueType const &
    value(indexType i, indexType j) const
    { return A(position(i,j)); }

    bool
    exists( indexType i, indexType j ) {
      return position(i,j) != SPARSE::sp_nnz;
    }

    valueType const &
    operator [] (indexType idx) const
    { SPARSE::test_nnz(idx); return A(idx); }

    valueType &
    operator [] (indexType idx)
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

    indexType         row    (void) const { return I(iter_ptr); }
    indexType         column (void) const { return iter_col-1;  }
    valueType const & value  (void) const { return A(iter_ptr); }
    valueType       & value  (void)       { return A(iter_ptr); }

    //! Assign the pointed element to \c rhs
    template <typename TS>
    void assign( TS & rhs ) const { rhs = A(iter_ptr); }

    //! perform the operation res += s * (A * x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_M_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      S_mul_Mt_mul_V( SPARSE::sp_ncols, C, I, A, s, res, x );
    }

    //! perform the operation res += s * (A ^ x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_Mt_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      S_mul_M_mul_V( SPARSE::sp_ncols, C, I, A, s, res, x );
    }

  };

  /*! \cond NODOC */
  SPARSELIB_MUL_STRUCTURES(CColMatrix)
  /*! \endcond */

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

  /*!
  
     The class <c> TridMatrix<T> </c> implement a sparse band matrix.
     It consists of a big matrix of \c valueType which contain the values
     of nonzero elements.  We call \c M this big matrix which
     represents an \f$ n\times m \f$ matrix which stores the rows of nonzero. 
     For example the following band matrix
     
     \htmlonly <TABLE><TR><TD> \endhtmlonly
     \f[
       \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c.c.c.c}{c.c.c.c.c.c.c}
         2 & -1 &    &    &    &    &    \\ 
         1 & 2  & -2 &    &    &    &    \\ 
           & 2  & 2  & -3 &    &    &    \\ 
           &    & 3  & 2  & -4 &    &    \\ 
           &    &    & 4  & 2  & -5 &    \\ 
           &    &    &    & 5  & 2  & -6 \\ 
           &    &    &    &    & 6  & 2 
        \end{BMAT}\right]
        \qquad
        \begin{tabular}{c|c|c}
           L & D & U \\
          \hline
          * & 2 & -1 \\
          1 & 2 & -2 \\
          2 & 2 & -3 \\
          3 & 2 & -4 \\ 
          4 & 2 & -5 \\ 
          5 & 2 & -6 \\ 
          6 & 2 & *  \\
          \hline
        \end{tabular}
     \f]
     \htmlonly </TD></TR></TABLE> \endhtmlonly
    
     where * means unused elements.  Notice that the index follows
     the C convention, starting from \c 0.
     The class <c> TridMatrix<T> </c> manage such a structure in a simple
     way for the user.  To define a <c> TridMatrix<T> </c> class you can
     use one of the following scripture

\code
  TridMatrix<double> tm;      // instance an empty TridMatrix<double> class
  TridMatrix<double> tm(100); // instance a 100 x 100 TridMatrix<double> class
  TridMatrix<double> tm(tm1); // instance a TridMatrix<double> class which is the copy
                              // of the TridMatrix<double> class tm1.
\endcode

    It is possible in any moment to change the sizes and the maximum
    number of nonzero by the \c resize method:
\code
  tm.resize(n);
\endcode

  */
  //! Tridiagonal Matrix

  template <typename T>
  class TridMatrix : public Sparse<T,TridMatrix<T> > {
    typedef TridMatrix<T>    MATRIX;
    typedef Sparse<T,MATRIX> SPARSE;
  public:
    typedef T valueType; //!< type of the element of the matrix

  private:

    template <typename VECTOR>
    void
    solve(VECTOR const & b, VECTOR & x) const {

      VECTOR ll(SPARSE::sp_nrows-1), dd(SPARSE::sp_nrows);

      // LU DEC
      indexType k = 0;
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

    Vector<T>      A;
    VectorSlice<T> L, D, U;

    mutable indexType iter_counter, iter_row, iter_col;

  public:

    //! contruct and exmpty \c 0x0 tridiagonal matrix
    TridMatrix( void ) {}

    //! contruct a <c> n x n </c> tridiagonal matrix
    TridMatrix( indexType ns ) { resize(ns); }

    //! contruct a copy of the tridiagonal matrix \c TM
    TridMatrix( TridMatrix<T> const & TM ) { resize(TM); }

    //! contruct a copy of the tridiagonal matrix \c TM
    TridMatrix<T> const &
    operator = ( TridMatrix<T> const & TM ) {
      if ( &TM == this ) return *this; // avoid copy to itself
      resize( A.nnz() );
      A = TM.A;
      return *this;
    }

    //! resize to a <c> n x n </c> tridiagonal matrix
    void
    resize( indexType ns ) {
      A.resize(3*ns-2);
      indexType kk = 0;
      L.slice(A, kk, kk + ns-1 ); kk += ns-1;
      D.slice(A, kk, kk + ns   ); kk += ns;
      U.slice(A, kk, kk + ns-1 );
      SPARSE::setup(ns,ns);
      SPARSE::sp_nnz = 3*ns-2;
    }

    valueType const &
    operator [] (indexType idx) const {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    valueType &
    operator [] (indexType idx) {
      SPARSE::test_nnz(idx);
      return A(idx);
    }

    void
    scaleRow( indexType nr, valueType const & val ) {
      SPARSE::test_row(nr);
      D(nr) *= val;
      if ( nr < SPARSE::sp_now ) U(nr)   *= val;
      if ( nr > 0              ) L(nr-1) *= val;
    }

    void
    scaleColumn( indexType nc, valueType const & val ) {
      SPARSE::test_col(nc);
      D(nc) *= val;
      if ( nc > 0               ) U(nc-1) *= val;
      if ( nc < SPARSE::sp_nrow ) L(nc)   *= val;
    }

    /*! \name Diagonal internal operation
     */
    //@{
    /*! \brief
        Set \c s to all the components of the diagonal,
        and \c 0 the others elements.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        (\bm{A}=2.1) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2.1 & 0    &      &     \\ 
            0   & 2.1  & 0    &     \\ 
                & 0    & 2.1  & 0   \\ 
                &      & 0    & 2.1 \\
         \end{BMAT}\right]
      \f]
      \htmlonly </TD></TR></TABLE> \endhtmlonly
    */

    TridMatrix<T> &
    operator = ( valueType const & s )
    { L = 0; U = 0; D = s; return *this; }

    /*! \brief
        Add \c s to all the components of the diagonal,
        and \c 0 the others elements.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        (\bm{A}+=2.1) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            4.1 & -1   &      &    \\ 
            1   & 4.1  & -1   &    \\ 
                & 1    & 4.1  & -1 \\ 
                &      & 1    & 4.1  \\
         \end{BMAT}\right]
        \f]    
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    TridMatrix<T> &
    operator += ( valueType const & s )
    { D += s; return *this; }

    /*! \brief
        Subtract \c s to all the components of the diagonal,
        and \c 0 the others elements.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        (\bm{A}-=2.1) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            -0.1 & -1   &      &    \\ 
            1   & -0.1  & -1   &    \\ 
                & 1    & -0.1  & -1 \\ 
                &      & 1    & -0.1  \\
         \end{BMAT}\right]
        \f]    
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    TridMatrix<T> &
    operator -= ( valueType const & s )
    { D -= s; return *this; }
    
    /*! \brief
        Multiply by \c s to all the nonzeros.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        (\bm{A}*=2) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            4 & -2 &    &    \\ 
            2 & 4  & -2 &    \\ 
              & 2  & 4  & -2 \\ 
              &    & 2  & 4  \\
         \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    TridMatrix<T> &
    operator *= ( valueType const & s )
    { A *= s; return *this; }
    
    /*! \brief
        Divide by \c s to all the nonzeros.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        (\bm{A}/=2) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            1   & -0.5 &      &    \\ 
            0.5 & 1    & -0.5 &    \\ 
                & 0.5  & 1    & -0.5 \\ 
                &      & 0.5  & 1  \\
         \end{BMAT}\right]
        \f]    
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    TridMatrix<T> &
    operator /= ( valueType const & s )
    { A /= s; return *this; }

    /*! \brief
        Set the diagonal values of the sparse
        matrix to \c v.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        \bm{v} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
        1 \\ 2\\ 3
        \end{BMAT}\right],
        \qquad
        (\bm{A}/=2) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            1 & 0  &   &   \\ 
            0 & 2  & 0 &   \\ 
              & 0  & 3 & 0 \\ 
              &    & 0 & 0 \\
         \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename VECTOR> inline
    TridMatrix<T> &
    operator = ( VectorBase<T,VECTOR> const & v )
    { L = 0; U = 0; D = v; return *this; }

    /*! \brief
        Add  \c v to the diagonal values of the sparse
        matrix.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        \bm{v} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
        1 \\ 2\\ 3
        \end{BMAT}\right],
        \qquad
        (\bm{A}\verb|+=|\bm{v}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            3 & -1 &    &    \\ 
            1 & 4  & -1 &    \\ 
              & 1  & 5  & -1 \\ 
              &    & 1  & 2  \\
         \end{BMAT}\right]
        \f]    
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename VECTOR> inline
    TridMatrix<T> &
    operator += ( VectorBase<T,VECTOR> const & v )
    { D += v; return *this; }

    /*! \brief
        Subtract  \c v to the diagonal values of the sparse
        matrix.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        \bm{v} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
        1 \\ 2\\ 3
        \end{BMAT}\right],
        \qquad
        (\bm{A}\verb|-=|\bm{v}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            3 & -1 &    &    \\ 
            1 & 4  & -1 &    \\ 
              & 1  & 5  & -1 \\ 
              &    & 1  & 2  \\
         \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename VECTOR> inline
    TridMatrix<T> &
    operator -= ( VectorBase<T,VECTOR> const & v )
    { D -= v; return *this; }

    /*! \brief
        Set the diagonal values of the sparse
        matrix to the vector expression \c e.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        \bm{a} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
        1 \\ 2\\ 3
        \end{BMAT}\right],\qquad
        \bm{b} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c.c}
        1 \\ -1 \\ 2 \\ 1
        \end{BMAT}\right],
        \qquad
        (\bm{A}=\bm{a}+2*\bm{b}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            3 & 0  &   &   \\ 
            0 & 0  & 0 &   \\ 
              & 0  & 7 & 0 \\ 
              &    & 0 & 0 \\
         \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename R> inline
    TridMatrix<T> &
    operator = ( VectorE<T,R> const & e )
    { L = 0; U = 0; D = e;; return *this; }

    /*! \brief
        Add to the diagonal values of the sparse
        matrix the vector expression \c e.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        \bm{a} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
        1 \\ 2\\ 3
        \end{BMAT}\right],\qquad
        \bm{b} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c.c}
        1 \\ -1 \\ 2 \\ 1
        \end{BMAT}\right],
        \qquad
        (\bm{A}\verb|+=|\bm{a}+2*\bm{b}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            5 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 9  & -1 \\ 
              &    & 1  & 2  \\
         \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename R> inline
    TridMatrix<T> &
    operator += ( VectorE<T,R> const & e )
    { D += e;; return *this; }

    /*! \brief
        Add to the diagonal values of the sparse
        matrix the vector expression \c e.  For example
        \htmlonly <TABLE><TR><TD> \endhtmlonly
        \f[
        \bm{A} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            2 & -1 &    &    \\ 
            1 & 2  & -1 &    \\ 
              & 1  & 2  & -1 \\ 
              &    & 1  & 2  \\
        \end{BMAT}\right],\qquad
        \bm{a} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c}
        1 \\ 2\\ 3
        \end{BMAT}\right],\qquad
        \bm{b} = \left[\begin{BMAT}(b,0.5cm,0.5cm){c}{c.c.c.c}
        1 \\ -1 \\ 2 \\ 1
        \end{BMAT}\right],
        \qquad
        (\bm{A}\verb|-=|\bm{a}+2*\bm{b}) = \left[\begin{BMAT}(b,0.5cm,0.5cm){c.c.c.c}{c.c.c.c}
            -1 & -1 &    &    \\ 
             1 & 2  & -1 &    \\ 
               & 1  & -5 & -1 \\ 
               &    & 1  & 2  \\
         \end{BMAT}\right]
        \f]
        \htmlonly </TD></TR></TABLE> \endhtmlonly
    */
    template <typename R> inline
    TridMatrix<T> &
    operator -= ( VectorE<T,R> const & e )
    { D -= e;; return *this; }

    //@}

    valueType const &
    operator () ( indexType i, indexType j ) const {
      SPARSE::test_index(i,j);
      indexType kk = ((j+1)-i); 
      SPARSETOOL_TEST(
        kk < 3,
        "TridMatrix(" << i << "," << j << ") index out of range"
      )
      return A( kk*SPARSE::sp_nrows+i-1);
    }

    valueType &
    operator () ( indexType i, indexType j ) {
      SPARSE::test_index(i,j);
      indexType kk = ((j+1)-i); 
      SPARSETOOL_TEST(
        kk < 3,
        "TridMatrix(" << i << "," << j << ") index out of range"
      )
      return A( kk*SPARSE::sp_nrows+i-1);
    }

    valueType const &
    value( indexType i, indexType j ) const {
      static const T zero(0);
      SPARSE::test_index(i,j);
      indexType kk = ((j+1)-i);
      if ( kk < 3 ) return A( kk*SPARSE::sp_nrows+i-1);
      else          return zero;
    }

    bool
    exists( indexType i, indexType j ) {
      return i < SPARSE::sp_nrows &&
             j < SPARSE::sp_nrows &&
             (i <= j+1 || j <= i+1 );
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Vector<valueType> const & getA(void) const { return A; } //!< return the value vector
    Vector<valueType>       & getA(void)       { return A; } //!< return the value vector
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    void
    Begin(void) const {
      iter_counter = 0;
      iter_row     = 0;
      iter_col     = 0;
    }

    void
    Next (void) const {
      ++iter_counter;
      iter_row = (iter_counter+1) % SPARSE::sp_nrows;
      iter_col = ( (iter_counter+1) / SPARSE::sp_nrows) +iter_row-1;
    }

    bool
    End(void) const
    { return iter_counter < SPARSE::sp_nnz; }

    indexType         row    (void) const { return iter_row; }
    indexType         column (void) const { return iter_col; }
    valueType const & value  (void) const { return A(iter_counter); }
    valueType       & value  (void)       { return A(iter_counter); }

    //! Assign the pointed element to \c rhs
    template <typename TS>
    void assign( TS & rhs ) const { rhs = A(iter_counter); }

    //! perform the operation res += s * (A * x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_M_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      SPARSETOOL_TEST(
        (void*)&res != (void*)&x,
        "add_S_mul_M_mul_V equal pointer"
      )
      SPARSETOOL_TEST(
        res.size() >= SPARSE::sp_nrows,
        "result vector too small"
      )
      res += s * D * x;
      res += s * U * x;
      SPARSELIB_LOOP( SPARSE::sp_nrows - 1, res(i+1) += s * L(i) * x(i+1) );
    }

    //! perform the operation res += s * (A ^ x)
    template <typename VRES, typename VB> inline
    void
    add_S_mul_Mt_mul_V(
      VectorBase<T,VRES>     & res,
      T const                & s,
      VectorBase<T,VB> const & x
    ) const {
      SPARSETOOL_TEST(
        (void*)&res != (void*)&x,
        "add_S_mul_Mt_mul_V equal pointer"
      )
      SPARSETOOL_TEST(
        res.size() >= SPARSE::sp_nrows,
        "result vector too small"
      )
      res += s * D * x;
      res += s * L * x;
      SPARSELIB_LOOP( SPARSE::sp_nrows - 1, res(i+1) += s * U(i) * x(i+1) );
    }

    //! solve the tridiagonal system \c A*x=b
    template <typename VECTOR> inline
    void
    ass_V_div_M(VECTOR & res, VECTOR const & b) const {
      solve(b, res);
    }

  };

  /*! \cond NODOC */
  SPARSELIB_MUL_STRUCTURES(TridMatrix)

  template <typename TM, typename T> inline
  Vector_V_div_M<Vector<T>,TridMatrix<TM> >
  operator / (Vector<T> const & a, TridMatrix<TM> const & M) {
    return Vector_V_div_M<Vector<T>,TridMatrix<TM> >(a,M);
  }
  /*! \endcond */

  /*
   *  ###       # #######
   *   #       #  #     #
   *   #      #   #     #
   *   #     #    #     #
   *   #    #     #     #
   *   #   #      #     #
   *  ### #       #######
   */

  //! \name I/O
  //@{

  //! Print the contents of vector \c v to the stream \c s
  template <typename T, typename VEC>
  inline
  void
  print( ostream & s, VectorBase<T,VEC> const & v ) {
    indexType sz1 = v.size() - 1;
    s << "[";
    for ( indexType i = 0; i < sz1; ++i ) s << v(i) << " , ";
    s << v(sz1) << "]";
  }

  //! Print the contents of vector \c v to the stream \c s
  template <typename T, typename VEC> inline
  ostream &
  operator << (ostream & s, VectorBase<T,VEC> const & v) {
    print(s,v);
    return s;
  }

  //! Print the contents of vector expression \c e to the stream \c s
  template <typename T, typename A>
  inline
  ostream &
  operator << (ostream & s, VectorE<T,A> const & e) {
    print(s,e);
    return s;
  }

  //! Read from the stream \c s to the vector \c v
  template <typename T> inline
  istream &
  operator >> (istream & s, Vector<T> & v) {
    typename Vector<T>::pointer pv=v.begin();
    while ( pv != v.end() ) s >> *pv++;
    return s;
  }

  //! Print the contents of the SparsePattern \c S to the stream \c s
  inline
  ostream &
  operator << (ostream & s, SparsePattern const & S) {
    for ( S.Begin(); S.End(); S.Next() )
      s << S.row() << " " << S.column() << "\n";
    return s;
  }

  //! Print the contents of the Sparse\<T,MAT\> object \c N to the stream \c s
  template <typename T, typename MAT> inline
  ostream &
  operator << (ostream & s, Sparse<T,MAT> const & M) {
    s << "\nsize    = " << M.numRows() << " x " << M.numCols()
      << "\nnnz     = " << M.nnz();
    if ( M.isOrdered() ) {
      indexType l,d,u;
      M.nnz( l, d, u );
      s << "  (diag=" << d << ", lower=" << l << ", upper=" << u << ")\nordered = YES\n";
    } else {
      s << "\nordered = NO\n";
    }
    indexType N = 1;
    indexType n = M.maxSize();
    while ( (n/=10) > 0 ) ++N; 
    for ( M.Begin(); M.End(); M.Next() ) {
      s << "("    << setw(N) << M.row()
        << ","    << setw(N) << M.column()
        << ") = " << setw(N) << M.value()
        << "\n";
    }
    return s;
  }
  //@}


}

namespace SparseToolLoad {

  using ::SparseTool::Vector;
  using ::SparseTool::VectorSlice;
  using ::SparseTool::SparsePattern;
  using ::SparseTool::CCoorMatrix;
  using ::SparseTool::CRowMatrix;
  using ::SparseTool::CColMatrix;
  using ::SparseTool::TridMatrix;
  
  using ::SparseTool::absval;
  using ::SparseTool::sin;
  using ::SparseTool::cos;
  using ::SparseTool::tan;
  using ::SparseTool::asin;
  using ::SparseTool::acos;
  using ::SparseTool::atan;
  using ::SparseTool::cosh;
  using ::SparseTool::sinh;
  using ::SparseTool::tanh;

  using ::SparseTool::sqrt;
  using ::SparseTool::ceil;
  using ::SparseTool::floor;
  using ::SparseTool::log;
  using ::SparseTool::log10;

  using ::SparseTool::pow;
  using ::SparseTool::atan2;

  using ::SparseTool::norm1;
  using ::SparseTool::norm2;
  using ::SparseTool::normi;
  using ::SparseTool::normp;
  using ::SparseTool::maxval;
  using ::SparseTool::minval;
  using ::SparseTool::maxabsval;
  using ::SparseTool::minabsval;
  using ::SparseTool::dot;
  using ::SparseTool::rdot;
  using ::SparseTool::sum;
  using ::SparseTool::prod;
  using ::SparseTool::dist;
  using ::SparseTool::dist2;

  // I/O
  using ::SparseTool::operator >>;
  using ::SparseTool::operator <<;

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
