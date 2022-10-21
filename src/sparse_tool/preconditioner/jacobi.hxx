#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_JACOBI_HH
#define SPARSETOOL_ITERATIVE_PRECO_JACOBI_HH

namespace Sparse_tool {

  /*
  //        #
  //        #   ##    ####   ####  #####  #
  //        #  #  #  #    # #    # #    # #
  //        # #    # #      #    # #####  #
  //  #     # ###### #      #    # #    # #
  //  #     # #    # #    # #    # #    # #
  //   #####  #    #  ####   ####  #####  #
  */
  //! Iterative JACOBI preconditioner
  template <typename T>
  class JACOBIpreconditioner : public Preco<JACOBIpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using JACOBIPRECO = JACOBIpreconditioner<T>;
    using PRECO       = Preco<JACOBIPRECO>;
    #endif

    using real_type = T; //!< type of the elements of the preconditioner

  private:

    real_type         omega, omega1;
    integer           maxIter;

    Vector<integer>   LU_R;
    Vector<integer>   LU_J;
    Vector<real_type> LU_A;

    Vector<real_type> D;

    mutable Vector<real_type> TMP, TMP1;

    Vector<integer> LUnnz;

    //! build incomplete LDU decomposition with specified pattern `P`
    template <typename MAT>
    void
    build_JACOBI(
      MAT       const & A,
      real_type const & _omega,
      integer           _maxIter
    ) {

      UTILS_ASSERT0(
        A.is_ordered(),
        "Sparse_tool: JACOBIpreconditioner::build_JACOBI\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        A.nrows() == A.ncols(),
        "Sparse_tool: ACOBIpreconditioner::build_JACOBI\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        A.nrows() > 0,
        "Sparse_tool: JACOBIpreconditioner::build_JACOBI\n"
        "empty matrix\n"
      );

      this -> omega   = _omega;
      this -> omega1  = 1.0/_omega-1.0;
      this -> maxIter = _maxIter;

      // step 0: compute necessary memory
      PRECO::pr_size = A.nrows();
      LUnnz.resize( PRECO::pr_size );
      D.resize( PRECO::pr_size );
      TMP.resize( PRECO::pr_size );
      TMP1.resize( PRECO::pr_size );

      LUnnz.setZero();
      D.setZero();

      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if ( i != j ) ++LUnnz(i);
      }

      // step 1: initialize structure
      LU_R.resize( PRECO::pr_size + 1 );
      LU_R(0) = 0;
      for ( integer i = 0; i < PRECO::pr_size; ++i ) LU_R(i+1) = LU_R(i) + LUnnz(i);

      // step 2: allocate memory
      LU_A.resize( LU_R(PRECO::pr_size) );
      LU_J.resize( LU_R(PRECO::pr_size) );
      LU_A.setZero();

      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        if ( i != j ) LU_J(LU_R(i)+(--LUnnz(i))) = j;
      }

      // step 4: sort structure
      for ( integer i = 0; i < PRECO::pr_size; ++i )
        std::sort( LU_J.data()+LU_R(i), LU_J.data()+LU_R(i+1) );

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        real_type val = A.value(); // (i,j);
        if ( i != j ) { // costruisco LU
          integer lo  = LU_R(i);
          integer hi  = LU_R(i+1);
          integer len = hi - lo;
          while ( len > 0 ) {
            integer half = len / 2;
            integer mid  = lo + half;
            if ( LU_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                   len = half;
          }
          LU_A(lo) = val;
        } else {
          D(i) = val;
        }
      }
      for ( integer i = 0; i < PRECO::pr_size; ++i )
        UTILS_ASSERT(
          D(i) != real_type(0),
          "Sparse_tool: JACOBIpreconditioner::D({}) = {}, size = {}\n",
          i, D(i), D.size()
        );
    }

  public:

    JACOBIpreconditioner(void) : Preco<JACOBIPRECO>() {}

    template <typename MAT>
    JACOBIpreconditioner( MAT const & M, real_type _omega, integer _maxIter ) : Preco<JACOBIPRECO>()
    { build_JACOBI( M, _omega, _maxIter ); }

    //! build the preconditioner from matrix `M` with pattern `P`
    template <typename MAT>
    void
    build( MAT const & M, real_type _omega, integer _maxIter )
    { build_JACOBI( M, _omega, _maxIter ); }

    //!
    //! Apply preconditioner to vector `b`
    //! and store result to vector `x`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & x, VECTOR const & b ) const {
      UTILS_ASSERT0(
        std::addressof(x) != std::addressof(b),
        "Sparse_tool: JACOBIpreconditioner::ass_preco(x,b)\n"
        "`x` and `b` cant be the same\n"
      );
      x = omega*(b.array()/D.array());
      for ( integer ii = 0; ii < maxIter; ++ii ) {

        // calcolo (D/omega)^{-1} [ ((1/omega-1)*D-L-U)*x + b ];
        integer   const * _pR = LU_R.data();
        integer   const * _pJ = LU_J.data();
        real_type const * _pA = LU_A.data();

        for ( integer k=0; k < PRECO::pr_size; ++k ) {
          real_type tt(0);
          for ( integer i_cnt = _pR[1] - _pR[0]; i_cnt > 0; --i_cnt )
            tt += *_pA++ * x(*_pJ++);
          TMP(k) = b(k) + omega1 * D(k) * x(k) - tt;
          ++_pR;
        }
        x = omega * (TMP.array()/D.array());
      }
    }

    //!
    //! Apply preconditioner to vector `x`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & x ) const {
      TMP1 = x;
      ass_preco( x, TMP1 );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,JACOBIpreconditioner<TP> >
  operator / (Vector<T> const & v, JACOBIpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,JACOBIpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::JACOBIpreconditioner;
}
#endif

#endif
