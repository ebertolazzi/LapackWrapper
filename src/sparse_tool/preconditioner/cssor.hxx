#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_CSSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_CSSOR_HH

namespace Sparse_tool {

  /*
  //   #####   #####   #####  ####### ######
  //  #     # #     # #     # #     # #     #
  //  #       #       #       #     # #     #
  //  #        #####   #####  #     # ######
  //  #             #       # #     # #   #
  //  #     # #     # #     # #     # #    #
  //   #####   #####   #####  ####### #     #
  */
  //! Iterative SSOR preconditioner
  /*
  //  Precondizionatore SSOR per sistema del tipo
  //
  //  /      \ / \   / \
  //  | A -B | |x| = |b|
  //  | B  A | |y|   |c|
  //  \      / \ /   \ /
  */
  template <typename T>
  class CSSORpreconditioner : public Preco<CSSORpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using CSSORPRECO = CSSORpreconditioner<T>;
    using PRECO      = Preco<CSSORPRECO>;
    #endif

    using real_type = T; //!< type of the elements of the preconditioner

  private:

    real_type omega, omega1;
    integer   maxIter;

    Vector<integer>   L_R;
    Vector<integer>   L_J;
    Vector<real_type> L_A;

    Vector<integer>   U_C;
    Vector<integer>   U_I;
    Vector<real_type> U_A;

    Vector<real_type> D, tmp;

    Vector<integer>   B_R;
    Vector<integer>   B_J;
    Vector<real_type> B_A;

    Vector<integer> Lnnz, Unnz, Bnnz;

    mutable Vector<real_type> br, bi, x, y;

    //! build incomplete LDU decomposition with specified pattern `P`
    template <typename MAT>
    void
    build_SOR( MAT const & A, real_type const & _omega, integer _maxIter ) {

      UTILS_ASSERT0(
        A.is_ordered(),
        "Sparse_tool: CSSORpreconditioner::build_SOR\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        A.nrows() == A.ncols(),
        "Sparse_tool: CSSORpreconditioner::build_SOR\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        A.nrows() > 0,
        "Sparse_tool: CSSORpreconditioner::build_SOR\n"
        "empty matrix\n"
      );

      this -> omega   = _omega;
      this -> omega1  = 1.0/_omega-1.0;
      this -> maxIter = _maxIter;

      // step 0: compute necessary memory
      PRECO::pr_size = A.nrows();
      Lnnz.resize( PRECO::pr_size );
      Unnz.resize( PRECO::pr_size );
      Bnnz.resize( PRECO::pr_size );
      D.resize( PRECO::pr_size );
      tmp.resize( PRECO::pr_size );

      Lnnz.setZero();
      Unnz.setZero();
      Bnnz.setZero();
      D.setZero();

      for ( A.Begin(); A.End(); A.Next() ) {
        integer i{ A.row() };
        integer j{ A.column() };
        ++Bnnz(i);
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: initialize structure
      L_R.resize( PRECO::pr_size + 1 );
      U_C.resize( PRECO::pr_size + 1 );
      B_R.resize( PRECO::pr_size + 1 );

      L_R(0) = U_C(0) = B_R(0) = 0;
      for ( integer i{0}; i < PRECO::pr_size; ++i ) {
        L_R(i+1) = L_R(i) + Lnnz(i);
        U_C(i+1) = U_C(i) + Unnz(i);
        B_R(i+1) = B_R(i) + Bnnz(i);
      }

      // step 2: allocate memory
      L_A.resize( L_R(PRECO::pr_size) );
      L_J.resize( L_R(PRECO::pr_size) );

      U_A.resize( U_C(PRECO::pr_size) );
      U_I.resize( U_C(PRECO::pr_size) );

      B_A.resize( B_R(PRECO::pr_size) );
      B_J.resize( B_R(PRECO::pr_size) );

      L_A.setZero();
      U_A.setZero();
      B_A.setZero();

      // step 3: fill structure
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i{ A.row() };
        integer j{ A.column() };
        B_J(B_R(i)+(--Bnnz(i))) = j;
        if      ( i > j ) L_J(L_R(i)+(--Lnnz(i))) = j;
        else if ( i < j ) U_I(U_C(j)+(--Unnz(j))) = i;
      }

      // step 4: sort structure
      for ( integer i{0}; i < PRECO::pr_size; ++i ) {
        std::sort( &L_J(L_R(i)), &L_J(L_R(i+1)) );
        std::sort( &U_I(U_C(i)), &U_I(U_C(i+1)) );
        std::sort( &B_J(B_R(i)), &B_J(B_R(i+1)) );
      }

      // step 5: insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer   i    { A.row() };
        integer   j    { A.column() };
        real_type rval { A.value().real() };
        real_type ival { A.value().imag() };
        { // costruisco B
          integer lo  { B_R(i)   };
          integer hi  { B_R(i+1) };
          integer len { hi - lo  };
          while ( len > 0 ) {
            integer half { len / 2   };
            integer mid  { lo + half };
            if ( B_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          B_A(lo) = ival;
        }
        if ( i > j ) { // costruisco L
          integer lo  { L_R(i)   };
          integer hi  { L_R(i+1) };
          integer len { hi - lo  };
          while ( len > 0 ) {
            integer half { len / 2   };
            integer mid  { lo + half };
            if ( L_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          L_A(lo) = rval;
        } else if ( i < j ) { // costruisco U
          integer lo  { U_C(j)   };
          integer hi  { U_C(j+1) };
          integer len { hi - lo  };
          while ( len > 0 ) {
            integer half { len / 2 };
            integer mid  { lo + half };
            if ( U_I(mid) < i ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          U_A(lo) = rval;
        } else {
          D(i) = rval;
        }
      }
      for ( integer i{0}; i < PRECO::pr_size; ++i )
        UTILS_ASSERT(
          D(i) != real_type(0),
          "Sparse_tool: CSSORpreconditioner::D({}) = {}, size = {}\n",
          i, D(i), D.size()
        );

      // step 6: allocate working vectors
      br.resize( PRECO::pr_size );
      bi.resize( PRECO::pr_size );
      x.resize( PRECO::pr_size );
      y.resize( PRECO::pr_size );
    }

  public:

    CSSORpreconditioner(void) : Preco<CSSORPRECO>() {}

    template <typename MAT>
    CSSORpreconditioner( MAT const & M, real_type _omega, integer _maxIter ) : Preco<CSSORPRECO>()
    { build_SOR( M, _omega, _maxIter ); }

    //! build the preconditioner from matrix `M` with pattern `P`
    template <typename MAT>
    void
    build( MAT const & M, real_type _omega, integer _maxIter )
    { build_SOR( M, _omega, _maxIter ); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & xc, VECTOR const & bc ) const {
      integer const * pC;  integer const * pI;  real_type const * pUA;
      integer const * pBR; integer const * pBJ; real_type const * pBA;
      integer const * pR;  integer const * pJ;  real_type const * pLA;

      integer k;

      // copia dato in ingresso
      for ( k=0; k < PRECO::pr_size; ++k ) {
        br(k) = bc(k).real();
        bi(k) = bc(k).imag();
      }

      x.setZero();
      y.setZero();
      for ( integer ii{0}; ii < maxIter; ++ii ) {

        // calcolo ((1/omega-1)*D-U)*x + B*y + b --------------------
        pC  = U_C.data();
        pI  = U_I.data();
        pUA = U_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type xk{ x(k) };
          x(k) = br(k) + omega1 * D(k) * xk;
          for ( integer i_cnt{pC[1] - pC[0]}; i_cnt > 0; --i_cnt )
            x(*pI++) -= *pUA++ * xk;
          ++pC;
        }

        // aggiungo B*y;
        pBR = B_R.data();
        pBJ = B_J.data();
        pBA = B_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type tt(0);
          for ( integer i_cnt{pBR[1] - pBR[0]}; i_cnt > 0; --i_cnt )
            tt += *pBA++ * y(*pBJ++);
          x(k) += tt;
          ++pBR;
        };

        // solve (D/omega+L) x(n+1) = x(n)
        pR  = L_R.data();
        pJ  = L_J.data();
        pLA = L_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type tt(0);
          for ( integer i_cnt{pR[1] - pR[0]}; i_cnt > 0; --i_cnt )
            tt += *pLA++ * x(*pJ++);
          x(k) = omega*(x(k)-tt)/D(k);
          ++pR;
        };

        // calcolo ((1/omega-1)*D-U)*y - B*x + c  ------------------------
        pC  = U_C.data();
        pI  = U_I.data();
        pUA = U_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type yk{ y(k) };
          y(k) = bi(k) + omega1 * D(k) * yk;
          for ( integer i_cnt{pC[1] - pC[0]}; i_cnt > 0; --i_cnt )
            y(*pI++) -= *pUA++ * yk;
          ++pC;
        }

        // sottraggo B*x;
        pBR = B_R.data();
        pBJ = B_J.data();
        pBA = B_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type tt(0);
          for ( integer i_cnt{pBR[1] - pBR[0]}; i_cnt > 0; --i_cnt )
            tt += *pBA++ * x(*pBJ++);
          y(k) -= tt;
          ++pBR;
        };

        // solve (D/omega+L) y(n+1) = y(n)
        pR  = L_R.data();
        pJ  = L_J.data();
        pLA = L_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type tt(0);
          for ( integer i_cnt{pR[1] - pR[0]}; i_cnt > 0; --i_cnt )
            tt += *pLA++ * y(*pJ++);
          y(k) = omega*(y(k)-tt)/D(k);
          ++pR;
        };
      }

      for ( integer ii{0}; ii < maxIter; ++ii ) {

        // calcolo ((1/omega-1)*D-L)*y - B*x + c -----------------------
        pR  = L_R.data() + PRECO::pr_size;
        pJ  = L_J.data() + *pR;
        pLA = L_A.data() + *pR;
        k = PRECO::pr_size;
        do {
          --k; --pR;
          real_type tt(0);
          for ( integer i_cnt{pR[1] - pR[0]}; i_cnt > 0; --i_cnt )
            tt += *--pLA * y(*--pJ);
          y(k) = omega1*D(k)*y(k)+bi(k)-tt;
        } while ( k > 0 );

        // sottraggo B*x;
        pBR = B_R.data();
        pBJ = B_J.data();
        pBA = B_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type tt(0);
          for ( integer i_cnt{pBR[1] - pBR[0]}; i_cnt > 0; --i_cnt )
            tt += *pBA++ * x(*pBJ++);
          y(k) -= tt;
          ++pBR;
        };

        // solve (D/omega+U) y(n+1) = y(n)
        pC  = U_C.data() + PRECO::pr_size;
        pI  = U_I.data() + *pC;
        pUA = U_A.data() + *pC;
        k = PRECO::pr_size;
        do {
          --k; --pC;
          real_type yk{ y(k)*omega/D(k) };
          for ( integer i_cnt{pC[1] - pC[0]}; i_cnt > 0; --i_cnt )
            y(*--pI) -= *--pUA * yk;
          y(k) = yk;
        } while ( k > 0 );

        // calcolo ((1/omega-1)*D-L)*x + B*y + b ----------------------
        pR  = L_R.data() + PRECO::pr_size;
        pJ  = L_J.data() + *pR;
        pLA = L_A.data() + *pR;
        k = PRECO::pr_size;
        do {
          --k; --pR;
          real_type tt(0);
          for ( integer i_cnt{pR[1] - pR[0]}; i_cnt > 0; --i_cnt )
            tt += *--pLA * x(*--pJ);
          x(k) = omega1*D(k)*x(k)+br(k)-tt;
        } while ( k > 0 );

        // aggiungo B*y;
        pBR = B_R.data();
        pBJ = B_J.data();
        pBA = B_A.data();
        for ( k=0; k < PRECO::pr_size; ++k ) {
          real_type tt(0);
          for ( integer i_cnt{pBR[1] - pBR[0]}; i_cnt > 0; --i_cnt )
            tt += *pBA++ * y(*pBJ++);
          x(k) += tt;
          ++pBR;
        };

        // solve (D/omega+U) x(n+1) = x(n)
        pC  = U_C.data() + PRECO::pr_size;
        pI  = U_I.data() + *pC;
        pUA = U_A.data() + *pC;
        k = PRECO::pr_size;
        do {
          --k; --pC;
          real_type xk{ x(k)*omega/D(k) };
          for ( integer i_cnt{pC[1] - pC[0]}; i_cnt > 0; --i_cnt )
            x(*--pI) -= *--pUA * xk;
          x(k) = xk;
        } while ( k > 0 );

      }

      // copia risulatato in uscita
      for ( k=0; k < PRECO::pr_size; ++k ) xc(k) = real_type(x(k),y(k));
    }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & xx ) const {
      tmp = xx;
      ass_preco( xx, tmp );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,CSSORpreconditioner<TP> >
  operator / (Vector<T> const & v, CSSORpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,CSSORpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::CSSORpreconditioner;
}
#endif

#endif
