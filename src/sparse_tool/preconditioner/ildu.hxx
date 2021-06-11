#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_ILDU_HH
#define SPARSETOOL_ITERATIVE_PRECO_ILDU_HH

namespace SparseTool {

  /*
  //  ###  #       ######  #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //   #   #       #     # #     #
  //  ###  ####### ######   #####
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class ILDUpreconditioner : public Preco<ILDUpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef ILDUpreconditioner<T> ILDUPRECO;
    typedef Preco<ILDUPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner

  private:

    Vector<integer> L_R;
    Vector<integer> L_J;
    Vector<real_type> L_A;

    Vector<integer> U_C;
    Vector<integer> U_I;
    Vector<real_type> U_A;

    Vector<real_type> W;
    Vector<real_type> D;

    Vector<integer> Lnnz, Unnz;

    //! build incomplete LDU decomposition with specified pattern `P` 
    template <typename MAT, typename PAT>
    void
    build_ILDU( MAT const & A, PAT const & P ) {

      SPARSETOOL_ASSERT(
        P.isOrdered(),
        "ILDUpreconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        P.numRows() == A.numRows() && P.numCols() == A.numCols(),
        "ILDUpreconditioner::build_LDU pattern do not match matrix size"
      )
      SPARSETOOL_ASSERT(
        P.numRows() == P.numCols(),
        "ILDUpreconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        P.numRows() > 0,
        "ILDUpreconditioner::build_LDU empty matrix"
      )

      // step 0: compute necessary memory
      PRECO::pr_size = A.numRows();
      Lnnz.resize( PRECO::pr_size );
      Unnz.resize( PRECO::pr_size );

      Lnnz.setZero();
      Unnz.setZero();

      for ( P.Begin(); P.End(); P.Next() ) {
        integer i = P.row();
        integer j = P.column();
        if      ( i > j ) ++Lnnz(i);
        else if ( i < j ) ++Unnz(j);
      }

      // step 1: initialize structure
      L_R.resize( PRECO::pr_size + 1 );
      U_C.resize( PRECO::pr_size + 1 );

      L_R(0) = U_C(0) = 0;
      for ( integer i = 0; i < PRECO::pr_size; ++i ) {
        L_R(i+1) = L_R(i) + Lnnz(i);
        U_C(i+1) = U_C(i) + Unnz(i);
      }

      // step 2: allocate memory
      L_A.resize( L_R(PRECO::pr_size) );
      L_J.resize( L_R(PRECO::pr_size) );

      U_A.resize( U_C(PRECO::pr_size) );
      U_I.resize( U_C(PRECO::pr_size) );

      D.resize( PRECO::pr_size );
      W.resize( PRECO::pr_size );

      D.fill(1);
      W.setZero();
      L_A.setZero();
      U_A.setZero();
      
      // step 3: fill structure
      for ( P.Begin(); P.End(); P.Next() ) {
        integer i = P.row();
        integer j = P.column();
        if      ( i > j ) { integer ii = --Lnnz(i); L_J(L_R(i)+ii) = j; }
        else if ( i < j ) { integer jj = --Unnz(j); U_I(U_C(j)+jj) = i; }
      }

      // step 4: sort structure
      for ( integer i = 0; i < PRECO::pr_size; ++i ) {
        std::sort( &L_J(L_R(i)), &L_J(L_R(i+1)) );
        std::sort( &U_I(U_C(i)), &U_I(U_C(i+1)) );
      }

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        typename MAT::real_type const val = A.value(); // (i,j);
        if ( i > j ) {
          integer lo  = L_R(i);
          integer hi  = L_R(i+1);
          integer len = hi - lo;
          while ( len > 0 ) {
            integer half = len / 2;
            integer mid  = lo + half;
            if ( L_J(mid) < j ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          L_A(lo) = val;
        } else if ( i < j ) {
          integer lo  = U_C(j);
          integer hi  = U_C(j+1);
          integer len = hi - lo;
          while ( len > 0 ) {
            integer half = len / 2;
            integer mid  = lo + half;
            if ( U_I(mid) < i ) { lo = mid + 1; len -= half + 1; }
            else                  len = half;
          }
          U_A(lo) = val;
        } else {
          D(i) = val;
        }
      }

      /*
      //
      //  +---+---+ +---+---+ +---+---+    +-----+-------+
      //  | L | 0 | | D | 0 | | U | u |    | LDU |  LDu  |
      //  +---+---+ +---+---+ +---+---+ =  +-----+-------+
      //  |l^T| 1 | | 0 | d | | 0 | 1 |    |l^TDU|d+l^TDu|
      //  +---+---+ +---+---+ +---+---+    +-----+-------+
      //
      //      +-----+-----+
      //      | M11 | M12 |
      //  M = +-----+-----+
      //      | M21 | M22 |
      //      +-----+-----+
      //
      //  l^T D U = M21 ==> U^T D^T l = M21^T
      //
      //  l^T = D^(-1)U^(-T) M21^T
      //  u   = D^(-1)L^(-1) M12
      //  d   = M22 - lDu
      //
      */

      // build LDU decomposition
      integer LRk1 = L_R(1);
      integer UCk1 = U_C(1);
      for ( integer k = 1; k < PRECO::pr_size; ++k ) {
        integer kk;
        integer LRk = LRk1; LRk1 = L_R(k+1);
        integer UCk = UCk1; UCk1 = U_C(k+1);

        //  W = M21^T  ---- l^T = D^(-1)U^(-T) M21^T
        for ( kk = LRk; kk < LRk1; ++kk ) W(L_J(kk)) = L_A(kk);

        //  W = U^(-T) W ---- l^T = D^(-1)U^(-T) M21^T
        #define LDU_FAST
        #ifdef LDU_FAST
        for ( kk = LRk; kk < LRk1; ++kk ) {
          integer j = L_J(kk);
        #else
        for ( integer j = 1; j < k; ++j ) {
        #endif
          integer UCj  = U_C(j);
          integer UCj1 = U_C(j+1);
          real_type bf   = 0;
          for ( integer jj = UCj; jj < UCj1; ++jj ) bf += W(U_I(jj))*U_A(jj);
          W(j) -= bf;
        }
        //  l^T = D^(-1) W;   W = 0---- l^T = D^(-1)U^(-T) M21^T
        #ifdef LDU_FAST
        for ( kk = LRk; kk < LRk1; ++kk ) { integer j = L_J(kk); L_A(kk) = W(j) / D(j); W(j) = 0; }
        #else
        for ( kk = LRk; kk < LRk1; ++kk ) L_A(kk) = W(L_J(kk)) / D(L_J(kk));
        W = 0;
        #endif
        
        //  W = M12  ---- u   = D^(-1)L^(-1) M12
        for ( kk = UCk; kk < UCk1; ++kk ) W(U_I(kk)) = U_A(kk);

        //  W = L^(-1) W  ---- u   = D^(-1)L^(-1) M12
        #ifdef LDU_FAST
        for ( kk = UCk; kk < UCk1; ++kk ) {
          integer i = U_I(kk);
        #else
        for ( integer i = 1; i < k; ++i ) {
        #endif
          integer LRi  = L_R(i);
          integer LRi1 = L_R(i+1);
          real_type bf = 0;
          for ( integer ii = LRi; ii < LRi1; ++ii ) bf += W(L_J(ii))*L_A(ii);
          W(i) -= bf;
        }

        real_type bf = 0;
        for ( kk = LRk; kk < LRk1; ++kk ) bf += L_A(kk) * W(L_J(kk));
        D(k) -= bf;

        #ifdef LDU_FAST
        for ( kk = UCk; kk < UCk1; ++kk ) { integer i = U_I(kk); U_A(kk) = W(i) / D(i); W(i) = 0; }
        #else
        for ( kk = UCk; kk < UCk1; ++kk ) U_A(kk) = W(U_I(kk)) / D(U_I(kk));
        W = 0;
        #endif

        SPARSETOOL_ASSERT( D(k) != real_type(0), "ILDUpreconditioner found D(" << k << ") == 0!" );
      }
    }

  public:

    ILDUpreconditioner(void) : Preco<ILDUPRECO>() {}
    
    template <typename MAT>
    ILDUpreconditioner( MAT const & M ) : Preco<ILDUPRECO>() 
    { build_ILDU( M, M ); }

    template <typename MAT, typename PRE>
    ILDUpreconditioner( MAT const & M, PRE const & P ) : Preco<ILDUPRECO>() 
    { build_ILDU(M,P); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M )
    { build_ILDU(M,M); }

    //! build the preconditioner from matrix `M` with pattern `P` 
    template <typename MAT, typename PRE>
    void
    build( MAT const & M, PRE const & P )
    { build_ILDU(M,P); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const {
      res = v;

      // solve L
      integer const *   pR  = L_R.data();
      integer const *   pJ  = L_J.data();
      real_type const * pLA = L_A.data();
      integer k;

      for ( k=1; k < PRECO::pr_size; ++k ) {
        ++pR;
        typename VECTOR::real_type tmp(0);
        for ( integer i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
          tmp += *pLA++ * res(*pJ++);
        res(k) -= tmp;
      };

      // solve D
      for ( k = 0; k < PRECO::pr_size; ++k ) res(k) /= D(k);

      // solve U
      integer const *   pC  = U_C.data() + PRECO::pr_size;
      integer const *   pI  = U_I.data() + *pC;
      real_type const * pUA = U_A.data() + *pC;

      do {
        typename VECTOR::real_type resk = res(--k);
        --pC;
        for ( integer i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
          res(*--pI) -= *--pUA * resk;
      } while ( k > 1 );

    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,ILDUpreconditioner<TP> >
  operator / (Vector<T> const & v, ILDUpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,ILDUpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::ILDUpreconditioner;
}
#endif

#endif
