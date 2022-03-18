#pragma once
#ifndef SPARSETOOL_ITERATIVE_HH
#define SPARSETOOL_ITERATIVE_HH

#include "sparse_tool.hh"
#include <iostream>

#include "preconditioner/id.hxx"
#include "preconditioner/diag.hxx"
#include "preconditioner/ildu.hxx"
#include "preconditioner/rildu.hxx"
#include "preconditioner/ilduk.hxx"
#include "preconditioner/sor.hxx"
#include "preconditioner/ssor.hxx"
#include "preconditioner/cssor.hxx"
#include "preconditioner/jacobi.hxx"
#include "preconditioner/hss_ssor.hxx"
#include "preconditioner/hss_ildu.hxx"
#include "preconditioner/hss_cgssor.hxx"
#include "preconditioner/hss_opoly.hxx"
#include "preconditioner/hss_opoly_ssor.hxx"
#include "preconditioner/hss_chebyshev.hxx"

#include "iterative/cg.hxx"
#include "iterative/cg_poly.hxx"
#include "iterative/gmres.hxx"
#include "iterative/bicgstab.hxx"
#include "iterative/cocg.hxx"
#include "iterative/cocr.hxx"

namespace Sparse_tool {

  //!
  //! Build incomplete SOR preconditioner.
  //!
  template <typename MAT, typename real_type>
  void
  separate_LDU(
    MAT const &         A,
    Vector<real_type> & L_A,
    Vector<integer>   & L_R,
    Vector<integer>   & L_J,
    Vector<real_type> & D,
    Vector<real_type> & U_A,
    Vector<integer>   & U_I,
    Vector<integer>   & U_C
  ) {
    UTILS_ASSERT0(
      A.is_ordered(),
      "Sparse_Tool: separate_LDU\n"
      "pattern must be ordered before use\n"
    );
    UTILS_ASSERT0(
      A.nrows() == A.ncols(),
      "Sparse_Tool: separate_LDU\n"
      "only square matrix allowed\n"
    );
    UTILS_ASSERT0(
      A.nrows() > 0,
      "Sparse_Tool: separate_LDU\n"
      "empty matrix\n"
    );

    // step 0: compute necessary memory
    integer nr = A.nrows();
    Vector<integer> Lnnz( nr ), Unnz( nr );

    D.resize( nr );

    Lnnz.setZero();
    Unnz.setZero();
    D.setZero();

    // cout nnz
    for ( A.Begin(); A.End(); A.Next() ) {
      integer i = A.row();
      integer j = A.column();
      if      ( i > j ) ++Lnnz(i);
      else if ( i < j ) ++Unnz(j);
    }

    // step 1: initialize structure
    L_R.resize( nr + 1 );
    U_C.resize( nr + 1 );

    L_R(0) = U_C(0) = 0;
    for ( integer i = 0; i < nr; ++i ) {
      L_R(i+1) = L_R(i) + Lnnz(i);
      U_C(i+1) = U_C(i) + Unnz(i);
    }

    // step 2: allocate memory
    L_A.resize( L_R( nr ) );
    L_J.resize( L_R( nr ) );

    U_A.resize( U_C( nr ) );
    U_I.resize( U_C( nr ) );

    L_A.setZero();
    U_A.setZero();

    // step 3: fill structure
    for ( A.Begin(); A.End(); A.Next() ) {
      integer   i = A.row();
      integer   j = A.column();
      real_type a = A.value();
      if ( i > j ) {
        integer ipos = L_R(i)+(--Lnnz(i));
        L_J(ipos) = j;
        L_A(ipos) = a;
      } else if ( i < j ) {
        integer ipos = U_C(j)+(--Unnz(j));
        U_I(ipos) = i;
        U_A(ipos) = a;
      } else {
        D(i) = a;
      }
    }

    // step 4: sort structure
    for ( integer i = 0; i < nr; ++i ) {
      QuickSortI( L_J.data()+L_R(i), L_A.data()+L_R(i), L_R(i+1)-L_R(i) );
      QuickSortI( U_I.data()+U_C(i), U_A.data()+U_C(i), U_C(i+1)-U_C(i) );
    }
  }

  //!
  //! \f[ \mathbf{x} = (\eta \mathbf{D}+\mathbf{L})^{-1} \mathbf{x} \f]
  //!
  template <typename real_type>
  void
  solve_DL(
    real_type                 eta,
    Vector<real_type> const & D,
    Vector<real_type> const & L_A,
    Vector<integer>   const & L_R,
    Vector<integer>   const & L_J,
    Vector<real_type>       & x
  ) {
    integer nr = D.size();
    // solve (eta*D+L) x(n+1) = x(n)
    integer   const * pR  = L_R.data();
    integer   const * pJ  = L_J.data();
    real_type const * pLA = L_A.data();
    for ( integer k=0; k < nr; ++k ) {
      real_type tt(0);
      for ( integer i_cnt = pR[1] - pR[0]; i_cnt > 0; --i_cnt )
        tt += *pLA++ * x(*pJ++);
      x(k) -= tt;
      x(k) /= eta*D(k);
      ++pR;
    };
  }

  //!
  //! \f[ \mathbf{x} = (\eta \mathbf{D}+\mathbf{U})^{-1} \mathbf{x} \f]
  //!
  template <typename real_type>
  void
  solve_DU(
    real_type                 eta,
    Vector<real_type> const & D,
    Vector<real_type> const & U_A,
    Vector<integer>   const & U_I,
    Vector<integer>   const & U_C,
    Vector<real_type>       & x
  ) {

    // calcolo b-(U+(1-1/omega)*D)*x
    integer nr  = D.size();
    integer nnz = U_I.size();
    integer   const * pC  = U_C.data()+nr;
    integer   const * pI  = U_I.data()+nnz;
    real_type const * pUA = U_A.data()+nnz;

    integer k = nr;
    while ( k-- > 0 ) {
      x(k) /= eta*D(k);
      --pC;
      real_type xk = x(k);
      for ( integer i_cnt = pC[1] - pC[0]; i_cnt > 0; --i_cnt )
        x(*--pI) -= *--pUA * xk;
    }
  }

  /*
  //  ### #       ######  #     #
  //   #  #       #     # #     # # ##### ###### #####
  //   #  #       #     # #     # #   #   #      #    #
  //   #  #       #     # #     # #   #   #####  #    #
  //   #  #       #     # #     # #   #   #      #####
  //   #  #       #     # #     # #   #   #      #   #
  //  ### ####### ######   #####  #   #   ###### #    #
  */
  //!
  //! Incomplete `LDU` preconditioner.
  //!
  template <typename T>
  class ILDUiterPreconditioner : public Preco<ILDUiterPreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef ILDUiterPreconditioner<T> ILDUITERPRECO;
    typedef Preco<ILDUITERPRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner

  private:

    integer                       maxIter;
    ILDUpreconditioner<real_type> P;
    CRowMatrix<real_type>         Mat;
    mutable Vector<real_type>     tmp;

    //!
    //! Build incomplete `LDU` decomposition with specified pattern `P`.
    //!
    template <typename MAT, typename PAT>
    void
    build_ILDUiter( MAT const & A, PAT const & PT, integer iter ) {
      maxIter = iter;
      P.build(A,PT);
      Mat = A;
      tmp.resize(A.nrows());
    }

  public:

    ILDUiterPreconditioner(void) : Preco<ILDUITERPRECO>() {}

    template <typename MAT>
    ILDUiterPreconditioner( MAT const & m_M, integer m_iter ) : Preco<ILDUITERPRECO>()
    { build_ILDUiter( m_M, m_M, m_iter ); }

    template <typename MAT, typename PRE>
    ILDUiterPreconditioner( MAT const & m_M, PRE const & m_P, integer m_iter ) : Preco<ILDUITERPRECO>()
    { build_ILDUiter(m_M,m_P,m_iter); }

    //!
    //! Build the preconditioner from matrix `M`.
    //!
    template <typename MAT>
    void
    build( MAT const & m_M, integer m_iter )
    { build_ILDUiter(m_M,m_M,m_iter); }

    //!
    //! Build the preconditioner from matrix `M` with pattern `P`.
    //!
    template <typename MAT, typename PRE>
    void
    build( MAT const & m_M, PRE const & m_P, integer m_iter )
    { build_ILDUiter(m_M,m_P,m_iter); }

    //!
    //! Apply preconditioner to vector `b`
    //! and store result to vector `x`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & x, VECTOR const & b ) const {
      VECTOR q(b.size()), r(x.size());
      x = b;
      for ( integer i = 0; i < maxIter; ++i ) {
        r = b-Mat*x;
        q = r / P;
        x += real_type(0.5)*q;
        //r -= real_type(0.1)*Mat*q;
      }
    }

    //!
    //! Apply preconditioner to vector `x`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & x ) const {
      tmp = x;
      ass_preco( x, tmp );
    }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,ILDUiterPreconditioner<TP> >
  operator / (Vector<T> const & v, ILDUiterPreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,ILDUiterPreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::ILDUiterPreconditioner;
}
#endif

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
