#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_OPOLY_SSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_OPOLY_SSOR_HH

namespace Sparse_tool {

  /*
  // #     #  #####   #####  ######  ####### #    #     #  #####   #####  ####### ######
  // #     # #     # #     # #     # #     # #     #   #  #     # #     # #     # #     #
  // #     # #       #       #     # #     # #      # #   #       #       #     # #     #
  // #######  #####   #####  ######  #     # #       #     #####   #####  #     # ######
  // #     #       #       # #       #     # #       #          #       # #     # #   #
  // #     # #     # #     # #       #     # #       #    #     # #     # #     # #    #
  // #     #  #####   #####  #       ####### ####### #     #####   #####  ####### #     #
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_OPOLY_SSOR_Preconditioner : public Preco<HSS_OPOLY_SSOR_Preconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef HSS_OPOLY_SSOR_Preconditioner<T> HSS_OPOLY_SSOR_PRECO;
    typedef Preco<HSS_OPOLY_SSOR_PRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner
    typedef typename T::value_type rreal_type; //!< type of the elements of the preconditioner

  private:

    integer neq, mdegree;
    SSORpreconditioner<rreal_type> preco;
    //IdPreconditioner<rreal_type>   preco;
    CCoorMatrix<rreal_type>        Amat;
    mutable Vector<rreal_type>     s0, s1, As1, y, tmp1, tmp2;

    //! build incomplete LDU decomposition with specified pattern `P`
    template <typename MAT>
    void
    build_HSS_OPOLY_SSOR( MAT const & A, integer m, rreal_type const & omega, integer iterssor ) {

      UTILS_ASSERT0(
        A.is_ordered(),
        "Sparse_tool: HSS_OPOLY_SSOR_Preconditioner::build_LDU\n"
        "pattern must be ordered before use\n"
      );
      UTILS_ASSERT0(
        A.nrows() == A.ncols(),
        "Sparse_tool: HSS_OPOLY_SSOR_Preconditioner::build_LDU\n"
        "only square matrix allowed\n"
      );
      UTILS_ASSERT0(
        A.nrows() > 0,
        "Sparse_tool: HSS_OPOLY_SSOR_Preconditioner::build_LDU\n"
        "empty matrix\n"
      );

      mdegree = m;
      neq     = A.nrows();
      Amat.resize(neq,neq,A.nnz());
      s0.resize(neq);
      s1.resize(neq);
      y.resize(neq);
      As1.resize(neq);
      tmp1.resize(neq);
      tmp2.resize(neq);

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        Amat.insert(i,j) = A.value().real() + A.value().imag();
      }

      Amat.internal_order();
      preco.build( Amat, omega, iterssor );
      //preco . build( Amat );
    }

    //! apply preconditioner to vector `v`  and store result in vector \c y
    void
    multiply_by_poly( Vector<rreal_type> & _y, Vector<rreal_type> const & v ) const {
      s0  = rreal_type(1.5)*v;
      As1 = Amat*v;
      preco.ass_preco(_y,As1);
      s1 = rreal_type(4.)*v - rreal_type(10./3.)*_y;
      for ( integer n = 2; n <= mdegree; ++n ) {
        rreal_type delta = ((3*n+6)*n+2.0)/((n+0.5)*(n+1));
        rreal_type a = -4+(6*n+10.0)/((n+2)*(n+2));
        rreal_type b = 2-delta;
        rreal_type c = -1+delta;
        As1 = Amat*s1;
        preco.ass_preco(_y,As1);
        _y = a*(_y-v)+b*s1+c*s0;
        s0 = s1;
        s1 = _y;
      }
    }

  public:

    HSS_OPOLY_SSOR_Preconditioner(void) : Preco<HSS_OPOLY_SSOR_PRECO>() {}

    template <typename MAT>
    HSS_OPOLY_SSOR_Preconditioner(
      MAT        const & M,
      integer          m,
      rreal_type const & omega,
      integer          iterssor
    ) : Preco<HSS_OPOLY_SSOR_PRECO>()
    { build_HSS_OPOLY_SSOR( M, m, omega, iterssor ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M, integer m, rreal_type const & omega, integer iterssor )
    { build_HSS_OPOLY_SSOR( M, m, omega, iterssor ); }

    //! apply preconditioner to vector `v`  and store result to vector \c y
    template <typename VECTOR>
    void
    ass_preco( VECTOR & _y, VECTOR const & v ) const {
      for ( integer k=0; k < neq; ++k ) tmp2(k) = v(k).real();
      preco.ass_preco( tmp1, tmp2 ); tmp1 *= 0.5;
      multiply_by_poly( tmp2, tmp1 );
      for ( integer k=0; k < neq; ++k ) { _y(k) = real_type(tmp2(k),_y(k).imag()); tmp2(k) = v(k).imag(); }
      preco.ass_preco( tmp1, tmp2 ); tmp1 *= 0.5;
      multiply_by_poly( tmp2, tmp1 );
      real_type cst = real_type(0.5,-0.5);
      for ( integer k=0; k < neq; ++k ) { _y(k) = real_type(_y(k).real(),tmp2(k)); _y(k) *= cst; }
    }
  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_OPOLY_SSOR_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_OPOLY_SSOR_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_OPOLY_SSOR_Preconditioner<TP> >(v,P); }
  #endif
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::HSS_OPOLY_SSOR_Preconditioner;
}
#endif

#endif
