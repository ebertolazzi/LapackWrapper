#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_CGSSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_CGSSOR_HH

namespace SparseTool {

  /*
  //  #     #  #####   #####   #####   #####   #####   #####  ####### ######  
  //  #     # #     # #     # #     # #     # #     # #     # #     # #     # 
  //  #     # #       #       #       #       #       #       #     # #     # 
  //  #######  #####   #####  #       #  ####  #####   #####  #     # ######  
  //  #     #       #       # #       #     #       #       # #     # #   #   
  //  #     # #     # #     # #     # #     # #     # #     # #     # #    #  
  //  #     #  #####   #####   #####   #####   #####   #####  ####### #     # 
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_CGSSOR_Preconditioner : public Preco<HSS_CGSSOR_Preconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef HSS_CGSSOR_Preconditioner<T> HSS_CGSSOR_PRECO;
    typedef Preco<HSS_CGSSOR_PRECO>      PRECO;
    #endif

    typedef T real_type; //!< type of the elements of the preconditioner
    typedef typename T::value_type rreal_type; //!< type of the elements of the preconditioner

    // variabili per CG
    CCoorMatrix<rreal_type>    Amat;
    mutable Vector<rreal_type> p, q, r, Ap, tmp1, tmp2;
    rreal_type                 epsi;

  private:

    integer neq, maxIter;
    SSORpreconditioner<rreal_type> preco;
    //ILDUpreconditioner<rreal_type> preco;

    //! build incomplete LDU decomposition with specified pattern `P` 
    template <typename MAT>
    void
    build_HSS_CGSSOR(
      MAT        const & A,
      integer          mIter,
      rreal_type const & rtol,
      rreal_type const & omega,
      integer          m
    ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "HSS_CGSSOR_Preconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "HSS_CGSSOR_Preconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "HSS_CGSSOR_Preconditioner::build_LDU empty matrix"
      )

      maxIter = mIter;
      epsi    = rtol;
      neq     = A.numRows();
      Amat . resize(neq,neq,A.nnz());
      p    . resize(neq);
      q    . resize(neq);
      r    . resize(neq);
      tmp1 . resize(neq);
      tmp2 . resize(neq);
      Ap   . resize(neq);

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        integer i = A.row();
        integer j = A.column();
        Amat.insert(i,j) = A.value().real() + A.value().imag();
      }

      Amat  . internalOrder();
      preco . build( Amat, omega, m );
    }

    //! apply preconditioner to vector `v`  and store result to vector `res`    void
    void
    assPrecoR( Vector<rreal_type> & x, Vector<rreal_type> const & b ) const {
      rreal_type rho, rho_1, normr0;
      x.setZero();
      r = b; // parto con x = 0
      preco.assPreco(p,r);
      normr0 = p.template lpNorm<Eigen::Infinity>();
      rho    = p.dot(r);

      for ( integer iter = 0; iter < maxIter; ++iter ) { 
        Ap = Amat * p;
        rreal_type alpha = rho/Ap.dot(p);
        x += alpha * p;
        r -= alpha * Ap;
        preco.assPreco(q,r);
        rreal_type resid = q.template lpNorm<Eigen::Infinity>()/normr0;
        //cout << "       SUB  iter = " << iter << " residual = " << resid << '\n';
        if ( resid <= epsi ) break;
        rho_1 = rho;
        rho   = q.dot(r);
        p = q + (rho / rho_1) * p;
      }
    }

  public:

    HSS_CGSSOR_Preconditioner(void) : Preco<HSS_CGSSOR_PRECO>() {}

    template <typename MAT>
    HSS_CGSSOR_Preconditioner(
      MAT        const & M,
      integer          mIter,
      rreal_type const & rtol,
      rreal_type const & omega,
      integer          m
    ) : Preco<HSS_CGSSOR_PRECO>()
    { build_HSS_CGSSOR( M, mIter, rtol, omega, m ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M, integer mIter, rreal_type const & rtol, rreal_type const & omega, integer m )
    { build_HSS_CGSSOR( M, mIter, rtol, omega, m ); }

    //! apply preconditioner to vector `v`  and store result to vector \c y
    template <typename VECTOR>
    void
    assPreco( VECTOR & _y, VECTOR const & v ) const {

      for ( integer k=0; k < neq; ++k ) tmp1(k) = v(k).real();
      assPrecoR(tmp2,tmp1);
      for ( integer k=0; k < neq; ++k ) _y(k) = real_type(tmp2(k), _y(k).imag());

      for ( integer k=0; k < neq; ++k ) tmp1(k) = v(k).imag();
      assPrecoR(tmp2,tmp1);
      for ( integer k=0; k < neq; ++k ) _y(k) = real_type(_y(k).real(),tmp2(k));

      real_type cst = real_type(0.5,-0.5);
      for ( integer k=0; k < neq; ++k ) _y(k) *= cst;
    }
  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_CGSSOR_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_CGSSOR_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_CGSSOR_Preconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::HSS_CGSSOR_Preconditioner;
}
#endif

#endif
