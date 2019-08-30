#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_CGSSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_CGSSOR_HH

using namespace std;

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

    //! \cond NODOC
    typedef HSS_CGSSOR_Preconditioner<T> HSS_CGSSOR_PRECO;
    typedef Preco<HSS_CGSSOR_PRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner
    typedef typename T::value_type rvalueType; //!< type of the elements of the preconditioner

    // variabili per CG
    CCoorMatrix<rvalueType>    Amat;
    mutable Vector<rvalueType> p, q, r, Ap, tmp1, tmp2;
    rvalueType                 epsi;

  private:

    indexType neq, maxIter;
    SSORpreconditioner<rvalueType> preco;
    //ILDUpreconditioner<rvalueType> preco;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_HSS_CGSSOR(
      MAT        const & A,
      indexType          mIter,
      rvalueType const & rtol,
      rvalueType const & omega,
      indexType          m
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
        indexType i = A.row();
        indexType j = A.column();
        Amat.insert(i,j) = A.value().real() + A.value().imag();
      }

      Amat  . internalOrder();
      preco . build( Amat, omega, m );
    }

    //! apply preconditioner to vector \c v and store result to vector \c res
    void
    assPrecoR( Vector<rvalueType> & x, Vector<rvalueType> const & b ) const {
      rvalueType rho, rho_1, normr0;
      x.setZero();
      r = b; // parto con x = 0
      preco.assPreco(p,r);
      normr0 = normi(p);
      rho    = dot(p,r);

      for ( indexType iter = 0; iter < maxIter; ++iter ) { 
        Ap = Amat * p;
        rvalueType alpha = rho/dot(Ap,p);
        x += alpha * p;
        r -= alpha * Ap;
        preco.assPreco(q,r);
        rvalueType resid = normi(q)/normr0;
        cout << "       SUB  iter = " << iter << " residual = " << resid << '\n';
        if ( resid <= epsi ) break;
        rho_1 = rho;
        rho   = dot(q,r);
        p = q + (rho / rho_1) * p;
      }
    }

  public:

    HSS_CGSSOR_Preconditioner(void) : Preco<HSS_CGSSOR_PRECO>() {}

    template <typename MAT>
    HSS_CGSSOR_Preconditioner(
      MAT        const & M,
      indexType          mIter,
      rvalueType const & rtol,
      rvalueType const & omega,
      indexType          m
    ) : Preco<HSS_CGSSOR_PRECO>()
    { build_HSS_CGSSOR( M, mIter, rtol, omega, m ); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M, indexType mIter, rvalueType const & rtol, rvalueType const & omega, indexType m )
    { build_HSS_CGSSOR( M, mIter, rtol, omega, m ); }

    //! apply preconditioner to vector \c v and store result to vector \c y
    template <typename VECTOR>
    void
    assPreco( VECTOR & _y, VECTOR const & v ) const {

      for ( indexType k=0; k < neq; ++k ) tmp1(k) = v(k).real();
      assPrecoR(tmp2,tmp1);
      for ( indexType k=0; k < neq; ++k ) _y(k) = valueType(tmp2(k), _y(k).imag());

      for ( indexType k=0; k < neq; ++k ) tmp1(k) = v(k).imag();
      assPrecoR(tmp2,tmp1);
      for ( indexType k=0; k < neq; ++k ) _y(k) = valueType(_y(k).real(),tmp2(k));

      valueType cst = valueType(0.5,-0.5);
      for ( indexType k=0; k < neq; ++k ) _y(k) *= cst;
    }
  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_CGSSOR_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_CGSSOR_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_CGSSOR_Preconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::HSS_CGSSOR_Preconditioner;
}

#endif
