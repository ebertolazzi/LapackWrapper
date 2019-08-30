#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_OPOLY_SSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_OPOLY_SSOR_HH

using namespace std;

namespace SparseTool {

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

    //! \cond NODOC
    typedef HSS_OPOLY_SSOR_Preconditioner<T> HSS_OPOLY_SSOR_PRECO;
    typedef Preco<HSS_OPOLY_SSOR_PRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner
    typedef typename T::value_type rvalueType; //!< type of the elements of the preconditioner

  private:

    indexType neq, mdegree;
    SSORpreconditioner<rvalueType> preco;
    //IdPreconditioner<rvalueType>   preco;
    CCoorMatrix<rvalueType>        Amat;
    mutable Vector<rvalueType>     s0, s1, As1, y, tmp1, tmp2;

    //! build incomplete LDU decomposition with specified pattern \c P
    template <typename MAT>
    void
    build_HSS_OPOLY_SSOR( MAT const & A, indexType m, rvalueType const & omega, indexType iterssor ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "HSS_OPOLY_SSOR_Preconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "HSS_OPOLY_SSOR_Preconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "HSS_OPOLY_SSOR_Preconditioner::build_LDU empty matrix"
      )
    
      mdegree = m;
      neq     = A.numRows();
      Amat . resize(neq,neq,A.nnz());
      s0   . resize(neq);
      s1   . resize(neq);
      As1  . resize(neq);
      y    . resize(neq);

      tmp1 . resize(neq);
      tmp2 . resize(neq);

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        Amat.insert(i,j) = A.value().real() + A.value().imag();
      }
      
      Amat.internalOrder();
      preco . build( Amat, omega, iterssor );
      //preco . build( Amat );
    }

    //! apply preconditioner to vector \c v and store result in vector \c y
    void
    mulPoly( Vector<rvalueType> & _y, Vector<rvalueType> const & v ) const {
      s0  = rvalueType(1.5)*v;
      As1 = Amat*v;
      preco.assPreco(_y,As1);
      s1 = rvalueType(4.)*v - rvalueType(10./3.)*_y;
      for ( indexType n = 2; n <= mdegree; ++n ) {
        rvalueType delta = ((3*n+6)*n+2.0)/((n+0.5)*(n+1));
        rvalueType a = -4+(6*n+10.0)/((n+2)*(n+2));
        rvalueType b = 2-delta;
        rvalueType c = -1+delta;
        As1 = Amat*s1;
        preco.assPreco(_y,As1);
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
      indexType          m,
      rvalueType const & omega,
      indexType          iterssor
    ) : Preco<HSS_OPOLY_SSOR_PRECO>()
    { build_HSS_OPOLY_SSOR( M, m, omega, iterssor ); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M, indexType m, rvalueType const & omega, indexType iterssor )
    { build_HSS_OPOLY_SSOR( M, m, omega, iterssor ); }

    //! apply preconditioner to vector \c v and store result to vector \c y
    template <typename VECTOR>
    void
    assPreco( VECTOR & _y, VECTOR const & v ) const {
      for ( indexType k=0; k < neq; ++k ) tmp2(k) = v(k).real();
      preco.assPreco( tmp1, tmp2 ); tmp1 *= 0.5;
      mulPoly( tmp2, tmp1 );
      for ( indexType k=0; k < neq; ++k ) { _y(k) = valueType(tmp2(k),_y(k).imag()); tmp2(k) = v(k).imag(); }
      preco.assPreco( tmp1, tmp2 ); tmp1 *= 0.5;
      mulPoly( tmp2, tmp1 );
      valueType cst = valueType(0.5,-0.5);
      for ( indexType k=0; k < neq; ++k ) { _y(k) = valueType(_y(k).real(),tmp2(k)); _y(k) *= cst; }
    }
  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_OPOLY_SSOR_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_OPOLY_SSOR_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_OPOLY_SSOR_Preconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::HSS_OPOLY_SSOR_Preconditioner;
}

#endif
