#pragma once
#ifndef SPARSETOOL_ITERATIVE_PRECO_HSS_SSOR_HH
#define SPARSETOOL_ITERATIVE_PRECO_HSS_SSOR_HH

namespace SparseTool {

  /*
  //    #     #  #####   #####  ####### ######
  //    #     # #     # #     # #     # #     #
  //    #     # #       #       #     # #     #
  //    #######  #####   #####  #     # ######
  //    #     #       #       # #     # #   #
  //    #     # #     # #     # #     # #    #
  //    #     #  #####   #####  ####### #     #
  */
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class HSS_SSOR_Preconditioner : public Preco<HSS_SSOR_Preconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef HSS_SSOR_Preconditioner<T> HSS_SSOR_PRECO;
    typedef Preco<HSS_SSOR_PRECO>      PRECO;
    #endif

    typedef T valueType; //!< type of the elements of the preconditioner
    typedef typename T::value_type rvalueType; //!< type of the elements of the preconditioner

  private:

    indexType neq;
    SSORpreconditioner<rvalueType> preco;

    //! build incomplete LDU decomposition with specified pattern `P` 
    template <typename MAT>
    void
    build_HSSOR( MAT const & A, rvalueType const & omega, indexType m ) {

      SPARSETOOL_ASSERT(
        A.isOrdered(),
        "HSS_SSOR_Preconditioner::build_LDU pattern must be ordered before use"
      )
      SPARSETOOL_ASSERT(
        A.numRows() == A.numCols(),
        "HSS_SSOR_Preconditioner::build_LDU only square matrix allowed"
      )
      SPARSETOOL_ASSERT(
        A.numRows() > 0,
        "HSS_SSOR_Preconditioner::build_LDU empty matrix"
      )
    
      neq = A.numRows();
      CCoorMatrix<rvalueType> Amat(neq,neq,A.nnz());

      // insert values
      for ( A.Begin(); A.End(); A.Next() ) {
        indexType i = A.row();
        indexType j = A.column();
        Amat.insert(i,j) = A.value().real() + A.value().imag();
      }
      
      Amat.internalOrder();
      preco.build(Amat,omega,m);
    }

  public:

    HSS_SSOR_Preconditioner(void) : Preco<HSS_SSOR_PRECO>() {}
    
    template <typename MAT>
    HSS_SSOR_Preconditioner( MAT const & M, rvalueType const & omega, indexType m ) : Preco<HSS_SSOR_PRECO>()
    { build_HSSOR( M, omega, m ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M, rvalueType const & omega, indexType m )
    { build_HSSOR( M, omega, m ); }

    //! apply preconditioner to vector `v`  and store result to vector `res`    template <typename VECTOR>
    void
    assPreco( VECTOR & _y, VECTOR const & v ) const {
      preco.assPreco(_y,v);
      valueType cst = valueType(0.5,-0.5);
      for ( indexType k=0; k < neq; ++k ) _y(k) *= cst;
    }
  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,HSS_SSOR_Preconditioner<TP> >
  operator / (Vector<T> const & v, HSS_SSOR_Preconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,HSS_SSOR_Preconditioner<TP> >(v,P); }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SparseToolLoad {
  using ::SparseTool::HSS_SSOR_Preconditioner;
}
#endif

#endif
