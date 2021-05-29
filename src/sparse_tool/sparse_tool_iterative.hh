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

namespace SparseTool {

  /*
  //  ### #       ######  #     #                       
  //   #  #       #     # #     # # ##### ###### #####  
  //   #  #       #     # #     # #   #   #      #    # 
  //   #  #       #     # #     # #   #   #####  #    # 
  //   #  #       #     # #     # #   #   #      #####  
  //   #  #       #     # #     # #   #   #      #   #  
  //  ### ####### ######   #####  #   #   ###### #    # 
  */
  //! Incomplete `LDU` preconditioner
  template <typename T>
  class ILDUiterPreconditioner : public Preco<ILDUiterPreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    typedef ILDUiterPreconditioner<T> ILDUITERPRECO;
    typedef Preco<ILDUITERPRECO>      PRECO;
    #endif

    typedef T valueType; //!< type of the elements of the preconditioner

  private:

    indexType                      maxIter;
    RILDUpreconditioner<valueType> P;
    CRowMatrix<valueType>          Mat;

    //!
    //! Build incomplete `LDU` decomposition with specified pattern `P`.
    //!
    template <typename MAT, typename PAT>
    void
    build_ILDUiter( MAT const & A, PAT const & PT, indexType iter ) {
      maxIter = iter;
      P.build(A,PT);
      Mat = A;
      //q.resize(A.numRows());
      //r.resize(A.numRows());
    }

  public:

    ILDUiterPreconditioner(void) : Preco<ILDUITERPRECO>() {}
    
    template <typename MAT>
    ILDUiterPreconditioner( MAT const & _M, indexType _iter ) : Preco<ILDUITERPRECO>() 
    { build_ILDUiter( _M, _M, _iter ); }

    template <typename MAT, typename PRE>
    ILDUiterPreconditioner( MAT const & _M, PRE const & _P, indexType _iter ) : Preco<ILDUITERPRECO>() 
    { build_ILDUiter(_M,_P,_iter); }

    //!
    //! Build the preconditioner from matrix `M`.
    //!
    template <typename MAT>
    void
    build( MAT const & _M, indexType _iter )
    { build_ILDUiter(_M,_M,_iter); }

    //!
    //! Build the preconditioner from matrix `M` with pattern `P`.
    //!
    template <typename MAT, typename PRE>
    void
    build( MAT const & _M, PRE const & _P, indexType _iter )
    { build_ILDUiter(_M,_P,_iter); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`.
    //!
    template <typename VECTOR>
    void
    assPreco( VECTOR & x, VECTOR const & b ) const {
      VECTOR q(b.size()), r(x.size());
      x = b;
      for ( indexType i = 0; i < maxIter; ++i ) {
        r = b-Mat*x;
        q = r / P;
        x += valueType(0.5)*q;
        //r -= valueType(0.1)*Mat*q;
      }
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
namespace SparseToolLoad {
  using ::SparseTool::ILDUiterPreconditioner;
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
