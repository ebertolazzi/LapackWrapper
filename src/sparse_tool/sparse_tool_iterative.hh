/*!

  \file     SparseTool_Iterative.hh
  \mainpage SparseTool: a Sparse Matrix Manager 
  \date     2011, July 21
  \version  1.0
  \note     based on SparseLib 11.0 (C 1999 -- 2008)

  \author Enrico Bertolazzi
 
  \par Affiliations:
       Dipartimento di Ingegneria Industriale<BR>
       Universita` degli Studi di Trento<BR>
       email: enrico.bertolazzi@unitn.it<BR>
  
  \section The Preconditioner Classes and Iterative algorithms
  
  The library consists of the following templated classes:
  - \c IdPreconditioner\<T\> which implements the identity preconditioner.
  - \c Dpreconditioner\<T\> which implements the diagonal preconditioner.
  - \c ILDUpreconditioner\<T\> which implement an incomplete \a LDU preconditioner.

  A set of template iterative solvers are available:
  
  - \c cg implementing the cojugate gradient solver
  - \c bicgstab implementing the Bi-conjugate stabilized 
        solver of Van Der Vorst.
  - \c gmres implementing generalized minimal residual 
       of Saad-Shultz
*/

/*!

  \page P2 Iterative solvers
  Here is an example of the use of the iterative solver:

\code
  CRowMatrix<double> A;
  Vector<double>     b, x;
  ILDUpreconditioner P;
  double             tolerance;
  int                maxIter, iter;
  .
  .
  .
  .
  double residual = cg(A, b, x, P, tolerance, maxIter, iter);
  
  double residual = bicgstab(A, b, x, P, tolerance, maxIter, iter);

  double residual = gmres(A, b, x, P, tolerance, maxSubIter, maxIter, iter);
\endcode

  In the example

  - \c A : is the coefficients matrix;
  - \c b : is the known vector;
  - \c x : is the vector which will contains the solution;
  - \c P : is the preconditioner object class;
  - \c tolerance : is the admitted tolerance;
  - \c maxSubIter : for \c gmres is the maximum number of
    iteration before restarting;
  - \c maxIter : is the maximum number of allowable iterations;
  - \c iter : is the number of iterations done;
  - \c residual : the residual of the approximated solution;
*/

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
  //! Incomplete \c LDU preconditioner
  template <typename T>
  class ILDUiterPreconditioner : public Preco<ILDUiterPreconditioner<T> > {
  public:

    //! \cond NODOC
    typedef ILDUiterPreconditioner<T> ILDUITERPRECO;
    typedef Preco<ILDUITERPRECO>      PRECO;

    //! \endcond
    typedef T valueType; //!< type of the elements of the preconditioner

  private:

    indexType                      maxIter;
    RILDUpreconditioner<valueType> P;
    CRowMatrix<valueType>          Mat;

    //! build incomplete LDU decomposition with specified pattern \c P
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

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & _M, indexType _iter )
    { build_ILDUiter(_M,_M,_iter); }

    //! build the preconditioner from matrix \c M with pattern \c P
    template <typename MAT, typename PRE>
    void
    build( MAT const & _M, PRE const & _P, indexType _iter )
    { build_ILDUiter(_M,_P,_iter); }

    //! apply preconditioner to vector \c v and store result to vector \c res
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

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,ILDUiterPreconditioner<TP> >
  operator / (Vector<T> const & v, ILDUiterPreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,ILDUiterPreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::ILDUiterPreconditioner;
}

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
