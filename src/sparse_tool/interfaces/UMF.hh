/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH UMFPACK  |
 |                                                                          |
 |  date         : 2008, 6 May                                              |
 |  version      : 1.0.1                                                    |
 |  file         : SparseTool_Umf.hh                                        |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_UMF_HH
#define SPARSETOOL_UMF_HH

#include "../sparse_tool.hh"
#include <umfpack.h>
#include <complex>

namespace SparseTool {

  template <typename T>
  class UMF {
  public:
    typedef int integer;
    typedef T   real_type;

  private:

    // for UMF =====================

    real_type Info[UMFPACK_INFO];
    real_type Control[UMFPACK_CONTROL];
    void    * Symbolic;
    void    * Numeric;

    // for UMF ===================== END

    CColMatrix<T> A_stored;

  public:

    UMF()
    : Symbolic(nullptr)
    , Numeric(nullptr)
    {
      /* get the default control parameters */
      umfpack_di_defaults( Control );

      /* change the default print level for this demo */
      /* (otherwise, nothing will print) */
      Control[UMFPACK_PRL] = 6;
    }

    ~UMF() {
      umfpack_di_free_numeric(&Numeric); Numeric = nullptr;
      umfpack_di_free_symbolic(&Symbolic); Symbolic = nullptr;
    }

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat ) {

      if ( Numeric != nullptr ) {
        umfpack_di_free_numeric(&Numeric);
        Numeric = nullptr;
      }
      if ( Symbolic != nullptr ) {
        umfpack_di_free_symbolic(&Symbolic);
        Symbolic = nullptr;
      }

      A_stored = Mat;

      int status = umfpack_di_symbolic(
        A_stored.numRows(), A_stored.numCols(),
        (int const *)&A_stored.getC().front(),
        (int const *)&A_stored.getI().front(),
        &A_stored.getA().front(),
        &Symbolic, Control, Info
      );

      if (status < 0) {
        umfpack_di_report_info( Control, Info ) ;
        umfpack_di_report_status( Control, status ) ;
      }

      status = umfpack_di_numeric(
        (int const *)&A_stored.getC().front(),
        (int const *)&A_stored.getI().front(),
        &A_stored.getA().front(),
        Symbolic, &Numeric, Control, Info
      );

      if ( status < 0 ) {
        umfpack_di_report_info( Control, Info );
        umfpack_di_report_status( Control, status );
      }

      // Free un-wanted storage
      umfpack_di_free_symbolic( &Symbolic );

      return status;
    }

    int
    solve(
      Vector<T> const & b,
      Vector<T>       & x,
      bool transpose = false
    ) {
      return this->solve( &b.front(), &x.front(), transpose );
    }

    int
    solve(
      real_type const b[],
      real_type       x[],
      bool transpose = false
    ) {
      int status = umfpack_di_solve (
        (transpose ? UMFPACK_At : UMFPACK_A),
        (int const *)&A_stored.getC().front(),
        (int const *)&A_stored.getI().front(),
        &A_stored.getA().front(),
        x, b, Numeric, Control, Info
      );
      if ( status < 0 ) {
        umfpack_di_report_info( Control, Info );
        umfpack_di_report_status( Control, status );
      }
      return status;
    }

  };

  /*
  //  ###  #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //  ###  #######  #####
  */
  /*! \class UMFpreconditioner
      \brief Incomplete \c LU preconditioner
   */
  template <typename T>
  class UMFpreconditioner : public Preco<UMFpreconditioner<T> > {
  public:

    // \cond NODOC
    typedef UMFpreconditioner<T> UMFPRECO;
    typedef Preco<UMFPRECO>      PRECO;
    // \endcond
    typedef T valueType; //!< type of the element of the preconditioner
    typedef typename return_trait<T>::valueType realType;

  private:

    mutable UMF<T> ILU;

  public:

    UMFpreconditioner(void) : Preco<UMFPRECO>() {}
    
    template <typename MAT>
    UMFpreconditioner( MAT const & M, realType const dropTolerance ) : Preco<UMFPRECO>() 
    { ILU.load( M, dropTolerance ); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M, realType const dropTolerance )
    { ILU.load( M, dropTolerance ); }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const
    { ILU.solve( v, res, false ); }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,UMFpreconditioner<TP> >
  operator / (Vector<T> const & v, UMFpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,UMFpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::UMF;
  using ::SparseTool::UMFpreconditioner;
}

#endif
