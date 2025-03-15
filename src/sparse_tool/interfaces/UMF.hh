/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH UMFPACK  |
 |                                                                          |
 |  date         : 2008, 6 May                                              |
 |  version      : 1.0.1                                                    |
 |  file         : UMF.hh                                                   |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Università degli Studi di Trento                         |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once
#ifndef SPARSETOOL_UMF_dot_HH
#define SPARSETOOL_UMF_dot_HH

#include "../sparse_tool.hh"
#include "../../lapack_wrapper_config.hh"

#include <umfpack.h>
#include <complex>

namespace Sparse_tool {

  template <typename T>
  class UMF {
  public:
    using integer   = int;
    using real_type = T;

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
        A_stored.nrows(), A_stored.ncols(),
        (int const *)A_stored.getC().data(),
        (int const *)A_stored.getI().data(),
        A_stored.getA().data(),
        &Symbolic, Control, Info
      );

      if (status < 0) {
        umfpack_di_report_info( Control, Info );
        umfpack_di_report_status( Control, status );
      }

      status = umfpack_di_numeric(
        (int const *)A_stored.getC().data(),
        (int const *)A_stored.getI().data(),
        A_stored.getA().data(),
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
        (int const *)A_stored.getC().data(),
        (int const *)A_stored.getI().data(),
        A_stored.getA().data(),
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
  //!
  //! Incomplete `LU` preconditioner.
  //!
  template <typename T>
  class UMFpreconditioner : public Preco<UMFpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using UMFPRECO = UMFpreconditioner<T>;
    using PRECO    = Preco<UMFPRECO>;
    #endif
    using real_type = T; //!< type of the element of the preconditioner
    using realType  = typename return_trait<T>::real_type;

  private:

    mutable UMF<T> ILU;

  public:

    UMFpreconditioner(void) : Preco<UMFPRECO>() {}

    template <typename MAT>
    UMFpreconditioner( MAT const & M, realType const dropTolerance ) : Preco<UMFPRECO>()
    { ILU.load( M, dropTolerance ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M, realType const dropTolerance )
    { ILU.load( M, dropTolerance ); }

    //!
    //! Apply preconditioner to vector `v`
    //! and store result to vector `res`
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & res, VECTOR const & v ) const
    { ILU.solve( v, res, false ); }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,UMFpreconditioner<TP> >
  operator / (Vector<T> const & v, UMFpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,UMFpreconditioner<TP> >(v,P); }
  #endif

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::UMF;
  using ::Sparse_tool::UMFpreconditioner;
}
#endif

#endif
