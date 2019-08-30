/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH PARDISO  |
 |                                                                          |
 |  date         : 2016, December 6                                         |
 |  version      : 1.0.                                                     |
 |  file         : SparseTool_pardiso.hh                                    |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_PARDISO_HH
#define SPARSETOOL_PARDISO_HH

#include "../sparse_tool.hh"
#include <numeric>
#include <algorithm>
#include <complex>

#ifndef F77NAME
  #define F77NAME(A) A##_
#endif

namespace SparseTool {

  // External Fortran routines defined in the PARDISO library.
  extern "C" {

    typedef int    integer;
    typedef double real;

    void
    F77NAME(pardisoinit)(
      void    *,
      integer *,
      integer *,
      integer *,
      real    *,
      integer *
    );

    void
    F77NAME(pardiso)(
      void    *,
      integer *,
      integer *,
      integer *,
      integer *,
      integer *,
      real    *,
      integer *,
      integer *,
      integer *,
      integer *,
      integer *,
      integer *,
      real    *,
      real    *,
      integer *,
      real    *
    );

    void
    F77NAME(pardiso_chkmatrix)(
      integer const *,
      integer const *,
      real    const *,
      integer const *,
      integer const *,
      integer *
    );

    void
    F77NAME(pardiso_chkvec)(
      integer const *,
      integer const *,
      real    const *,
      integer       *
    );

    void
    F77NAME(pardiso_printstats)(
      integer const *,
      integer const *,
      real    const *,
      integer const *,
      integer const *,
      integer const *,
      real    const *,
      integer       *
    );

  };

  // -----------------------------------------------------------------
  // static integer const real_structsym       =  1;
  // static integer const real_symmetric_pdef  =  2;
  // static integer const real_symmetric_indef = -2;
  // static integer const real_nonsym          = 11;
  // static integer const cmpx_structsym       =  3;
  // static integer const cmpx_hermitian_pdef  =  4;
  // static integer const cmpx_hermitian_indef = -4;
  // static integer const cmpx_symmetric       =  6;
  // static integer const cmpx_nonsym          = 13;

  // Class PardisoInfo.
  // -----------------------------------------------------------------
  // An object of this class stores all the parameters and internal data
  // structures that are common to all PARDISO subroutines.
  template <typename real_or_complex>
  class Pardiso {
  protected:

    integer mtype;
    // mtype =  1 real and structurally symmetric, supernode pivoting
    // mtype =  2 real and symmetric positive definite
    // mtype = -2 real and symmetric indefinite, diagonal or Bunch-Kaufman pivoting
    // mtype = 11 real and nonsymmetric, complete supernode pivoting
    // mtype =  3 complex and structurally symmetric, supernode pivoting
    // mtype =  4 complex and hermitian positive definite
    // mtype = -4 complex and hermitian indefinite, diagonal or Bunch-Kaufman pivoting
    // mtype =  6 complex and symmetric
    // mtype = 13 complex and nonsymmetric, supernode pivoting
    
    //  phase = 11  Analysis
    //  phase = 12  Analysis, Numerical Factorization
    //  phase = 13  Analysis, Numerical Factorization, Solve, Iterative Refinement
    //  phase = 22  Numerical Factorization
    //  phase = -22 Selected Inversion
    //  phase = 23  Numerical Factorization, Solve, Iterative Refinement
    //  phase = 33  Solve, Iterative Refinement
    //  phase = -1  Release all internal memory for all matrices
    //  phase =  0  Release memory for matrix number MNUM

    integer maxfct; // Maximum number of numerical factorizations.
    integer mnum;   // Which factorization to use.
    integer msglvl;

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void *pt[64];

    /* Pardiso control parameters. */
    integer iparm[64];
    real    dparm[64];

    integer                 n;
    Vector<real_or_complex> A;
    Vector<integer>         R;
    Vector<integer>         J;
    
    char const *
    error_to_string( integer error ) const {
      switch ( error ) {
        case  0:   return "No error.";
        case -1:   return "Input inconsistent.";
        case -2:   return "Not enough memory.";
        case -3:   return "Reordering problem.";
        case -4:   return "Zero pivot, numerical fact. or iterative refinement problem.";
        case -5:   return "Unclassified (internal) error.";
        case -6:   return "Preordering failed (matrix types 11, 13 only).";
        case -7:   return "Diagonal matrix problem.";
        case -8:   return "32-bit integer overflow problem.";
        case -10:  return "No license file pardiso.lic found.";
        case -11:  return "License is expired.";
        case -12:  return "Wrong username or hostname.";
        case -100: return "Reached maximum number of Krylov-subspace iteration in iterative solver.";
        case -101: return "No sufficient convergence in Krylov-subspace iteration within 25 iterations.";
        case -102: return "Error in Krylov-subspace iteration.";
        case -103: return "Break-Down in Krylov-subspace iteration.";
      }
      return "Unkonwn error.";
    }

  public:

    Pardiso( integer _mtype, integer solver )
    : mtype(_mtype)
    , maxfct(1)
    , mnum(1)
    , msglvl(0)
    {
      int error = 0;
      F77NAME(pardisoinit)( pt, &mtype, &solver, iparm, dparm, &error );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, Error N. " << error <<
        "\n" << error_to_string(error)
      );
      iparm[2] = 1; // Set the number of processors.
    }

    ~Pardiso() {
      free();
    }
    
    template <typename MT>
    void
    load( Sparse<real,MT> const & M ) {
      load( M, all_ok() );
    }

    template <typename MT, typename Compare>
    void
    load( Sparse<real,MT> const & M, Compare cmp ) {
      SPARSETOOL_ASSERT(
        M.numRows() == M.numCols(),
        "Pardiso interface, matrix must be square numRows = " <<
        M.numRows() << " numCols = " <<  M.numCols()
      );

      // step 0: Count nonzero
      n = M.numRows();
      integer nnz = 0;
      for ( M.Begin(); M.End(); M.Next() )
        if ( cmp(M.row(), M.column() ) ) ++nnz;

      A.resize( nnz + 1 );
      J.resize( nnz );
      R.resize( M.numRows() + 1 );
      R = 0;

      // step 1: Evaluate not zero pattern
      for ( M.Begin(); M.End();  M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) ++R(i);
      }
      for ( indexType k = 0; k < n; ++k ) R(k+1) += R(k);

      // step 2: Fill matrix
      for ( M.Begin(); M.End(); M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) {
          indexType ii = --R(i);
          M.assign(A(ii));
          J(ii) = j;
        }
      }

      SPARSETOOL_TEST(
        R(0) == 0 && R(n) == nnz,
        "Pardiso interface, load failed"
      )
      // step 3: internalOrder matrix
      std::vector<int> index;
      indexType ii, rk, rk1;
      for ( ii = 0, rk = R(0); ii < n; ++ii, rk = rk1 ) {
        rk1 = R(ii+1);
        // skip empty rows or row with 1 element
        if ( rk1 > rk+1 )
          QuickSortI<integer,real>( &J(rk), &A(rk), rk1-rk );
      }

      // step 4: translate to 1-based address
      R += 1;
      J += 1;
    }

    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    void
    check_matrix() const {
      integer error = 0;
      F77NAME(pardiso_chkmatrix)( &mtype, &n,
                                  &A.front(), &R.front(), &J.front(), &error );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, in consistency matrix ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }

    /* -------------------------------------------------------------------- */
    /* ..  pardiso_chkvec(...)                                              */
    /*     Checks the given vectors for infinite and NaN values             */
    /*     Input parameters (see PARDISO user manual for a description):    */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    void
    check_rhs( integer nrhs, real_or_complex const b[] ) const {
      integer error = 0;
      F77NAME(pardiso_chkvec)( &n, &nrhs, b, &error );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, in right hand side ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }

    /* -------------------------------------------------------------------- */
    /* .. pardiso_printstats(...)                                           */
    /*    prints information on the matrix to STDOUT.                       */
    /*    Use this functionality only for debugging purposes                */
    /* -------------------------------------------------------------------- */
    void
    stats( integer nrhs, real_or_complex const b[] ) const {
      integer error = 0;
      F77NAME(pardiso_printstats)(
        &mtype, &n,
        &A.front(), &R.front(), &J.front(),
        &nrhs, b, &error
      );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, in printstats ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }
 
    /* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */
    void
    reorder( integer * perm = nullptr ) {
      integer phase = 11;
      integer error = 0;
      integer idum  = 0; // Integer placeholder.
      double  ddum[2];   // Double placeholder.

      // Check whether the user has provided a permutation vector.
      iparm[4] = perm == nullptr ? 0 : 1;

      F77NAME(pardiso)(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, &A.front(), &R.front(), &J.front(), perm, &idum,
        iparm, &msglvl, ddum, ddum, &error, dparm
      );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, during symbolic factorization ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
      // printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
      // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
    }

    /* -------------------------------------------------------------------- */
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */
    void
    factorize() {
      integer phase = 22;
      integer error = 0;
      integer idum  = 0; // Integer placeholder.
      double  ddum[2];   // Double placeholder.
      // iparm[1]  = 0; /* 0 = minimum degree, 2 = nested dissection */
      // iparm[11] = 0; /* 1 = solve transposed system */
      // iparm[17] = -1; /* return the number of nnz in factors */
      // iparm[18] = -1; /* Gflops for the factorization */
      // iparm[20] = -1; /* Pivoting for symmetric indefinite matrices, if 1 1 × 1 and 2 × 2 Bunch-Kaufman pivoting */
      // iparm[21]; /* Inertia: Number of positive eigenvalues. */
      // iparm[22]; /* Inertia: Number of negative eigenvalues. */
      // iparm[32] = 1; /* compute determinant */
      F77NAME(pardiso)(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, &A.front(), &R.front(), &J.front(), &idum, &idum,
        iparm, &msglvl, ddum, ddum, &error, dparm
      );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, during numerical factorization ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }

    /* -------------------------------------------------------------------- */
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */

    // Ask PARDISO to solve the system(s) of equations AX = B, returning
    // the solution(s) in the matrix X.
    void
    solve( integer nrhs, real_or_complex const b[], real_or_complex x[] ) {
      integer phase = 33;
      integer idum  = 0; // Integer placeholder.
      iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
      integer error = 0;
      F77NAME(pardiso)(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, (real*)&A.front(), &R.front(), &J.front(), &idum, &nrhs,
        iparm, &msglvl, (real*)b, (real*)x, &error, dparm
      );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, during solution ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }
    
    void
    solve( Vector<real_or_complex> const & b, Vector<real_or_complex> & x ) {
      x.resize( b.size() );
      solve( 1, &b.front(), &x.front() );
    }

    /* -------------------------------------------------------------------- */    
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */
    void
    free() {
      integer idum  = 0; // Integer placeholder.
      real    ddum[2];   // Double placeholder.
      integer phase = -1; /* Release internal memory. */
      integer error = 0;
      F77NAME(pardiso)(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, ddum, &R.front(), &J.front(), &idum, &idum,
        iparm, &msglvl, ddum, ddum, &error, dparm
      );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, during free memory ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }
  };

  class PardisoRealU : public Pardiso<real> {
  public:

    using Pardiso::check_matrix;
    using Pardiso::check_rhs;
    using Pardiso::stats;
    using Pardiso::load;
    using Pardiso::reorder;
    using Pardiso::factorize;
    using Pardiso::solve;

    PardisoRealU() : Pardiso(11,0)
    {}

  };

  // symmetric but not positive definite
  class PardisoRealS : public Pardiso<real> {
  public:

    using Pardiso::check_matrix;
    using Pardiso::check_rhs;
    using Pardiso::stats;
    using Pardiso::load;
    using Pardiso::reorder;
    using Pardiso::factorize;
    using Pardiso::solve;

    PardisoRealS() : Pardiso(-2,0)
    {}

  };

  // symmetric AND positive definite
  class PardisoRealSPD : public Pardiso<real> {
  public:

    using Pardiso::check_matrix;
    using Pardiso::check_rhs;
    using Pardiso::stats;
    using Pardiso::load;
    using Pardiso::reorder;
    using Pardiso::factorize;
    using Pardiso::solve;

    PardisoRealSPD() : Pardiso(2,0)
    {}

  };

  class PardisoComplexU : public Pardiso<std::complex<real> > {
  public:
  
    typedef Pardiso<std::complex<real> > PARDISO;

    using PARDISO::check_matrix;
    using PARDISO::check_rhs;
    using PARDISO::stats;
    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    PardisoComplexU() : Pardiso(13,0)
    {}

    void
    solve( integer nrhs, std::complex<real> const b[], std::complex<real> x[] )
    { PARDISO::solve( nrhs, b, x ); }
    
    void
    solve(
      Vector<std::complex<real> > const & b,
      Vector<std::complex<real> >       & x
    )
    { PARDISO::solve( b, x ); }

  };

  // symmetric but not positive definite
  class PardisoComplexS : public Pardiso<std::complex<real> > {
  public:
  
    typedef Pardiso<std::complex<real> > PARDISO;

    using PARDISO::check_matrix;
    using PARDISO::check_rhs;
    using PARDISO::stats;
    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    PardisoComplexS() : Pardiso(6,0)
    {}

    void
    solve( integer nrhs, std::complex<real> const b[], std::complex<real> x[] )
    { PARDISO::solve( nrhs, b, x ); }
    
    void
    solve(
      Vector<std::complex<real> > const & b,
      Vector<std::complex<real> >       & x
    )
    { PARDISO::solve( b, x ); }

  };

  // hermitian but not positive definite
  class PardisoComplexH : public Pardiso<std::complex<real> > {
  public:
  
    typedef Pardiso<std::complex<real> > PARDISO;

    using PARDISO::check_matrix;
    using PARDISO::check_rhs;
    using PARDISO::stats;
    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    PardisoComplexH() : Pardiso(-4,0)
    {}

    void
    solve( integer nrhs, std::complex<real> const b[], std::complex<real> x[] )
    { PARDISO::solve( nrhs, b, x ); }
    
    void
    solve(
      Vector<std::complex<real> > const & b,
      Vector<std::complex<real> >       & x
    )
    { PARDISO::solve( b, x ); }

  };

  // hermitian AND positive definite
  class PardisoComplexSPD : public Pardiso<std::complex<real> > {
  public:
  
    typedef Pardiso<std::complex<real> > PARDISO;

    using PARDISO::check_matrix;
    using PARDISO::check_rhs;
    using PARDISO::stats;
    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    PardisoComplexSPD() : Pardiso(4,0)
    {}

    void
    solve( integer nrhs, std::complex<real> const b[], std::complex<real> x[] )
    { PARDISO::solve( nrhs, b, x ); }
    
    void
    solve(
      Vector<std::complex<real> > const & b,
      Vector<std::complex<real> >       & x
    )
    { PARDISO::solve( b, x ); }

  };

}

namespace SparseToolLoad {
  using ::SparseTool::PardisoRealU;
  using ::SparseTool::PardisoRealS;
  using ::SparseTool::PardisoRealSPD;
  using ::SparseTool::PardisoComplexU;
  using ::SparseTool::PardisoComplexS;
  using ::SparseTool::PardisoComplexH;
  using ::SparseTool::PardisoComplexSPD;
}

#endif
