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

#ifndef SPARSETOOL_MKL_PARDISO_HH
#define SPARSETOOL_MKL_PARDISO_HH

#include "../sparse_tool.hh"
#include <numeric>
#include <algorithm>
#include <complex>
#include <mkl_pardiso.h>

#ifndef F77NAME
  #define F77NAME(A) A##_
#endif

namespace SparseTool {

  typedef double real;

  // External Fortran routines defined in the PARDISO library.

  // -----------------------------------------------------------------
  // static MKL_INT const real_structsym       =  1;
  // static MKL_INT const real_symmetric_pdef  =  2;
  // static MKL_INT const real_symmetric_indef = -2;
  // static MKL_INT const real_nonsym          = 11;
  // static MKL_INT const cmpx_structsym       =  3;
  // static MKL_INT const cmpx_hermitian_pdef  =  4;
  // static MKL_INT const cmpx_hermitian_indef = -4;
  // static MKL_INT const cmpx_symmetric       =  6;
  // static MKL_INT const cmpx_nonsym          = 13;

  // Class PardisoInfo.
  // -----------------------------------------------------------------
  // An object of this class stores all the parameters and internal data
  // structures that are common to all PARDISO subroutines.
  template <typename real_or_complex>
  class mkl_pardiso {
  protected:

    MKL_INT mtype;
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

    MKL_INT maxfct; // Maximum number of numerical factorizations.
    MKL_INT mnum;   // Which factorization to use.
    MKL_INT msglvl;

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    MKL_INT pt[2*64];

    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    //real    dparm[64];

    MKL_INT                 n;
    Vector<real_or_complex> A;
    Vector<MKL_INT>         R;
    Vector<MKL_INT>         J;
    Vector<MKL_INT>         perm;

    char const *
    error_to_string( MKL_INT error ) const {
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

    mkl_pardiso( MKL_INT _mtype, MKL_INT solver )
    : mtype(_mtype)
    , maxfct(1)
    , mnum(1)
    , msglvl(0)
    {
      pardisoinit( pt, &mtype, iparm );
    }

    ~mkl_pardiso() {
      free();
    }
    
    template <typename MT>
    void
    load( Sparse<real_or_complex,MT> const & M ) {
      load( M, all_ok() );
    }

    template <typename MT, typename Compare>
    void
    load( Sparse<real_or_complex,MT> const & M, Compare cmp ) {
      SPARSETOOL_ASSERT(
        M.numRows() == M.numCols(),
        "Pardiso interface, matrix must be square numRows = " <<
        M.numRows() << " numCols = " <<  M.numCols()
      );

      // step 0: Count nonzero
      n = M.numRows();
      MKL_INT nnz = 0;
      for ( M.Begin(); M.End(); M.Next() )
        if ( cmp(M.row(), M.column() ) ) ++nnz;

      A.resize( nnz + 1 );
      J.resize( nnz );
      R.resize( M.numRows() + 1 );
      perm.resize( M.numRows() );
      R = 0;

      // step 1: Evaluate not zero pattern
      for ( M.Begin(); M.End();  M.Next() ) {
        indexType i = M.row();
        indexType j = M.column();
        if ( cmp(i,j) ) ++R(i);
      }
      for ( indexType k = 0; k < n; ++k ) {
        R(k+1) += R(k);
        perm(k) = k+1;
      }

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
          QuickSortI<MKL_INT,real_or_complex>( &J(rk), &A(rk), rk1-rk );
      }

      // step 4: translate to 1-based address
      R += 1;
      J += 1;
    }

    /* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */
    void
    reorder( MKL_INT * perm = nullptr ) {
      MKL_INT phase = 11;
      MKL_INT error = 0;

      // Check whether the user has provided a permutation vector.
      iparm[4] = perm == nullptr ? 0 : 1;

      pardiso(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, (real*)&A.front(), &R.front(), &J.front(),
        perm, nullptr, iparm, &msglvl,
        nullptr, nullptr, &error
      );

      SPARSETOOL_ASSERT(
        error == 0,
        "mkl_pardiso interface, during symbolic factorization ERROR N. " << error <<
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

      MKL_INT phase = 12;
      MKL_INT error = 0;
      MKL_INT nrhs  = 0;
      // iparm[1]  = 0; /* 0 = minimum degree, 2 = nested dissection */
      // iparm[11] = 0; /* 1 = solve transposed system */
      // iparm[17] = -1; /* return the number of nnz in factors */
      // iparm[18] = -1; /* Gflops for the factorization */
      // iparm[20] = -1; /* Pivoting for symmetric indefinite matrices, if 1 1 × 1 and 2 × 2 Bunch-Kaufman pivoting */
      // iparm[21]; /* Inertia: Number of positive eigenvalues. */
      // iparm[22]; /* Inertia: Number of negative eigenvalues. */
      // iparm[32] = 1; /* compute determinant */
      pardiso(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, (real*)&A.front(), &R.front(), &J.front(),
        &perm.front(), &nrhs, iparm, &msglvl,
        nullptr, nullptr, &error
      );

      SPARSETOOL_ASSERT(
        error == 0,
        "mkl_pardiso interface, during numerical factorization ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }

    /* -------------------------------------------------------------------- */
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */

    // Ask PARDISO to solve the system(s) of equations AX = B, returning
    // the solution(s) in the matrix X.
    void
    solve( MKL_INT nrhs, real_or_complex const b[], real_or_complex x[] ) {
      MKL_INT phase = 33;
      iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
      MKL_INT error = 0;
      pardiso(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, (real*)&A.front(), &R.front(), &J.front(),
        &perm.front(), &nrhs, iparm, &msglvl,
        (real*)b, (real*)x, &error
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
      MKL_INT phase = -1; /* Release internal memory. */
      MKL_INT error = 0;
      MKL_INT nrhs  = 0;
      pardiso(
        pt, &maxfct, &mnum, &mtype, &phase,
        &n, nullptr, &R.front(), &J.front(),
        nullptr, &nrhs, iparm, &msglvl,
        nullptr, nullptr, &error
      );
      SPARSETOOL_ASSERT(
        error == 0,
        "Pardiso interface, during free memory ERROR N. " << error <<
        "\n" << error_to_string(error)
      );
    }
  };

  class mkl_PardisoRealU : public mkl_pardiso<real> {
  public:

    using mkl_pardiso<real>::load;
    using mkl_pardiso<real>::reorder;
    using mkl_pardiso<real>::factorize;
    using mkl_pardiso<real>::solve;

    mkl_PardisoRealU() : mkl_pardiso<real>(11,0)
    {}

  };

  // symmetric but not positive definite
  class mkl_PardisoRealS : public mkl_pardiso<real> {
  public:

    //using Pardiso<real>::check_matrix;
    using mkl_pardiso<real>::load;
    using mkl_pardiso<real>::reorder;
    using mkl_pardiso<real>::factorize;
    using mkl_pardiso<real>::solve;

    mkl_PardisoRealS() : mkl_pardiso<real>(-2,0)
    {}

  };

  // symmetric AND positive definite
  class mkl_PardisoRealSPD : public mkl_pardiso<real> {
  public:

    using mkl_pardiso<real>::load;
    using mkl_pardiso<real>::reorder;
    using mkl_pardiso<real>::factorize;
    using mkl_pardiso<real>::solve;

    mkl_PardisoRealSPD() : mkl_pardiso<real>(2,0)
    {}

  };

  class mkl_PardisoComplexU : public mkl_pardiso<std::complex<real> > {
  public:
  
    typedef mkl_pardiso<std::complex<real> > PARDISO;

    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    mkl_PardisoComplexU() : PARDISO(13,0)
    {}

    void
    solve(
      MKL_INT                  nrhs,
      std::complex<real> const b[],
      std::complex<real>       x[]
    )
    { PARDISO::solve( nrhs, b, x ); }
    
    void
    solve(
      Vector<std::complex<real> > const & b,
      Vector<std::complex<real> >       & x
    )
    { PARDISO::solve( b, x ); }

  };

  // symmetric but not positive definite
  class mkl_PardisoComplexS : public mkl_pardiso<std::complex<real> > {
  public:
  
    typedef mkl_pardiso<std::complex<real> > PARDISO;

    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    mkl_PardisoComplexS() : mkl_pardiso(6,0)
    {}

    void
    solve(
      MKL_INT                  nrhs,
      std::complex<real> const b[],
      std::complex<real>       x[]
    )
    { PARDISO::solve( nrhs, b, x ); }
    
    void
    solve(
      Vector<std::complex<real> > const & b,
      Vector<std::complex<real> >       & x
    )
    { PARDISO::solve( b, x ); }

  };

  // hermitian but not positive definite
  class mkl_PardisoComplexH : public mkl_pardiso<std::complex<real> > {
  public:
  
    typedef mkl_pardiso<std::complex<real> > PARDISO;

    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    mkl_PardisoComplexH() : mkl_pardiso(-4,0)
    {}

    void
    solve(
      MKL_INT                  nrhs,
      std::complex<real> const b[],
      std::complex<real>       x[]
    )
    { PARDISO::solve( nrhs, b, x ); }
    
    void
    solve(
      Vector<std::complex<real> > const & b,
      Vector<std::complex<real> >       & x
    )
    { PARDISO::solve( b, x ); }

  };

  // hermitian AND positive definite
  class mkl_PardisoComplexSPD : public mkl_pardiso<std::complex<real> > {
  public:
  
    typedef mkl_pardiso<std::complex<real> > PARDISO;

    using PARDISO::load;
    using PARDISO::reorder;
    using PARDISO::factorize;

    mkl_PardisoComplexSPD() : mkl_pardiso(4,0)
    {}

    void
    solve(
      MKL_INT                  nrhs,
      std::complex<real> const b[],
      std::complex<real>       x[]
    )
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
  using ::SparseTool::mkl_PardisoRealU;
  using ::SparseTool::mkl_PardisoRealS;
  using ::SparseTool::mkl_PardisoRealSPD;
  using ::SparseTool::mkl_PardisoComplexU;
  using ::SparseTool::mkl_PardisoComplexS;
  using ::SparseTool::mkl_PardisoComplexH;
  using ::SparseTool::mkl_PardisoComplexSPD;
}

#endif
