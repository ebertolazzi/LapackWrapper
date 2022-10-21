/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH SuperLU  |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : SuperLU.hh                                               |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once
#ifndef SPARSETOOL_SUPERLU_dot_HH
#define SPARSETOOL_SUPERLU_dot_HH

#include "../sparse_tool.hh"
#include "../../lapack_wrapper_config.hh"

#include <complex>

/*
 * SUPERLU interface http://crd.lbl.gov/~xiaoye/SuperLU/
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// workaround for SUPERLU INCLUSION
#define dgemm_  dgemm_BUGGED
#define dtrsv_  dtrsv_BUGGED
#define dtrsm_  dtrsm_BUGGED
#define dgemv_  dgemv_BUGGED
#define sgemm_  sgemm_BUGGED
#define strsv_  strsv_BUGGED
#define strsm_  strsm_BUGGED
#define sgemv_  sgemv_BUGGED
#define xerbla_ xerbla_BUGGED

#ifdef LAPACK_WRAPPER_USE_SYSTEM_SUPERLU
  #include <superlu/slu_sdefs.h>
  #include <superlu/slu_ddefs.h>
#else
  #include "superlu/slu_sdefs.h"
  #include "superlu/slu_ddefs.h"
#endif

#undef dgemm_
#undef dtrsv_
#undef dtrsm_
#undef dgemv_
#undef sgemm_
#undef strsv_
#undef strsm_
#undef sgemv_
#undef xerbla_

#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace Sparse_tool {

  /*\
   |   ____                        _    _   _
   |  / ___| _   _ _ __   ___ _ __| |  | | | |
   |  \___ \| | | | '_ \ / _ \ '__| |  | | | |
   |   ___) | |_| | |_) |  __/ |  | |__| |_| |
   |  |____/ \__,_| .__/ \___|_|  |_____\___/
   |              |_|
  \*/

  template <typename T> struct RealType { using Type = T; };

  template <> struct RealType<std::complex<float> >  { using Type = float; };
  template <> struct RealType<std::complex<double> > { using Type = double; };

  template <typename T>
  class SuperLU {
  public:
    using integer   = int;
    using real_type = typename RealType<T>::Type;

  private:

    integer nRow, nCol, nnz;

    // for SuperLU =====================
    Vector<integer> slu_perm_r; // row permutations from partial pivoting
    Vector<integer> slu_perm_c; // column permutation vector
    Vector<integer> slu_etree;

    superlu_options_t     slu_options;
    mutable SuperLUStat_t slu_stats;
    mutable SuperMatrix   slu_A, slu_AC, slu_L, slu_U; // messo mutable per zittire warning

    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    mutable GlobalLU_t    slu_glu;
    #endif

    // for SuperLU ===================== END

    int
    SuperLU_solve(
      trans_t         trans,
      SuperMatrix   & L,
      SuperMatrix   & U,
      integer         perm_c[],
      integer         perm_r[],
      SuperMatrix   & B,
      SuperLUStat_t & stat
    );

    int
    SuperLU_factor(
      superlu_options_t & options,
      SuperMatrix       & A,
      integer             relax,
      integer             panel_size,
      integer             etree[],
      T                   work[],
      integer             lwork,
      integer             perm_c[],
      integer             perm_r[],
      SuperMatrix       & L,
      SuperMatrix       & U,
      #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      GlobalLU_t        & Glu,
      #endif
      SuperLUStat_t     & stat
    );

    void
    Create_CompCol_Matrix(
      SuperMatrix & A,
      integer       m,
      integer       n,
      integer       nnz,
      T             nzval[],
      integer       rowind[],
      integer       colptr[],
      Stype_t       stype,
      Dtype_t       dtype,
      Mtype_t       mtype
    );

    void
    Create_Dense_Matrix(
      SuperMatrix & X,
      integer       m,
      integer       n,
      T             x[],
      integer       ldx,
      Stype_t       stype,
      Dtype_t       dtype,
      Mtype_t       mtype
    );

  public:

    SuperLU() {}
    ~SuperLU() {}

    int
    load( CColMatrix<T> const & A ) {
      this->nRow = A.nrows();
      this->nCol = A.ncols();
      this->nnz  = A.nnz();
      // Create matrix A in the format expected by SuperLU.
      T const   * AA = static_cast<T const *>(A.getA().data());
      int const * II = reinterpret_cast<int const *>(A.getI().data());
      int const * CC = reinterpret_cast<int const *>(A.getC().data());
      Create_CompCol_Matrix(
        slu_A, A.nrows(), A.ncols(), A.nnz(),
        const_cast<double *>(AA),
        const_cast<int *>(II),
        const_cast<int *>(CC),
        SLU_NC, SLU_D, SLU_GE
      );
      return factorize();
    }

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat ) {
      CColMatrix<T> A(Mat);
      return load( A );
    }

    int
    factorize() {

      set_default_options(&slu_options);

      // Initialize the statistics variables.
      StatInit(&slu_stats);

      slu_perm_r.resize( nRow ); /* row permutations from partial pivoting */
      slu_perm_c.resize( nCol ); /* column permutation vector */
      slu_etree.resize( std::max( nRow, nCol ) );

      /*
       * Get column permutation vector perm_c[], according to permc_spec:
       *   ColPerm = 0: natural ordering
       *   ColPerm = 1: minimum degree on structure of A'*A
       *   ColPerm = 2: minimum degree on structure of A'+A
       *   ColPerm = 3: approximate minimum degree for unsymmetric matrices
       */
      //cout << "get_perm_c.\n";
      get_perm_c(
        slu_options.ColPerm,
        &slu_A,
        slu_perm_c.data()
      );
      //cout << "sp_preorder.\n";
      sp_preorder(
        &slu_options,
        &slu_A,
        slu_perm_c.data(),
        slu_etree.data(),
        &slu_AC
      );

      integer panel_size = sp_ienv(1);
      integer relax      = sp_ienv(2);
      //cout << "dgstrf.\n";
      int info = this->SuperLU_factor(
        slu_options, slu_AC, relax, panel_size,
        slu_etree.data(), nullptr, 0,
        slu_perm_c.data(),
        slu_perm_r.data(),
        slu_L, slu_U,
      #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
        slu_glu,
      #endif
        slu_stats
      );

      // Free un-wanted storage
      Destroy_SuperMatrix_Store(&slu_A);
      Destroy_CompCol_Permuted(&slu_AC);
      StatFree(&slu_stats);
      return info;
    }
    int
    solve(
      Vector<T> const & b,
      Vector<T>       & x,
      bool transpose = false
    ) {
      return this->solve( b.data(), x.data(), transpose );
    }

    int
    solve(
      real_type const b[],
      real_type       x[],
      bool transpose = false
    ) {
      std::copy_n( b, slu_L.nrow, x );
      return solve( x, transpose );
    }

    int
    solve( real_type x[], bool transpose = false ) {

      SuperMatrix BX;

      // Initialize the statistics variables.
      StatInit(&slu_stats);

      int nrow = slu_L.nrow;

      Create_Dense_Matrix(
        BX, nrow, 1, x, nrow, SLU_DN, SLU_D, SLU_GE
      );

      // Solve the system A*X=B, overwriting B with X.
      int info = this->SuperLU_solve(
        (transpose ? TRANS : NOTRANS),
        slu_L, slu_U,
        slu_perm_c.data(),
        slu_perm_r.data(),
        BX, slu_stats
      );

      Destroy_SuperMatrix_Store( &BX );
      StatFree(&slu_stats);

      return info;
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
  //! Incomplete `LU` preconditioner
  //!
  template <typename T>
  class SuperLUpreconditioner : public Preco<SuperLUpreconditioner<T> > {
  public:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    using SLUPRECO = SuperLUpreconditioner<T>;
    using PRECO    = Preco<SLUPRECO>;
    #endif
    using real_type = T; //!< type of the element of the preconditioner
    using realType  = typename RealType<T>::Type;

  private:

    mutable SuperLU<T> ILU;

  public:

    SuperLUpreconditioner(void) : Preco<SLUPRECO>() {}

    template <typename MAT>
    SuperLUpreconditioner( MAT const & M, real_type dropTolerance ) : Preco<SLUPRECO>()
    { ILU.load( M, dropTolerance ); }

    //! build the preconditioner from matrix `M`.
    template <typename MAT>
    void
    build( MAT const & M, real_type dropTolerance )
    { ILU.load( M, dropTolerance ); }

    //!
    //! Apply preconditioner to vector `v`  and store result to vector `res`
    //!
    template <typename VECTOR>
    void
    ass_preco( VECTOR & res, VECTOR const & v ) const
    { ILU.solve( v, res, false ); }

  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,SuperLUpreconditioner<TP> >
  operator / (Vector<T> const & v, SuperLUpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,SuperLUpreconditioner<TP> >(v,P); }
  #endif

  template <>
  inline
  void
  SuperLU<double>::Create_CompCol_Matrix(
    SuperMatrix &            A,
    SuperLU<double>::integer m,
    SuperLU<double>::integer n,
    SuperLU<double>::integer numnonzeros,
    double                   nzval[],
    SuperLU<double>::integer rowind[],
    SuperLU<double>::integer colptr[],
    Stype_t                  stype,
    Dtype_t                  dtype,
    Mtype_t                  mtype
  ) {
    dCreate_CompCol_Matrix(
      &A, m, n, numnonzeros, nzval, rowind, colptr, stype, dtype, mtype
    );
  }

  template <>
  inline
  void
  SuperLU<float>::Create_CompCol_Matrix(
    SuperMatrix &           A,
    SuperLU<float>::integer m,
    SuperLU<float>::integer n,
    SuperLU<float>::integer numnonzeros,
    float                   nzval[],
    SuperLU<float>::integer rowind[],
    SuperLU<float>::integer colptr[],
    Stype_t                 stype,
    Dtype_t                 dtype,
    Mtype_t                 mtype
  ) {
    sCreate_CompCol_Matrix(
      &A, m, n, numnonzeros, nzval, rowind, colptr, stype, dtype, mtype
    );
  }

  template <>
  inline
  void
  SuperLU<double>::Create_Dense_Matrix(
    SuperMatrix &            X,
    SuperLU<double>::integer m,
    SuperLU<double>::integer n,
    double                   x[],
    SuperLU<double>::integer ldx,
    Stype_t                  stype,
    Dtype_t                  dtype,
    Mtype_t                  mtype
  ) {
    dCreate_Dense_Matrix( &X, m, n, x, ldx, stype, dtype, mtype );
  }

  template <>
  inline
  void
  SuperLU<float>::Create_Dense_Matrix(
    SuperMatrix &           X,
    SuperLU<float>::integer m,
    SuperLU<float>::integer n,
    float                   x[],
    SuperLU<float>::integer ldx,
    Stype_t                 stype,
    Dtype_t                 dtype,
    Mtype_t                 mtype
  ) {
    sCreate_Dense_Matrix( &X, m, n, x, ldx, stype, dtype, mtype );
  }

  template <>
  inline
  int
  SuperLU<float>::SuperLU_solve(
    trans_t                 trans,
    SuperMatrix &           L,
    SuperMatrix &           U,
    SuperLU<float>::integer perm_c[],
    SuperLU<float>::integer perm_r[],
    SuperMatrix &           B,
    SuperLUStat_t &         stat
  ) {
    int info;
    sgstrs( trans, &L, &U, perm_c, perm_r, &B, &stat, &info );
    return info;
  }

  template <>
  inline
  int
  SuperLU<double>::SuperLU_solve(
    trans_t                  trans,
    SuperMatrix &            L,
    SuperMatrix &            U,
    SuperLU<double>::integer perm_c[],
    SuperLU<double>::integer perm_r[],
    SuperMatrix &            B,
    SuperLUStat_t &          stat
  ) {
    int info;
    dgstrs( trans, &L, &U, perm_c, perm_r, &B, &stat, &info );
    return info;
  }

  template <>
  inline
  int
  SuperLU<float>::SuperLU_factor(
    superlu_options_t &     options,
    SuperMatrix &           A,
    SuperLU<float>::integer relax,
    SuperLU<float>::integer panel_size,
    SuperLU<float>::integer etree[],
    float                   work[],
    SuperLU<float>::integer lwork,
    SuperLU<float>::integer perm_c[],
    SuperLU<float>::integer perm_r[],
    SuperMatrix &           L,
    SuperMatrix &           U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    GlobalLU_t &            Glu,
    #endif
    SuperLUStat_t &         stat
  ) {
    int info;
    sgstrf(
      &options, &A, relax, panel_size, etree, work, lwork,
      perm_c, perm_r, &L, &U,
      #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      &Glu,
      #endif
      &stat, &info
    );
    return info;
  }

  template <>
  inline
  int
  SuperLU<double>::SuperLU_factor(
    superlu_options_t &      options,
    SuperMatrix &            A,
    SuperLU<double>::integer relax,
    SuperLU<double>::integer panel_size,
    SuperLU<double>::integer etree[],
    double                   work[],
    SuperLU<double>::integer lwork,
    SuperLU<double>::integer perm_c[],
    SuperLU<double>::integer perm_r[],
    SuperMatrix &            L,
    SuperMatrix &            U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    GlobalLU_t &             Glu,
    #endif
    SuperLUStat_t &          stat
  ) {
    int info;
    dgstrf(
      &options, &A, relax, panel_size, etree, work, lwork,
      perm_c, perm_r, &L, &U,
      #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      &Glu,
      #endif
      &stat, &info
    );
    return info;
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::SuperLU;
  using ::Sparse_tool::SuperLUpreconditioner;
}
#endif

#endif
