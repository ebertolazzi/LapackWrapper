/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2019                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: eig.hxx
///

namespace lapack_wrapper {

  /*\
  :|:   _____ _                            _
  :|:  | ____(_) __ _  ___ _ ____   ____ _| |_   _  ___  ___
  :|:  |  _| | |/ _` |/ _ \ '_ \ \ / / _` | | | | |/ _ \/ __|
  :|:  | |___| | (_| |  __/ | | \ V / (_| | | |_| |  __/\__ \
  :|:  |_____|_|\__, |\___|_| |_|\_/ \__,_|_|\__,_|\___||___/
  :|:           |___/
  \*/

  //! Sparse Matrix Structure
  template <typename T>
  class Eigenvalues {
  public:
    using real_type    = T;
    using complex_type = std::complex<T>;
    using MatW         = MatrixWrapper<T>;
    using Sparse       = SparseCCOOR<T>;

  private:
    Malloc<real_type> m_mem{"Eigenvalues::m_mem"};

    integer     m_N{0};
    integer     m_Lwork{0};
    real_type * m_Re{nullptr};
    real_type * m_Im{nullptr};
    real_type * m_Work{nullptr};
    real_type * m_A_saved{nullptr};

    void allocate( integer N );
    void compute( );

  public:

    Eigenvalues() = default;

    Eigenvalues( integer NRC, real_type const data[], integer ldData ) {
      this->setup( NRC, data, ldData );
    }

    explicit
    Eigenvalues( MatW const & M ) {
      this->setup( M );
    }

    Eigenvalues(
      integer         NRC,
      integer         nnz,
      real_type const values[],
      integer   const row[],
      integer   const col[]
    ) {
      this->setup( NRC, nnz, values, row, col );
    }

    void
    setup( integer NRC, real_type const data[], integer ldData );

    void
    setup( MatW const & M );

    void
    setup(
      integer         NRC,
      integer         nnz,
      real_type const values[],
      integer   const row[],
      integer   const col[]
    );

    void getEigenvalue( integer n, real_type & re, real_type & im ) const;
    void getEigenvalue( integer n, complex_type & eig ) const;

    void getEigenvalues( std::vector<real_type> & re, std::vector<real_type> & im ) const;
    void getEigenvalues( std::vector<complex_type> & eigs ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class Eigenvectors {
  public:
    using real_type    = T;
    using complex_type = std::complex<T>;
    using MatW         = MatrixWrapper<T>;
    using Sparse       = SparseCCOOR<T>;

  private:
    Malloc<real_type> m_mem{"Eigenvectors::mem_real"};

    integer     m_N{0};
    integer     m_Lwork{0};
    real_type * m_Re{nullptr};
    real_type * m_Im{nullptr};
    real_type * m_A_saved{nullptr};
    real_type * m_VL{nullptr};
    real_type * m_VR{nullptr};
    real_type * m_Work{nullptr};

    void allocate( integer N );
    void compute( );

  public:

    Eigenvectors() = default;

    Eigenvectors( integer NRC, real_type const A[], integer ldA ) {
      this->setup( NRC, A, ldA );
    }

    explicit
    Eigenvectors( MatW const & A ) {
      this->setup( A );
    }

    Eigenvectors(
      integer         NRC,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[]
    ) {
      this->setup( NRC, A_nnz, A_values, A_row, A_col );
    }

    void
    setup( integer NRC, real_type const A[], integer ldA );

    void setup( MatW const & A );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[]
    );

    void getEigenvalue( integer n, real_type & re, real_type & im ) const;
    void getEigenvalue( integer n, complex_type & eig ) const;

    void getEigenvalues( std::vector<real_type> & re, std::vector<real_type> & im ) const;
    void getEigenvalues( std::vector<complex_type> & eigs ) const;

    void getLeftEigenvector( std::vector<std::vector<complex_type> > & vecs ) const;
    void getRightEigenvector( std::vector<std::vector<complex_type> > & vecs ) const;
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class GeneralizedEigenvalues {
  public:
    using real_type    = T;
    using complex_type = std::complex<T>;
    using MatW         = MatrixWrapper<T>;
    using Sparse       = SparseCCOOR<T>;

  private:
    Malloc<real_type> m_mem{"GeneralizedEigenvalues::mem_real"};

    integer     m_N{0};
    integer     m_Lwork{0};
    real_type * m_alphaRe{nullptr};
    real_type * m_alphaIm{nullptr};
    real_type * m_beta{nullptr};
    real_type * m_Work{nullptr};
    real_type * m_A_saved{nullptr};
    real_type * m_B_saved{nullptr};

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvalues() = default;

    GeneralizedEigenvalues(
      integer NRC,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    ) {
      this->setup( NRC, A, ldA, B, ldB );
    }

    GeneralizedEigenvalues( MatW const & A, MatW const & B ) {
      this->setup( A, B );
    }

    GeneralizedEigenvalues(
      integer         NRC,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      real_type const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    ) {
      this->setup(
        NRC,
        A_nnz, A_values, A_row, A_col,
        B_nnz, B_values, B_row, B_col
      );
    }

    void
    setup(
      integer NRC,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      real_type const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void getEigenvalue( integer n, real_type & re, real_type & im ) const;
    void getEigenvalue( integer n, complex_type & eig ) const;

    void getEigenvalues( std::vector<real_type> & re, std::vector<real_type> & im ) const;
    void getEigenvalues( std::vector<complex_type> & eigs ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class GeneralizedEigenvectors {
  public:
    using real_type    = T;
    using complex_type = std::complex<T>;
    using MatW         = MatrixWrapper<T>;
    using Sparse       = SparseCCOOR<T>;

  private:
    Malloc<real_type> m_mem_real{"GeneralizedEigenvectors::mem_real"};
    Malloc<integer>   m_mem_int{"GeneralizedEigenvectors::mem_int"};

    integer     m_N{0};
    integer     m_Lwork{0};
    integer     m_ilo{0};
    integer     m_ihi{0};
    real_type   m_abnorm{0};
    real_type   m_bbnorm{0};
    real_type * m_alphaRe{nullptr};
    real_type * m_alphaIm{nullptr};
    real_type * m_beta{nullptr};
    real_type * m_A_saved{nullptr};
    real_type * m_B_saved{nullptr};
    real_type * m_VL{nullptr};
    real_type * m_VR{nullptr};
    real_type * m_lscale{nullptr};
    real_type * m_rscale{nullptr};
    real_type * m_rconde{nullptr};
    real_type * m_rcondv{nullptr};
    real_type * m_Work{nullptr};
    integer   * m_iWork{nullptr};
    integer   * m_bWork{nullptr};

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvectors() = default;

    GeneralizedEigenvectors(
      integer NRC,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    ) {
      this->setup( NRC, A, ldA, B, ldB );
    }

    GeneralizedEigenvectors( MatW const & A, MatW const & B ) {
      this->setup( A, B );
    }

    GeneralizedEigenvectors(
      integer         NRC,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      real_type const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    ) {
      this->setup(
        NRC,
        A_nnz, A_values, A_row, A_col,
        B_nnz, B_values, B_row, B_col
      );
    }

    void
    setup(
      integer NRC,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      real_type const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void getEigenvalue( integer n, real_type & re, real_type & im ) const;
    void getEigenvalue( integer n, complex_type & eig ) const;

    void getEigenvalues( std::vector<real_type> & re, std::vector<real_type> & im ) const;
    void getEigenvalues( std::vector<complex_type> & eigs ) const;

    void getLeftEigenvector( std::vector<std::vector<complex_type> > & vecs ) const;
    void getRightEigenvector( std::vector<std::vector<complex_type> > & vecs ) const;

    real_type         balancedAnorm1()    const { return m_abnorm; }
    real_type         balancedBnorm1()    const { return m_bbnorm; }
    real_type const * getLscale()         const { return m_lscale; }
    real_type const * getRscale()         const { return m_rscale; }
    real_type const * RcondEigenvalues()  const { return m_rconde; }
    real_type const * RcondEigenvectors() const { return m_rcondv; }
  };

}

///
/// eof: eig.hxx
///
