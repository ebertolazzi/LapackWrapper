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
 |      Universita` degli Studi di Trento                                   |
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
    Malloc<real_type> m_mem;

    integer     m_N;
    integer     m_Lwork;
    real_type * m_Re;
    real_type * m_Im;
    real_type * m_Work;
    real_type * m_A_saved;

    void allocate( integer N );
    void compute( );

  public:

    Eigenvalues();
    Eigenvalues( integer NRC, real_type const data[], integer ldData );
    explicit
    Eigenvalues( MatW const & M );
    Eigenvalues(
      integer         NRC,
      integer         nnz,
      real_type const values[],
      integer   const row[],
      integer   const col[]
    );

    void setup( integer NRC, real_type const data[], integer ldData );
    void setup( MatW const & M );

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
    Malloc<real_type> m_mem;

    integer     m_N;
    integer     m_Lwork;
    real_type * m_Re;
    real_type * m_Im;
    real_type * m_A_saved;
    real_type * m_VL;
    real_type * m_VR;
    real_type * m_Work;

    void allocate( integer N );
    void compute( );

  public:

    Eigenvectors();

    Eigenvectors( integer NRC, real_type const A[], integer ldA );

    explicit
    Eigenvectors( MatW const & A );

    Eigenvectors(
      integer         NRC,
      integer         A_nnz,
      real_type const A_values[],
      integer   const A_row[],
      integer   const A_col[]
    );

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
    Malloc<real_type> m_mem;

    integer     m_N;
    integer     m_Lwork;
    real_type * m_alphaRe;
    real_type * m_alphaIm;
    real_type * m_beta;
    real_type * m_Work;
    real_type * m_A_saved;
    real_type * m_B_saved;

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvalues();

    GeneralizedEigenvalues(
      integer NRC,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    );

    GeneralizedEigenvalues( MatW const & A, MatW const & B );

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
    );

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
    Malloc<real_type> m_mem_real;
    Malloc<integer>   m_mem_int;

    integer     m_N;
    integer     m_Lwork;
    integer     m_ilo;
    integer     m_ihi;
    real_type   m_abnorm;
    real_type   m_bbnorm;
    real_type * m_alphaRe;
    real_type * m_alphaIm;
    real_type * m_beta;
    real_type * m_A_saved;
    real_type * m_B_saved;
    real_type * m_VL;
    real_type * m_VR;
    real_type * m_lscale;
    real_type * m_rscale;
    real_type * m_rconde;
    real_type * m_rcondv;
    real_type * m_Work;
    integer   * m_iWork;
    integer   * m_bWork;

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvectors();

    GeneralizedEigenvectors(
      integer NRC,
      real_type const A[], integer ldA,
      real_type const B[], integer ldB
    );

    GeneralizedEigenvectors( MatW const & A, MatW const & B );

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
    );

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
