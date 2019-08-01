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
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;

    integer     N, Lwork;
    valueType * Re;
    valueType * Im;
    valueType * Work;
    valueType * A_saved;

    void allocate( integer N );
    void compute( );

  public:

    Eigenvalues();
    Eigenvalues( integer NRC, valueType const data[], integer ldData );
    Eigenvalues( MatW const & M );
    Eigenvalues(
      integer         NRC,
      integer         nnz,
      valueType const values[],
      integer   const row[],
      integer   const col[]
    );

    void setup( integer NRC, valueType const data[], integer ldData );
    void setup( MatW const & M );

    void
    setup(
      integer         NRC,
      integer         nnz,
      valueType const values[],
      integer   const row[],
      integer   const col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class Eigenvectors {
  public:
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;

    integer     N, Lwork;
    valueType * Re;
    valueType * Im;
    valueType * A_saved;
    valueType * VL;
    valueType * VR;
    valueType * Work;

    void allocate( integer N );
    void compute( );

  public:

    Eigenvectors();

    Eigenvectors( integer NRC, valueType const A[], integer ldA );

    Eigenvectors( MatW const & A );

    Eigenvectors(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[]
    );

    void
    setup( integer NRC, valueType const A[], integer ldA );

    void setup( MatW const & A );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

    void getLeftEigenvector( std::vector<std::vector<complexType> > & vecs ) const;
    void getRightEigenvector( std::vector<std::vector<complexType> > & vecs ) const;
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class GeneralizedEigenvalues {
  public:
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;

    integer     N, Lwork;
    valueType * alphaRe;
    valueType * alphaIm;
    valueType * beta;
    valueType * Work;
    valueType * A_saved;
    valueType * B_saved;

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvalues();

    GeneralizedEigenvalues(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    GeneralizedEigenvalues( MatW const & A, MatW const & B );

    GeneralizedEigenvalues(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void
    setup(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  class GeneralizedEigenvectors {
  public:
    typedef T                valueType;
    typedef std::complex<T>  complexType;
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

  private:
    Malloc<valueType> mem_real;
    Malloc<integer>   mem_int;

    integer     N, Lwork, ilo, ihi;
    valueType   abnorm, bbnorm;
    valueType * alphaRe;
    valueType * alphaIm;
    valueType * beta;
    valueType * A_saved;
    valueType * B_saved;
    valueType * VL;
    valueType * VR;
    valueType * lscale;
    valueType * rscale;
    valueType * rconde;
    valueType * rcondv;
    valueType * Work;
    integer   * iWork;
    integer   * bWork;

    void allocate( integer N );
    void compute( );

  public:

    GeneralizedEigenvectors();

    GeneralizedEigenvectors(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    GeneralizedEigenvectors( MatW const & A, MatW const & B );

    GeneralizedEigenvectors(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void
    setup(
      integer NRC,
      valueType const A[], integer ldA,
      valueType const B[], integer ldB
    );

    void setup( MatW const & A, MatW const & B );

    void
    setup(
      integer         NRC,
      integer         A_nnz,
      valueType const A_values[],
      integer   const A_row[],
      integer   const A_col[],
      integer         B_nnz,
      valueType const B_values[],
      integer   const B_row[],
      integer   const B_col[]
    );

    void getEigenvalue( integer n, valueType & re, valueType & im ) const;
    void getEigenvalue( integer n, complexType & eig ) const;

    void getEigenvalues( std::vector<valueType> & re, std::vector<valueType> & im ) const;
    void getEigenvalues( std::vector<complexType> & eigs ) const;

    void getLeftEigenvector( std::vector<std::vector<complexType> > & vecs ) const;
    void getRightEigenvector( std::vector<std::vector<complexType> > & vecs ) const;

    valueType balancedAnorm1() const { return this->abnorm; }
    valueType balancedBnorm1() const { return this->bbnorm; }
    valueType const * getLscale() const { return this->lscale; }
    valueType const * getRscale() const { return this->rscale; }
    valueType const * RcondEigenvalues()  const { return this->rconde; }
    valueType const * RcondEigenvectors() const { return this->rcondv; }
  };

}

///
/// eof: eig.hxx
///
