/*!

  \file     SparseToolMatrixMarket.hh
  \date     2011, July 21
  \version  1.1
  \note     based on SparseLib 11.0 (C 1999 -- 2008)

  \author Enrico Bertolazzi
 
  \par Affiliations:
       Dipartimento di Ingegneria Industriale<br>
       Universita` degli Studi di Trento<br>
       enrico.bertolazzi\@unitn.it
 
*/

#ifndef SPARSETOOL_MATRIX_MARKET_HH
#define SPARSETOOL_MATRIX_MARKET_HH

#include <string.h>

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
  #define sscanf sscanf_s
#endif

namespace SparseTool {

  using ::std::ifstream;
  using ::std::ofstream;

  /*
  //  #     #                              
  //  ##   ##   ##   ##### #####  # #    # 
  //  # # # #  #  #    #   #    # #  #  #  
  //  #  #  # #    #   #   #    # #   ##   
  //  #     # ######   #   #####  #   ##   
  //  #     # #    #   #   #   #  #  #  #  
  //  #     # #    #   #   #    # # #    # 
  //                                       
  //  #     #                                   
  //  ##   ##   ##   #####  #    # ###### ##### 
  //  # # # #  #  #  #    # #   #  #        #   
  //  #  #  # #    # #    # ####   #####    #   
  //  #     # ###### #####  #  #   #        #   
  //  #     # #    # #   #  #   #  #        #   
  //  #     # #    # #    # #    # ######   #   
  */

  //! type of storage for sparse matrix
  typedef enum { MM_COORDINATE = 0, MM_ARRAY = 1 } CoorType;

  //! type of the element of the sparse matrix
  typedef enum { MM_PATTERN = 0, MM_INTEGER = 1, MM_REAL = 2, MM_COMPLEX = 3 } ValueType;

  //! the type of the matrix
  typedef enum { MM_GENERAL = 0, MM_SYMMETRIC = 1, MM_SKEW_SYMMETRIC = 2, MM_HERMITIAN = 3 } MatrixType;
  
  /*! \cond NODOC */

  static char const *cC[] = { "coordinate", "array" };
  static char const *cV[] = { "pattern", "integer", "real", "complex" };
  static char const *cM[] = { "general", "symmetric", "skew-symmetric", "Hermitian" };
  
  static char const MM_HEADER[] = "% Generated with saveToMatrixMarket of toolkit SparseTool\n"
                                  "% by Enrico Bertolazzi\n"
                                  "%--------------------------------------------------------\n";
  /*! \endcond */

  /*!
     A standard way to store and exchange large sparse matrix if to save
     matrix accordingly to some <a> file exchange format </a>.
     For sparse matrices there are two widely used format the Harwell-Boeing
     (HB) Sparse Matrix Format and Matrix Market (MM) Sparse Matrix Format.
     
     The HB format is strongly depended on \c FORTRAN language and is 
     difficult to manage in other languages unless using a sophisticated parser.
     The alternative is to use dedicated \c FORTRAN routine to manage such a format.
     The design of \c SparseTool is to use a unique header \c SparseTool.hh
     which contains the whole library so that HB format is not supported.
     The MM format is easier to manage so that a simple support is included in 
     the toolkit by the class \c MatrixMarket.
     In any case there are free software for convert from one format to another
     (see e.g. \c http://bebop.cs.berkeley.edu/ ).

     The following code shows how \c MatrixMarket should be used:

\code
MatrixMarket mm; // define the object mm to manage Matrix Market file
CCoorMatrix<double> A;      // an empty sparse CCOOR matrix
SparsePattern sp;           // an empty sparse pattern
mm.read("hor__131.mtx");  // read matrix and store in mm object
mm.load_pattern( sp );    // extract the pattern
mm.load_matrix( A );      // copy loaded matrix in sparse matrix A
\endcode

   The class \c MatrixMarket is not a \c template class 
   because the value type of the nonzeros is specified by the format.
   The available type in MM format are \c int, \c double
   and <c> complex<double> </c>. 
   The representation of sparse matrices in MM format is of type CCOOR
   so that \c MatrixMarket store in this form the loaded matrix.
   Full matrices are stored in column major order.
   Sometimes it can be useful to access directly the structure of the 
   loaded matrix with the following methods:

   In the library it is not provided a \c write method to MM
   format. This choice is due to the fact that write a MM file is very easy 
   using the iterators and the only complication is to write 
   the first few rows of the header file.
  
   */

  //! Interface with Matrix Market exchange format
  class MatrixMarket {

    char       line[128], str[5][128];
    unsigned   nRows, nCols, numNnz, numLine;

    CoorType   cType; //! the type of matrix data: sparse or full
    ValueType  vType; //! the type of the data: ral integer or pattern
    MatrixType mType; //! the matrix type

    bool
    getLine( istream & stream ) {
      stream.getline( line, 128 );
      ++numLine;
      return line[0] == '%';
    }

  public:

    //! initialize an empty MatrixMarket object
    MatrixMarket() { }

    //! print to the stream object \c s some information about the last loaded matrix. 
    void
    info( ostream & stream ) const {
      stream << "\nRows:        " << nRows
             << "\nColumns:     " << nCols
             << "\nNon zeros:   " << numNnz
             << "\nADDRESSING:  " << cC[cType]
             << "\nDATA:        " << cV[vType]
             << "\nMATRIX TYPE: " << cM[mType]
             << "\n\n";
    }
    
    /*! \brief
     *  read the file in MM format from the opened stream \c s.
     *  The result is stored in the class instance
     */
    void
    readHeader( istream & stream ) {

      SPARSETOOL_ASSERT(
        stream.good(),
        "In reading Matrix Market File, cannot read header!"
      );

      numLine = 0;
      cType   = MM_COORDINATE;
      vType   = MM_REAL;
      mType   = MM_GENERAL;

      while ( this -> getLine( stream ) ) {
        #if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
          unsigned nstr = sscanf_s( line,
                                    "%s%s%s%s%s",
                                    str[0], unsigned(_countof(str[0])),
                                    str[1], unsigned(_countof(str[1])),
                                    str[2], unsigned(_countof(str[2])),
                                    str[3], unsigned(_countof(str[3])),
                                    str[4], unsigned(_countof(str[4])) );
        #else
          unsigned nstr = sscanf(line, "%s%s%s%s%s", str[0], str[1], str[2], str[3], str[4] );
        #endif

        if ( strcmp( str[0], "%%MatrixMarket" ) != 0 ||
             strcmp( str[1], "matrix" )         != 0 ) continue;

        // coordinate or array ?
        for ( unsigned i = 2; i < nstr; ++i ) {
          char const * p = str[i];
          if      ( strcmp( p, "coordinate")     == 0 ) cType = MM_COORDINATE;
          else if ( strcmp( p, "array")          == 0 ) cType = MM_ARRAY;   

          else if ( strcmp( p, "pattern")        == 0 ) vType = MM_PATTERN;
          else if ( strcmp( p, "integer")        == 0 ) vType = MM_INTEGER;
          else if ( strcmp( p, "real")           == 0 ) vType = MM_REAL;        
          else if ( strcmp( p, "complex")        == 0 ) vType = MM_COMPLEX;        

          else if ( strcmp( p, "general")        == 0 ) mType = MM_GENERAL;
          else if ( strcmp( p, "symmetric")      == 0 ) mType = MM_SYMMETRIC;
          else if ( strcmp( p, "skew-symmetric") == 0 ) mType = MM_SKEW_SYMMETRIC;
          else if ( strcmp( p, "Hermitian")      == 0 ) mType = MM_HERMITIAN;
          
          else SPARSETOOL_ASSERT( false, "bad MatrixMarket file format" );
        }
      }
      
      // check consistency of type
      if ( mType == MM_HERMITIAN ) {
        SPARSETOOL_ASSERT(
          vType == MM_COMPLEX,
          "inconsistent file data:\n" <<
          "a matrix with entries of type: " << cV[vType] <<
          "can be a matrix of type: " << cM[mType]
        );
      }
      if ( mType == MM_SKEW_SYMMETRIC ) {
        SPARSETOOL_ASSERT(
          vType == MM_REAL || vType == MM_COMPLEX,
          "inconsistent file data:\n" <<
          "a matrix with entries of type: " << cV[vType] <<
          "can be a matrix of type: " << cM[mType]
        );
      }

      switch ( cType ) {
      case MM_COORDINATE:
        sscanf( line, "%u%u%u", &nRows, &nCols, &numNnz );
        break;
      case MM_ARRAY:
        sscanf( line, "%u%u", &nRows, &nCols );
        numNnz = nRows * nCols;
        break;
      default:
        break;
      }
    }

    void
    read( istream & stream, SparsePattern & sp ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_COORDINATE,
        "MatrixMarket: file must be a coordinate file!"
      );

      sp.resize( nRows, nCols, numNnz );

      for ( unsigned kk = 0; kk < numNnz; ++kk ) {
        while ( getLine( stream ) ); // skip comments
        unsigned i, j;
        sscanf( line, "%u%u", &i, &j );
        --i; --j; // zero base index

        SPARSETOOL_ASSERT(
          i <= nRows && j <= nCols,
          "In reading Matrix Market File, bad pattern index on line " <<
          numLine << "\nRead: <<" << line << ">>"
        );

        sp.insert( i, j );
        if ( mType != MM_GENERAL ) sp.insert( j, i );
      }

      sp.internalOrder();

    }

    void
    read( istream & stream, CCoorMatrix<indexType> & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_COORDINATE,
        "MatrixMarket: file must be a coordinate file!"
      );
      SPARSETOOL_ASSERT(
        vType == MM_INTEGER,
        "MatrixMarket: file must have integer entries!"
      );

      mat.resize( nRows, nCols, numNnz );

      for ( unsigned kk = 0; kk < numNnz; ++kk ) {
        while ( getLine( stream ) ); // skip comments
        unsigned i, j;
        int      a;
        sscanf( line, "%u%u%d", &i, &j, &a );
        --i; --j; // zero base index

        SPARSETOOL_ASSERT(
          i <= nRows && j <= nCols,
          "In reading Matrix Market File, bad pattern index on line " <<
          numLine << "\nRead<<" << line << ">>"
        );

        mat.insert(i,j) = a;
        switch ( mType ) {
          case MM_SYMMETRIC:      mat.insert(j,i) =  a; break;
          case MM_SKEW_SYMMETRIC: mat.insert(j,i) = -a; break;
          default: break;
        }
      }
      
      mat.internalOrder();

    }

    template <typename T>
    void
    read( istream & stream, CCoorMatrix<std::complex<T> > & mat ) {

      readHeader( stream );

      static char const *fmts[4] = { "%u%u", "%u%u%d", "%u%u%lf", "%u%u%lf%lf" };

      SPARSETOOL_ASSERT(
        cType == MM_COORDINATE,
        "MatrixMarket: file must be a coordinate file!"
      );
      SPARSETOOL_ASSERT(
        vType != MM_PATTERN,
        "MatrixMarket: try to read a matrix from a pattern only file!"
      );

      mat.resize( nRows, nCols, numNnz );

      for ( unsigned kk = 0; kk < numNnz; ++kk ) {
        while ( getLine( stream ) ); // skip comments
        unsigned i, j;
        double   re, im = 0;
        sscanf( line, fmts[vType], &i, &j, &re, &im );
        --i; --j; // zero base index

        SPARSETOOL_ASSERT(
          i <= nRows && j <= nCols,
          "In reading Matrix Market File, bad pattern index on line " <<
          numLine << "\nRead<<" << line << ">>"
        );

        mat.insert(i,j) = std::complex<T>(re,im);
        switch ( mType ) {
          case MM_SYMMETRIC:      mat.insert(j,i) = std::complex<T>(re,im);   break;
          case MM_SKEW_SYMMETRIC: mat.insert(j,i) = std::complex<T>(-re,-im); break;
          case MM_HERMITIAN:      mat.insert(j,i) = std::complex<T>(re,-im);  break;
          default: break;
        }
      }

      mat.internalOrder();

    }

    template <typename T>
    void
    read( istream & stream, CCoorMatrix<T> & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_COORDINATE,
        "MatrixMarket: file must be a coordinate file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_PATTERN,
        "MatrixMarket: try to read a matrix from a pattern only file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_COMPLEX,
        "MatrixMarket: try to read a complex matrix into a non complex one!, data is " << *this
      );

      mat.resize( nRows, nCols, numNnz );

      for ( unsigned kk = 0; kk < numNnz; ++kk ) {
        while ( getLine( stream ) ); // skip comments
        unsigned i, j;
        double   a;
        sscanf( line, "%u%u%lf", &i, &j, &a );
        --i; --j; // zero base index

        SPARSETOOL_ASSERT(
          i <= nRows && j <= nCols,
          "In reading Matrix Market File, bad pattern index on line " <<
          numLine << "\nRead<<" << line << ">>, data is " << *this
        );

        mat.insert(i,j) = a;
        switch ( mType ) {
          case MM_SYMMETRIC:      mat.insert(j,i) =  a; break;
          case MM_SKEW_SYMMETRIC: mat.insert(j,i) = -a; break;
          default: break;
        }
      }

      mat.internalOrder();

    }

    void
    readFull( istream & stream, Vector<indexType> & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_ARRAY,
        "MatrixMarket: file must be an array file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType == MM_INTEGER,
        "MatrixMarket: try to read an integer matrix into a non integer one!, data is " << *this
      );
      mat.clear(); mat.resize( numNnz );
      for ( unsigned kk = 0; kk < numNnz; ++kk ) {
        int a;
        sscanf(line, "%d", &a);
        mat.push_back( a );
      }
    }

    template <typename T>
    void
    readFull( istream & stream, Vector<std::complex<T> > & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_ARRAY,
        "MatrixMarket: file must be an array file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_PATTERN,
        "MatrixMarket: try to read a pattern from a full matrix!, data is " << *this
      );
      static char const *fmts[4] = { "", "%lf", "%lf", "%lf%lf" };
      mat.clear(); mat.resize( numNnz );
      for ( unsigned kk = 0; kk < numNnz; ++kk ) {
        double re, im = 0;
        sscanf(line, fmts[mType], &re, &im );
        mat.push_back( std::complex<T>(re,im) );
      }
    }

    template <typename T>
    void
    readFull( istream & stream, Vector<T> & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_ARRAY,
        "MatrixMarket: file must be an array file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_PATTERN,
        "MatrixMarket: try to read a pattern from a full matrix!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_COMPLEX,
        "MatrixMarket: try to read a a complex matrix into a non complex one!, data is " << *this
      );
      //static char const *fmts[4] = { "", "%lf", "%lf", "%lf%lf" };
      mat.clear(); mat.resize( numNnz );
      for ( unsigned kk = 0; kk < numNnz; ++kk ) {
        double a;
        sscanf(line, "%lf", &a );
        mat.push_back( a );
      }
    }

    //! read the red matrix in the template \c M1 class.
    template<typename MAT>
    void
    read( istream & stream, MAT & M ) {
      CCoorMatrix<typename MAT::valueType> M1;
      read( stream, M1 );
      M.resize( M1 );
    }
        
    //! read the red matrix in the template \c M1 class.
    template<typename MAT>
    void
    read( char const fname[], MAT & M ) {
      ifstream file(fname);
      SPARSETOOL_ASSERT(
        file.is_open(),
        "In reading Matrix Market File, cannot open " << fname << '\n'
      );
      read( file, M );
      file.close();
    }

    ///////////////////////////////////////////////////////////////////

    unsigned numRows() const { return nRows;  } //!< number of rows of loaded matrix
    unsigned numCols() const { return nCols;  } //!< number of columns of loaded matrix
    unsigned nnz()     const { return numNnz; } //!< total number of nonzeros

    /*! \brief
     *  type of storage:
     *  - \c 0 by coordinate
     *  - \c 1 full matrix in column major order
     */
    CoorType coor_type() const { return cType; }
    /*! \brief 
     *  value type:
     *  - \c 0 no values, only pattern,
     *  - \c 1 type \c int
     *  - \c 2 type \c double
     *  - \c 3 type <c> complex<double> </c>
     */
    ValueType value_type() const { return vType; }
    /*! \brief 
     *  matrix type:
     *  - \c 0 general matrix
     *  - \c 1 symmetric matrix
     *  - \c 2 skew symmetrix matrix
     *  - \c 3 Hermitian matrix
     */
    MatrixType matrix_type() const { return mType; }

    //! print to the stream object \c stream some information about the last loaded matrix. 
    friend ostream & operator << ( ostream & stream, MatrixMarket const & mm )
    { mm.info( stream ); return stream; }
  };

  /*!
   *  Save a matrix to a file in MatrixMarket format
   *  \param fname the name of the file to save
   *  \param A     sparse matrix to save
   *  \param vType \copydoc SparseTool::CoorType
   *  \param mType \copydoc SparseTool::MatrixType
   */

  template <typename T,typename M>
  static
  void
  MatrixMarketSaveToFile(
    std::string               const & fname,
    Sparse<std::complex<T>,M> const & A,
    ValueType                         vType,
    MatrixType                        mType
  ) {

    std::ofstream file;
    file.open( fname.c_str() );
    file.precision(15);
    file
      << "%%MatrixMarket matrix coordinate "
      << cV[vType]   << ' '
      << cM[mType]   << '\n'
      << MM_HEADER
      << A.numRows() << ' '
      << A.numCols() << ' '
      << A.nnz()     << '\n';

    switch ( vType ) {
      case MM_PATTERN:
        for ( A.Begin(); A.End(); A.Next() )
          file
            << A.row()+1    << '\t'
            << A.column()+1 << '\n';
        break;
      case MM_INTEGER:
      case MM_REAL:
        SPARSETOOL_ERR("Cannot write a complex matrix with non complex type"); 
        break;
      case MM_COMPLEX:
        for ( A.Begin(); A.End(); A.Next() )
          file
            << A.row()+1        << '\t'
            << A.column()+1     << '\t'
            << A.value().real() << '\t'
            << A.value().imag() << '\n';
        break;
      default:
        break;
    }

    file.close();
  }

  /*!
   * Save a matrix to a file in MatrixMarket format
   * \param fname the name of the file to save
   * \param A     sparse matrix to save
   * \param vType \copydoc SparseTool::CoorType
   * \param mType \copydoc SparseTool::MatrixType
   */

  template <typename T,typename M>
  static
  void
  MatrixMarketSaveToFile(
    std::string const & fname,
    Sparse<T,M> const & A,
    ValueType           vType,
    MatrixType          mType
  ) {

    std::ofstream file;
    file.open( fname.c_str() );
    file.precision(15);

    file
      << "%%MatrixMarket matrix coordinate "
      << cV[vType]   << ' '
      << cM[mType]   << '\n'
      << MM_HEADER
      << A.numRows() << ' '
      << A.numCols() << ' '
      << A.nnz()     << '\n';

    switch ( vType ) {
      case MM_PATTERN:
        for ( A.Begin(); A.End(); A.Next() )
          file
            << A.row()+1    << '\t'
            << A.column()+1 << '\n';
        break;
      case MM_INTEGER:
      case MM_REAL:
        for ( A.Begin(); A.End(); A.Next() )
          file
            << A.row()+1    << '\t'
            << A.column()+1 << '\t'
            << A.value()    << '\n';
        break;
      case MM_COMPLEX:
        SPARSETOOL_ERR("Cannot write a NON complex matrix with a complex type"); 
        break;
      default:
        break;
    }

    file.close();
  }
    
  /*!
   * Save a vector to a file in MatrixMarket format
   * \param fname the name of the file to save
   * \param V     vector to save
   */

  template <typename T>
  static
  void
  MatrixMarketSaveToFile(
    std::string              const & fname,
    Vector<std::complex<T> > const & V
  ) {

    std::ofstream file;
    file.open( fname.c_str() );
    file.precision(15);

    file
      << "%%MatrixMarket matrix array complex general\n"
      << MM_HEADER
      << V.size() << " 1\n";
    for ( indexType i = 0; i < V.size(); ++i )
      file << V(i).real() << ' ' << V(i).imag() << '\n';

    file.close();
  }

  /*!
   * Save a vector to a file in MatrixMarket format
   * \param fname the name of the file to save
   * \param V     vector to save
   */

  template <typename T>
  static
  void
  MatrixMarketSaveToFile(
    std::string const & fname,
    Vector<T>   const & V
  ) {

    std::ofstream file;
    file.open( fname.c_str() );
    file.precision(15);

    file << "%%MatrixMarket matrix array real general\n"
         << MM_HEADER
         << V.size() << " 1\n";
    for ( indexType i = 0; i < V.size(); ++i )
      file << V(i) << '\n';

    file.close();
  }

  /*!
   * Save a pattern to a file in MatrixMarket format
   * \param fname the name of the file to save
   * \param A     sparse matrix to save
   * \param mType \copydoc SparseTool::MatrixType
   */

  static
  inline
  void
  MatrixMarketSaveToFile(
    std::string   const & fname,
    SparsePattern const & A,
    MatrixType    const   mType
  ) {
    std::ofstream file;
    file.open( fname.c_str() );
    file.precision(15);
    file
      << "%%MatrixMarket matrix coordinate pattern "
      << cM[mType]   << '\n'
      << MM_HEADER
      << A.numRows() << ' '
      << A.numCols() << ' '
      << A.nnz()     << '\n';
    for ( A.Begin(); A.End(); A.Next() )
      file << A.row()+1 << '\t' << A.column()+1 << '\n';
    file.close();
  }

}

namespace SparseToolLoad {

  using ::SparseTool::MatrixMarket;
  using ::SparseTool::MatrixMarketSaveToFile;

  using ::SparseTool::MM_COORDINATE;
  using ::SparseTool::MM_ARRAY;
  using ::SparseTool::MM_PATTERN;
  using ::SparseTool::MM_INTEGER;
  using ::SparseTool::MM_REAL;
  using ::SparseTool::MM_COMPLEX;
  using ::SparseTool::MM_GENERAL;
  using ::SparseTool::MM_SYMMETRIC;
  using ::SparseTool::MM_SKEW_SYMMETRIC;
  using ::SparseTool::MM_HERMITIAN;

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
