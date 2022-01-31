/*
  \file     sparse_tool_matrix_market.hh
  \date     2011, July 21
  \version  1.1
  \note     based on SparseLib 11.0 (C 1999 -- 2008)

  \author Enrico Bertolazzi
 
  Affiliations:
  Dipartimento di Ingegneria Industriale
  Universita` degli Studi di Trento
  enrico.bertolazzi@unitn.it
*/

#pragma once
#ifndef SPARSETOOL_MATRIX_MARKET_HH
#define SPARSETOOL_MATRIX_MARKET_HH

#include "sparse_tool.hh"

#include <Utils_zstr.hh>
#include <string.h>

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
  #define sscanf sscanf_s
#endif

namespace Sparse_tool {

  using ::std::ifstream;
  using ::std::ofstream;
  using ::std::string;

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

  //!
  //! Type of storage for sparse matrix.
  //!
  typedef enum { MM_COORDINATE = 0, MM_ARRAY = 1 } CoorType;

  //!
  //! Type of the element of the sparse matrix.
  //!
  typedef enum { MM_PATTERN = 0, MM_INTEGER = 1, MM_REAL = 2, MM_COMPLEX = 3 } real_type;

  //!
  //! The type of the matrix.
  //!
  typedef enum { MM_GENERAL = 0, MM_SYMMETRIC = 1, MM_SKEW_SYMMETRIC = 2, MM_HERMITIAN = 3 } MatrixType;
  
  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  static char const *cC[] = { "coordinate", "array" };
  static char const *cV[] = { "pattern", "integer", "real", "complex" };
  static char const *cM[] = { "general", "symmetric", "skew-symmetric", "Hermitian" };
  
  static char const MM_HEADER[] = "% Generated with saveToMatrixMarket of toolkit Sparse_tool\n"
                                  "% by Enrico Bertolazzi\n"
                                  "%--------------------------------------------------------\n";
  #endif

  //!
  //! Interface with Matrix Market exchange format
  //!
  class MatrixMarket {

    char    m_line[128];
    char    m_str[5][128];
    integer m_nrows;
    integer m_ncols;
    integer m_nnz;
    integer m_num_line;

    //! the type of matrix data: sparse or full
    CoorType   cType;
    //! the type of the data: ral integer or pattern
    real_type  vType;
    MatrixType mType; //! the matrix type

    bool
    getLine( istream_type & stream ) {
      stream.getline( m_line, 128 );
      ++m_num_line;
      return m_line[0] == '%';
    }

  public:

    //!
    //! Initialize an empty MatrixMarket object.
    //!
    MatrixMarket() { }

    //!
    //! Print to the stream object `s` 
    //! some information about the last loaded matrix. 
    //!
    void
    info( ostream_type & stream ) const {
      fmt::print( stream,
        "Rows:        {}\n"
        "Columns:     {}\n"
        "Non zeros:   {}\n"
        "ADDRESSING:  {}\n"
        "DATA:        {}\n"
        "MATRIX TYPE: {}\n",
        m_nrows, m_ncols, m_nnz,
        cC[cType], cV[vType], cM[mType]
      );
    }

    //!
    //! read the file in MM format from the opened stream `s`.
    //! The result is stored in the class instance
    //!
    void
    readHeader( istream_type & stream ) {

      UTILS_ASSERT0(
        stream.good(),
        "Sparse_tool::MatrixMarket::readHeader\n"
        "In reading Matrix Market File, cannot read header!\n"
      );

      m_num_line = 0;
      cType   = MM_COORDINATE;
      vType   = MM_REAL;
      mType   = MM_GENERAL;

      while ( this -> getLine( stream ) ) {
        #if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
          unsigned nstr = sscanf_s( m_line,
                                    "%s%s%s%s%s",
                                    m_str[0], unsigned(_countof(m_str[0])),
                                    m_str[1], unsigned(_countof(m_str[1])),
                                    m_str[2], unsigned(_countof(m_str[2])),
                                    m_str[3], unsigned(_countof(m_str[3])),
                                    m_str[4], unsigned(_countof(m_str[4])) );
        #else
          unsigned nstr = sscanf(m_line, "%s%s%s%s%s", m_str[0], m_str[1], m_str[2], m_str[3], m_str[4] );
        #endif

        if ( strcmp( m_str[0], "%%MatrixMarket" ) != 0 ||
             strcmp( m_str[1], "matrix" )         != 0 ) continue;

        // coordinate or array ?
        for ( unsigned i = 2; i < nstr; ++i ) {
          char const * p = m_str[i];
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
          
          else UTILS_ASSERT0( false, "Sparse_tool::MatrixMarket::readHeader\nbad MatrixMarket file format" );
        }
      }

      // check consistency of type
      if ( mType == MM_HERMITIAN ) {
        UTILS_ASSERT(
          vType == MM_COMPLEX,
          "Sparse_tool::MatrixMarket::readHeader\n"
          "inconsistent file data, a matrix with entries of type: {}\n"
          "can be a matrix of type: {}\n",
          cV[vType], cM[mType]
        );
      }
      if ( mType == MM_SKEW_SYMMETRIC ) {
        UTILS_ASSERT(
          vType == MM_REAL || vType == MM_COMPLEX,
          "Sparse_tool::MatrixMarket::readHeader\n"
          "inconsistent file data, a matrix with entries of type: {}\n"
          "can be a matrix of type: {}\n",
          cV[vType], cM[mType]
        );
      }

      switch ( cType ) {
      case MM_COORDINATE:
        sscanf( m_line, "%u%u%u", &m_nrows, &m_ncols, &m_nnz );
        break;
      case MM_ARRAY:
        sscanf( m_line, "%u%u", &m_nrows, &m_ncols );
        m_nnz = m_nrows * m_ncols;
        break;
      //default:
      //  break;
      }
    }

    void
    read( istream_type & stream, SparsePattern & sp ) {

      readHeader( stream );

      UTILS_ASSERT0(
        cType == MM_COORDINATE,
        "Sparse_tool::MatrixMarket::read\n"
        "file must be a coordinate file!"
      );

      sp.resize( m_nrows, m_ncols, m_nnz );

      integer kk = 0;
      while( kk < m_nnz ) {
        UTILS_ASSERT0(
          stream.good(), 
          "Sparse_tool::MatrixMarket::read\n"
          "Failed in reading Matrix Market File\n"
        );

        if ( getLine( stream ) ) continue; // skip comments and empty line

        integer i, j;
        sscanf( m_line, "%d%d", &i, &j );

        UTILS_ASSERT(
          i >= 1 && j >= 1 && i <= m_nrows && j <= m_ncols,
          "Sparse_tool::MatrixMarket::read\n"
          "In reading Matrix Market File, bad pattern ({},{}) index on line {}\n"
          "Read: <<{}>>\n",
          i, j, m_num_line, m_line
        );

        --i; --j; // zero base index
        sp.insert( i, j );
        if ( mType != MM_GENERAL ) sp.insert( j, i );
        ++kk; // next nnz
      }

      sp.internalOrder();
    }

    void
    read( istream_type & stream, CCoorMatrix<integer> & mat ) {

      readHeader( stream );

      UTILS_ASSERT0(
        cType == MM_COORDINATE,
        "Sparse_tool::MatrixMarket::read\n"
        "file must be a coordinate file!\n"
      );
      UTILS_ASSERT0(
        vType == MM_INTEGER,
        "Sparse_tool::MatrixMarket::read\n"
        "file must have integer entries!\n"
      );

      mat.resize( m_nrows, m_ncols, m_nnz );

      integer kk = 0;
      while ( kk < m_nnz ) {
        UTILS_ASSERT0(
          stream.good(), 
          "Sparse_tool::MatrixMarket::read\n"
          "Failed in reading Matrix Market File\n"
        );

        if ( getLine( stream ) ) continue; // skip comments

        integer i, j, a;
        sscanf( m_line, "%d%d%d", &i, &j, &a );

        UTILS_ASSERT(
          i >= 1 && j >= 1 && i <= m_nrows && j <= m_ncols,
          "Sparse_tool::MatrixMarket::read\n"
          "In reading Matrix Market File, bad pattern ({},{}) index on line {}\n"
          "Read: <<{}>>\n",
          i, j, m_num_line, m_line
        );

        --i; --j; // zero base index
        mat.insert(i,j) = a;
        switch ( mType ) {
          case MM_SYMMETRIC:      mat.insert(j,i) =  a; break;
          case MM_SKEW_SYMMETRIC: mat.insert(j,i) = -a; break;
          default: break;
        }
        ++kk;
      }

      mat.internalOrder();
    }

    template <typename T>
    void
    read( istream_type & stream, CCoorMatrix<std::complex<T> > & mat ) {

      readHeader( stream );

      static char const *fmts[4] = { "%d%d", "%d%d%d", "%d%d%lf", "%d%d%lf%lf" };

      UTILS_ASSERT0(
        cType == MM_COORDINATE,
        "Sparse_tool::MatrixMarket::read\n"
        "file must be a coordinate file!\n"
      );
      UTILS_ASSERT0(
        vType != MM_PATTERN,
        "Sparse_tool::MatrixMarket::read\n"
        "try to read a matrix from a pattern only file!\n"
      );

      mat.resize( m_nrows, m_ncols, m_nnz );

      integer kk = 0;
      while( kk < m_nnz ) {
        UTILS_ASSERT0(
          stream.good(), 
          "Sparse_tool::MatrixMarket::read\n"
          "Failed in reading Matrix Market File\n"
        );

        if ( getLine( stream ) ) continue; // skip comments

        integer i, j;
        double  re, im = 0;
        sscanf( m_line, fmts[vType], &i, &j, &re, &im );

        UTILS_ASSERT(
          i >= 1 && j >= 1 && i <= m_nrows && j <= m_ncols,
          "Sparse_tool::MatrixMarket::read\n"
          "In reading Matrix Market File, bad pattern ({},{}) index on line {}\n"
          "Read: <<{}>>\n",
          i, j, m_num_line, m_line
        );
        --i; --j; // zero base index

        mat.insert(i,j) = std::complex<T>(re,im);
        switch ( mType ) {
          case MM_SYMMETRIC:      mat.insert(j,i) = std::complex<T>(re,im);   break;
          case MM_SKEW_SYMMETRIC: mat.insert(j,i) = std::complex<T>(-re,-im); break;
          case MM_HERMITIAN:      mat.insert(j,i) = std::complex<T>(re,-im);  break;
          default: break;
        }
        ++kk;
      }

      mat.internalOrder();

    }

    template <typename T>
    void
    read( istream_type & stream, CCoorMatrix<T> & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_COORDINATE,
        "Sparse_tool::MatrixMarket::read\n"
        "file must be a coordinate file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_PATTERN,
        "Sparse_tool::MatrixMarket::read\n"
        "try to read a matrix from a pattern only file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_COMPLEX,
        "Sparse_tool::MatrixMarket::read\n"
        "try to read a complex matrix into a non complex one!, data is " << *this
      );

      mat.resize( m_nrows, m_ncols, m_nnz );

      integer kk = 0;
      while ( kk < m_nnz ) {

        UTILS_ASSERT0(
          stream.good(), 
          "Sparse_tool::MatrixMarket::read\n"
          "Failed in reading Matrix Market File\n"
        );

        if ( getLine( stream ) ) continue; // skip comments

        integer i, j;
        double  a;
        sscanf( m_line, "%d%d%lf", &i, &j, &a );

        UTILS_ASSERT(
          i >= 1 && j >= 1 && i <= m_nrows && j <= m_ncols,
          "Sparse_tool::MatrixMarket::read\n"
          "In reading Matrix Market File, bad pattern ({},{}) index on line {}\n"
          "Read: <<{}>>\n",
          i, j, m_num_line, m_line
        );
        --i; --j; // zero base index

        mat.insert(i,j) = a;
        switch ( mType ) {
          case MM_SYMMETRIC:      mat.insert(j,i) =  a; break;
          case MM_SKEW_SYMMETRIC: mat.insert(j,i) = -a; break;
          default: break;
        }
        ++kk;
      }

      mat.internalOrder();

    }

    void
    readFull( istream_type & stream, Vector<integer> & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_ARRAY,
        "Sparse_tool::MatrixMarket::readFull\n"
        "file must be an array file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType == MM_INTEGER,
        "Sparse_tool::MatrixMarket::readFull\n"
        "try to read an integer matrix into a non integer one!, data is " << *this
      );
      mat.resize( m_nnz );
      for ( integer kk = 0; kk < m_nnz; ++kk ) {
        int a;
        sscanf(m_line, "%d", &a);
        mat(kk) = a;
      }
    }

    template <typename T>
    void
    readFull( istream_type & stream, Vector<std::complex<T> > & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_ARRAY,
        "Sparse_tool::MatrixMarket::readFull\n"
        "file must be an array file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_PATTERN,
        "Sparse_tool::MatrixMarket::readFull\n"
        "try to read a pattern from a full matrix!, data is " << *this
      );
      static char const *fmts[4] = { "", "%lf", "%lf", "%lf%lf" };
      mat.clear(); mat.resize( m_nnz );
      for ( integer kk = 0; kk < m_nnz; ++kk ) {
        double re, im = 0;
        sscanf(m_line, fmts[mType], &re, &im );
        mat.push_back( std::complex<T>(re,im) );
      }
    }

    template <typename T>
    void
    readFull( istream_type & stream, Vector<T> & mat ) {

      readHeader( stream );

      SPARSETOOL_ASSERT(
        cType == MM_ARRAY,
        "Sparse_tool::MatrixMarket::readFull\n"
        "file must be an array file!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_PATTERN,
        "Sparse_tool::MatrixMarket::readFull\n"
        "try to read a pattern from a full matrix!, data is " << *this
      );
      SPARSETOOL_ASSERT(
        vType != MM_COMPLEX,
        "Sparse_tool::MatrixMarket::readFull\n"
        "try to read a a complex matrix into a non complex one!, data is " << *this
      );
      //static char const *fmts[4] = { "", "%lf", "%lf", "%lf%lf" };
      mat.clear(); mat.resize( m_nnz );
      for ( integer kk = 0; kk < m_nnz; ++kk ) {
        double a;
        sscanf(m_line, "%lf", &a );
        mat.push_back( a );
      }
    }

    //!
    //! Read the matrix in the template `MAT` class.
    //!
    template<typename MAT>
    void
    read( istream_type & stream, MAT & M ) {
      CCoorMatrix<typename MAT::real_type> M1;
      read( stream, M1 );
      M.resize( M1 );
    }

    //!
    //! Read the matrix in the template `MAT` class.
    //!
    template<typename MAT>
    void
    read( char const fname[], MAT & M ) {
      integer len = strlen(fname);
      string  msg = fmt::format(
        "Sparse_tool::MatrixMarket::read\n"
        "In reading Matrix Market File, cannot open {}\n",
        fname
      );
      if ( len > 4 && strcmp(fname+len-3,".gz") == 0 ) {
        ifstream file(fname,std::ios::binary);
        UTILS_ASSERT( file.is_open(), msg );
        zstr::istream fz(file);
        read( fz, M );
        file.close();
      } else {
        ifstream file(fname);
        UTILS_ASSERT( file.is_open(), msg );
        read( file, M );
        file.close();
      }
    }

    //!
    //! Read the matrix in the template `MAT` class.
    //!
    template<typename MAT>
    void
    read( string const & fname, MAT & M ) {
      read( fname.c_str(), M );
    }

    ///////////////////////////////////////////////////////////////////

    integer numRows() const { return m_nrows; } //!< number of rows of loaded matrix
    integer numCols() const { return m_ncols; } //!< number of columns of loaded matrix
    integer nnz()     const { return m_nnz;   } //!< total number of nonzeros

    //!
    //! Type of storage:
    //! - `0` by coordinate
    //! - `1` full matrix in column major order
    //!
    CoorType coor_type() const { return cType; }
    //! 
    //! Value type:
    //! - `0` no values, only pattern,
    //! - `1` type `int`
    //! - `2` type `double`
    //! - `3` type `complex<double>`
    //!
    real_type value_type() const { return vType; }
    //! 
    //! Matrix type:
    //! - `0` general matrix
    //! - `1` symmetric matrix
    //! - `2` skew symmetrix matrix
    //! - `3` Hermitian matrix
    //!
    MatrixType matrix_type() const { return mType; }

    //!
    //! Print to the stream object `stream` some
    //! information about the last loaded matrix.
    //!
    friend ostream_type & operator << ( ostream_type & stream, MatrixMarket const & mm )
    { mm.info( stream ); return stream; }
  };

  //!
  //! Save a matrix to a file in MatrixMarket format.
  //!
  //! \param fname the name of the file to save
  //! \param A     sparse matrix to save
  //! \param vType \copydoc Sparse_tool::CoorType
  //! \param mType \copydoc Sparse_tool::MatrixType
  //!
  template <typename T,typename M>
  static
  void
  MatrixMarketSaveToFile(
    std::string               const & fname,
    Sparse<std::complex<T>,M> const & A,
    real_type                         vType,
    MatrixType                        mType
  ) {

    std::ofstream file;
    file.open( fname.c_str() );
    file.precision(15);
    fmt::print( file,
      "%%MatrixMarket matrix coordinate {} {}\n"
      "{}"
      "{} {} {}\n",
      cV[vType], cM[mType], MM_HEADER,
      A.numRows(), A.numCols(), A.nnz()
    );

    switch ( vType ) {
      case MM_PATTERN:
        for ( A.Begin(); A.End(); A.Next() )
          fmt::print( file, "{}\t{}\n", A.row()+1, A.column()+1 );
        break;
      case MM_INTEGER:
      case MM_REAL:
        UTILS_ERROR("Sparse_tool::MatrixMarketSaveToFile: Cannot write a complex matrix with non complex type");
        break;
      case MM_COMPLEX:
        for ( A.Begin(); A.End(); A.Next() )
          fmt::print(
            file, "{}\t{}\t{}\t{}\n",
            A.row()+1, A.column()+1, A.value().real(), A.value().imag()
          );
        break;
    //  default:
    //    break;
    }

    file.close();
  }

  //! 
  //! Save a matrix to a file in MatrixMarket format.
  //!
  //! \param fname the name of the file to save
  //! \param A     sparse matrix to save
  //! \param vType \copydoc Sparse_tool::CoorType
  //! \param mType \copydoc Sparse_tool::MatrixType
  //! 
  template <typename T,typename M>
  static
  void
  MatrixMarketSaveToFile(
    std::string const & fname,
    Sparse<T,M> const & A,
    real_type           vType,
    MatrixType          mType
  ) {

    std::ofstream file;
    file.open( fname.c_str() );
    file.precision(15);

    fmt::print( file,
      "%%MatrixMarket matrix coordinate {} {}\n"
      "{}"
      "{} {} {}\n",
      cV[vType], cM[mType], MM_HEADER,
      A.numRows(), A.numCols(), A.nnz()
    );

    switch ( vType ) {
      case MM_PATTERN:
        for ( A.Begin(); A.End(); A.Next() )
          fmt::print( file, "{}\t{}\n", A.row()+1, A.column()+1 );
        break;
      case MM_INTEGER:
      case MM_REAL:
        for ( A.Begin(); A.End(); A.Next() )
          fmt::print( file, "{}\t{}\t{}\n", A.row()+1, A.column()+1, A.value() );
        break;
      case MM_COMPLEX:
        UTILS_ERROR("Sparse_tool::MatrixMarketSaveToFile: Cannot write a NON complex matrix with a complex type"); 
        break;
    //  default:
    //    break;
    }

    file.close();
  }
    
  //! 
  //! Save a vector to a file in MatrixMarket format.
  //!
  //! \param fname the name of the file to save
  //! \param V     vector to save
  //! 
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

    fmt::print( file,
      "%%MatrixMarket matrix array complex general\n"
      "{}"
      "{} 1\n",
      MM_HEADER, V.size()
    );
    for ( integer i = 0; i < V.size(); ++i )
      fmt::print( file, "{} {}\n", V(i).real(), V(i).imag() );

    file.close();
  }

  //! 
  //! Save a vector to a file in MatrixMarket format.
  //!
  //! \param fname the name of the file to save
  //! \param V     vector to save
  //! 
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

    fmt::print( file,
      "%%MatrixMarket matrix array real general\n"
      "{}"
      "{} 1\n",
      MM_HEADER, V.size()
    );
    for ( integer i = 0; i < V.size(); ++i )
      fmt::print( file, "{}\n", V(i) );

    file.close();
  }

  //! 
  //! Save a pattern to a file in MatrixMarket format.
  //!
  //! \param fname the name of the file to save
  //! \param A     sparse matrix to save
  //! \param mType \copydoc Sparse_tool::MatrixType
  //! 
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
    fmt::print( file,
      "%%MatrixMarket matrix coordinate pattern {}\n"
      "{}"
      "{} {} {}\n",
      cM[mType], MM_HEADER,
      A.numRows(), A.numCols(), A.nnz()
    );
    for ( A.Begin(); A.End(); A.Next() )
      fmt::print( file, "{}\t{}\n", A.row()+1, A.column()+1 );
    file.close();
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {

  using ::Sparse_tool::MatrixMarket;
  using ::Sparse_tool::MatrixMarketSaveToFile;

  using ::Sparse_tool::MM_COORDINATE;
  using ::Sparse_tool::MM_ARRAY;
  using ::Sparse_tool::MM_PATTERN;
  using ::Sparse_tool::MM_INTEGER;
  using ::Sparse_tool::MM_REAL;
  using ::Sparse_tool::MM_COMPLEX;
  using ::Sparse_tool::MM_GENERAL;
  using ::Sparse_tool::MM_SYMMETRIC;
  using ::Sparse_tool::MM_SKEW_SYMMETRIC;
  using ::Sparse_tool::MM_HERMITIAN;

}
#endif

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
