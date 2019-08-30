/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MA41     |
 |                                                                          |
 |  date         : 2011, 18 July                                            |
 |  version      : 1.0.                                                     |
 |  file         : MA41.hxx                                                 |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_MA48_HH
#define SPARSETOOL_MA48_HH

#include "../sparse_tool.hh"
#include "../../HSL/hsl.h"

namespace SparseTool {

  template <typename T>
  class MA48 {
  public:
    typedef T            real_type;
    typedef HSL::integer integer;

  private:

    bool      verbose;
    integer   nRow, nCol, nnz;
    integer   ICNTL[20];
    integer   INFO[20];
    real_type CNTL[10];
    real_type RINFO[10];

    // the matrix
    vector<real_type> A;
    vector<integer>   IRN, JCN;

    // work arrays
    vector<real_type> A_work;
    vector<integer>   IRN_work, JCN_work;
    vector<real_type> W;
    vector<integer>   IW, KEEP;

    void
    msg_info() const {
      std::cout
        << "\nEstimated minimum space needed for factors = " << INFO[2]
        << "\nEstimated space needed for factor          = " << INFO[3]
        << "\nEstimate of the rank of the matrix         = " << INFO[4]
        << "\nEntries dropped from the data structure    = " << INFO[5]
        << "\nOrder of the largest non-triangular block  = " << INFO[6]
        << "\nSum of orders of non-triangular blocks     = " << INFO[7]
        << "\n# of entries in the non-triangular blocks  = " << INFO[8]
        << "\nstructural rank of the matrix.             = " << INFO[9]
        << "\n# of multiple entries in the input matrix  = " << INFO[10]
        << "\n# of entries with out-of-range indices     = " << INFO[11]
        << "\n";
    }

    void
    msg_infor() const {
      std::cout
        << "\nNFLOPs for LU factorization    = " << RINFO[0]
        << "\nNFLOPs for assembly process    = " << RINFO[1]
        << "\nNFLOPs for elimination process = " << RINFO[2];
      if ( ICNTL[10] > 0 )
        std::cout
          << "\nInfinity norm of the matrix    = " << RINFO[3]
          << "\nInfinity norm of the solution  = " << RINFO[4]
          << "\nNorm of scaled residual        = " << RINFO[5];
      std::cout << "\n";
    }

    void
    msg_error() const {
      if ( INFO[0] >= 0 ) return;
      std::cerr
        << "\n***** ERROR *****"
        << "\nnRow = " << this->nRow
        << "\nnCol = " << this->nCol
        << "\nnnz  = " << this->nnz
        << "\n";
      switch ( INFO[0] ) {
      case -1:
        std::cerr << "Value of N is out of range N = " << INFO[1] << "\n";
        break;
      case -2:
        std::cerr << "Value of NE is out of range NE = " << INFO[1] << "\n";
        break;
      case -3:
        std::cerr << "JOB has wrong value or analisys was not performed prior to factorization. JOB = " << INFO[1] << "\n";
        break;
      case -4:
        std::cerr << "Error in permutation error\n";
        break;
      case -5:
        std::cerr << "Not enought space to preprocess the input matrix\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      case -6:
        std::cerr << "The matrix is structurally singular\n"
                  << "Estimated rank = " << INFO[1] << "\n";
        break;
      case -7:
        std::cerr << "Error from analysis\n"
                  << "MAXIS should be increased to at least " << INFO[1] << "\n";
        break;
      case -8:
        std::cerr << "Error from numerical factorization\n"
                  << "MAXIS should be increased to at least " << INFO[1] << "\n";
        break;
      case -9:
        std::cerr << "Error from numerical factorization\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      case -10:
        std::cerr << "The matrix is numerically singular\n"
                  << "Estimated rank = " << INFO[1] << "\n";
        break;
      case -11:
        std::cerr << "Error from the solution phase\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      case -12:
        std::cerr << "Not enought space to postprocess the solution\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      }
      std::cerr << "***** ERROR *****\n";
    }

    int
    factorize() {
      if ( this->verbose ) std::cout << "ma48::factorize() perform analisys\n";

      // analysis phase
      A_work.resize( 3*A.size() );
      std::fill( A_work.begin(), A_work.end(), 0 );
      std::copy( A.begin(), A.end(), A_work.begin() );

      IRN_work.resize( 3*A.size() );
      std::fill( IRN_work.begin(), IRN_work.end(), 0 );
      std::copy( IRN.begin(), IRN.end(), IRN_work.begin() );

      JCN_work.resize( 3*A.size() );
      std::fill( JCN_work.begin(), JCN_work.end(), 0 );
      std::copy( JCN.begin(), JCN.end(), JCN_work.begin() );

      HSL::ma48a<T>(
        this->nRow, this->nCol, this->nnz, 1,
        integer(A_work.size()), &A_work.front(),
        &IRN_work.front(), &JCN_work.front(), &KEEP.front(),
        this->CNTL, this->ICNTL, &IW.front(),
        this->INFO, this->RINFO
      );

      for (
        unsigned kkk = 0;
        this->INFO[3] > A_work.size() && kkk < 10;
        ++kkk
      ) {
        size_t newsz = size_t(this->INFO[3]);
        // Increase internal workspace:
        A_work.resize( newsz );
        std::fill( A_work.begin(), A_work.end(), 0 );
        std::copy( A.begin(), A.end(), A_work.begin() );

        IRN_work.resize( newsz );
        std::fill( IRN_work.begin(), IRN_work.end(), 0 );
        std::copy( IRN.begin(), IRN.end(), IRN_work.begin() );

        JCN_work.resize( newsz );
        std::fill( JCN_work.begin(), JCN_work.end(), 0 );
        std::copy( JCN.begin(), JCN.end(), JCN_work.begin() );

        HSL::ma48a<T>(
          this->nRow, this->nCol, this->nnz, 1,
          integer(A_work.size()), &A_work.front(),
          &IRN_work.front(), &JCN_work.front(), &KEEP.front(),
          this->CNTL, this->ICNTL, &IW.front(),
          this->INFO, this->RINFO
        );
      }

      if ( this->verbose ) msg_error();

      if ( this->INFO[0] < 0 ) return this->INFO[0];

      // factorization phase
      // -------------------
      HSL::ma48b<real_type>(
        this->nRow, this->nCol, this->nnz, 1,
        integer(A_work.size()), &A_work.front(),
        &IRN_work.front(), &JCN_work.front(), &KEEP.front(),
        this->CNTL, this->ICNTL,
        &W.front(), &IW.front(),
        this->INFO, this->RINFO
      );

      if ( this->verbose ) {
        msg_info();
        msg_error();
        std::cout << "done\n";
      }
      return this->INFO[0];
    }

  public:

    MA48() {}
    ~MA48() {}

    void
    init(void) {
      IRN.clear();
      JCN.clear();
      A.clear();
    }

    void
    insert(
      integer   i,
      integer   j,
      real_type a
    ) {
      SPARSETOOL_ASSERT(
        i >= 0 && i < this->nRow && j >= 0 && j < this->nCol,
        "MA48::insert( " << i << ", " << j << ", a ) out of matrix"
      )
      IRN.push_back(i+1);
      JCN.push_back(j+1);
      A.push_back(a);
    }

    int
    setup( bool v ) {
      this->verbose = v;
      if ( this->verbose ) std::cout << "ma48::setup()\n";
      // set pars
      HSL::ma48i<real_type>(this->CNTL, this->ICNTL);
      //this->ICNTL[6] = 0; // handle structurally rank deficient matrices.
      this->KEEP.resize(this->nRow+9*this->nCol+7);
      this->IW.resize(6*this->nRow+3*this->nCol);
      this->W.resize(this->nRow);
      if ( this->verbose ) std::cout << "ma48::setup() done\n";
      return factorize();
    }

    void
    solve( real_type X[], real_type const RHS[], bool transposed ) {
      real_type ERROR[3];
      integer   NM = std::max(this->nRow,this->nCol);
      W.resize(2*NM);
      HSL::ma48c<real_type>(
        this->nRow, this->nCol, (transposed ? 1 : 0), 1,
        integer(A_work.size()), &A_work.front(),
        &IRN_work.front(), &KEEP.front(),
        this->CNTL, this->ICNTL,
        RHS, X, ERROR,
        &W.front(), &IW.front(), this->INFO
      );

      if ( this->verbose ) {
        msg_infor();
        msg_error();
      }
    }

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat, bool v = false ) {
      this->init();
      IRN.reserve( Mat.nnz() );
      JCN.reserve( Mat.nnz() );
      A.reserve( Mat.nnz() );
      this->nRow = Mat.numRows();
      this->nCol = Mat.numCols();
      this->nnz  = Mat.nnz();
      for ( Mat.Begin(); Mat.End(); Mat.Next() )
        this->insert( Mat.row(), Mat.column(), Mat.value() );

      return this->setup( v );
    }

    int
    solve(
      Vector<T> const & b,
      Vector<T>       & x,
      bool transpose = false
    ) {
      this->solve( &x.front(), &b.front(), transpose );
      return 0;
    }

  };
}

namespace SparseToolLoad {
  using ::SparseTool::MA48;
}

#endif
