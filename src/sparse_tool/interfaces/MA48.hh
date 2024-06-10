/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MA41     |
 |                                                                          |
 |  date         : 2011, 18 July                                            |
 |  version      : 1.0.                                                     |
 |  file         : MA48.hh                                                  |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once
#ifndef SPARSETOOL_MA48_dot_HH
#define SPARSETOOL_MA48_dot_HH

#include "../sparse_tool.hh"
#include "../../HSL/hsl.h"
#include "../../lapack_wrapper_config.hh"

namespace Sparse_tool {

  template <typename T>
  class MA48 {
  public:
    using real_type = T;
    using integer   = HSL::integer;

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

    void msg_info() const;
    void msg_infor() const;
    void msg_error() const;
    int  factorize();

  public:

    MA48() {}
    ~MA48() {}

    void init();
    void insert( integer i, integer j, real_type a );
    int  setup( bool v );
    void solve( real_type X[], real_type const RHS[], bool transposed );

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat, bool v = false ) {
      this->init();
      IRN.reserve( Mat.nnz() );
      JCN.reserve( Mat.nnz() );
      A.reserve( Mat.nnz() );
      this->nRow = Mat.nrows();
      this->nCol = Mat.ncols();
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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::MA48;
}
#endif

#endif
