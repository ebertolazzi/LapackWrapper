/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MA41     |
 |                                                                          |
 |  date         : 2011, 18 July                                            |
 |  version      : 1.0.                                                     |
 |  file         : MA41.hh                                                  |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Università degli Studi di Trento                         |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once
#ifndef SPARSETOOL_MA41_dot_HH
#define SPARSETOOL_MA41_dot_HH

#include "../sparse_tool.hh"
#include "../../HSL/hsl.h"
#include "../../lapack_wrapper_config.hh"

namespace Sparse_tool {

  template <typename T>
  class MA41 {
  public:
    using real_type = T;
    using integer   = HSL::integer;

  private:
    enum {
      ANALISYS                             = 1,
      FACTORIZATION                        = 2,
      SOLVE                                = 3,
      ANALISYS_and_FACTORIZATION           = 4,
      FACTORIZATION_and_SOLVE              = 5,
      ANALISYS_and_FACTORIZATION_and_SOLVE = 6
    };

    bool      verbose;
    integer   N, NE;
    integer   KEEP[50];
    integer   ICNTL[20];
    integer   INFO[20];
    real_type CNTL[10];
    real_type RINFO[20];

    // work arrays
    vector<integer>   IS;
    vector<real_type> COLSCA, ROWSCA, S;

    // the matrix
    vector<real_type> A;
    vector<integer>   IRN, JCN;

    void msg_info() const;
    void msg_infor() const;
    void msg_error() const;
    void load_matrix( integer nr );
    void symbfac();

  public:

    MA41() {}
    ~MA41() {}

    void init();
    void insert( integer i, integer j, real_type a );
    void setup( integer nr, bool v = false );

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat ) {
      this->init();
      this->IRN.reserve( Mat.nnz() );
      this->JCN.reserve( Mat.nnz() );
      this->A.reserve( Mat.nnz() );
      this->N = Mat.nrows();
      for ( Mat.Begin(); Mat.End(); Mat.Next() )
        this->insert( Mat.row(), Mat.column(), Mat.value() );
      this->setup( Mat.nrows() );
      return 0;
    }

    int solve( real_type RHS[], bool transpose = false );

    int solve( T const b[], T x[], bool transpose = false );

    int
    solve(
      Vector<T> const & b,
      Vector<T>       & x,
      bool transpose = false
    ) {
      return this->solve( &b.front(), &x.front(), transpose  );
    }

  };
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::MA41;
}
#endif

#endif
