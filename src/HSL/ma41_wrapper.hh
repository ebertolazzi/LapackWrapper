/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  MA41                                                                    |
 |                                                                          |
 |  date         : April 29, 2004                                           |
 |  release      : 0.1                                                      |
 |  release_date : 2004, April 29                                           |
 |  file         : ma41.hh                                                  |
 |  author       : Enrico Bertolazzi                                        |
 |  email        : enrico.bertolazzi@ing.unitn.it                           |
 |                                                                          |
 |  This program is free software; you can redistribute it and/or modify    |
 |  it under the terms of the GNU General Public License as published by    |
 |  the Free Software Foundation; either version 2, or (at your option)     |
 |  any later version.                                                      |
 |                                                                          |
 |  This program is distributed in the hope that it will be useful,         |
 |  but WITHOUT ANY WARRANTY; without even the implied warranty of          |
 |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           |
 |  GNU General Public License for more details.                            |
 |                                                                          |
 |  You should have received a copy of the GNU General Public License       |
 |  along with this program; if not, write to the Free Software             |
 |  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               |
 |                                                                          |
 |  Copyright (C) 1999                                                      |
 |                                                                          |
 |            ___    ____  ___   _   _        ___    ____  ___   _   _      |
 |           /   \  /     /   \  \  /        /   \  /     /   \  \  /       |
 |          /____/ /__   /____/   \/        /____/ /__   /____/   \/        |
 |         /   \  /     /  \      /        /   \  /     /  \      /         |
 |        /____/ /____ /    \    /        /____/ /____ /    \    /          |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Meccanica e Strutturale                  |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |                                                                          |
\*--------------------------------------------------------------------------*/
#pragma once

#ifndef MA41_dot_HH
#define MA41_dot_HH

#include "hsl_solver.hh"
#include <vector>

namespace lapack_wrapper {
  
  template < typename real>
  class MA41 {
  public:
    // the matrix
    typedef std::vector<real>    vector_value;
    typedef std::vector<integer> vector_integer;

  private:
    enum {
      ANALISYS                             = 1,
      FACTORIZATION                        = 2,
      SOLVE                                = 3,
      ANALISYS_and_FACTORIZATION           = 4,
      FACTORIZATION_and_SOLVE              = 5,
      ANALISYS_and_FACTORIZATION_and_SOLVE = 6
    };

    bool    verbose;
    integer N, NE;
    integer KEEP[50], ICNTL[20], INFO[20];
    real    CNTL[10], RINFO[20];

    // work arrays
    vector_integer IS;
    vector_value   COLSCA, ROWSCA, S;

    // the matrix
    vector_value   A;
    vector_integer IRN, JCN;
  
    void msg_info( ostream_type & s ) const;
    void msg_infor( ostream_type & s ) const;
    void msg_error( ostream_type & s ) const;

    // private functionalities
    void load_matrix( integer nr );
    void symbfac();

  private:
    MA41( MA41 const & ) = delete;

  public:
    MA41() : verbose(false) {}

    void
    init(void) {
      IRN.clear() ;
      JCN.clear() ;
      A.clear() ;
    }

    void
    insert( integer i, integer j, real a ) {
      IRN.push_back(i+1) ;
      JCN.push_back(j+1) ;
      A.push_back(a) ;
    }

    void
    setup( integer nr, bool v = false) {
      verbose = v ;
      load_matrix(nr);
      symbfac();
    }

    void solve( real RHS[] );
  };
}

#endif
