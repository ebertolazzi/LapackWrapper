/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MUMPS    |
 |                                                                          |
 |  date         : 2011, 17 July                                            |
 |  version      : 1.0                                                      |
 |  file         : MUMPS.hh                                                 |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email : enrico.bertolazzi@unitn.it                       |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once
#ifndef SPARSETOOL_MUMPS_dot_HH
#define SPARSETOOL_MUMPS_dot_HH

#include "../sparse_tool.hh"
#include "../../lapack_wrapper_config.hh"

#include <smumps_c.h>
#include <dmumps_c.h>
#include <cmumps_c.h>
#include <zmumps_c.h>
#include <complex>

namespace Sparse_tool {

  extern int myid;

  void init_MUMPS( int argc, char *argv[] );
  void end_MUMPS();

  static int const JOB_END           { -2 };
  static int const JOB_INIT          { -1 };
  static int const JOB_ANALYSIS      {  1 };
  static int const JOB_FACTORIZATION {  2 };
  static int const JOB_SOLVE         {  3 };
  static int const USE_COMM_WORLD    { -987654 };

  static int const MATRIX_GE   { 0 };
  static int const MATRIX_SPD  { 1 };
  static int const MATRIX_SYMM { 2 };

  template <typename T> struct MUMPSType {};

  template <> struct MUMPSType<float>                { typedef SMUMPS_STRUC_C stuctType; typedef float                scalarType; };
  template <> struct MUMPSType<double>               { typedef DMUMPS_STRUC_C stuctType; typedef double               scalarType; };
  template <> struct MUMPSType<mumps_complex>        { typedef CMUMPS_STRUC_C stuctType; typedef mumps_complex        scalarType; };
  template <> struct MUMPSType<mumps_double_complex> { typedef ZMUMPS_STRUC_C stuctType; typedef mumps_double_complex scalarType; };
  template <> struct MUMPSType<complex<float> >      { typedef CMUMPS_STRUC_C stuctType; typedef mumps_complex        scalarType; };
  template <> struct MUMPSType<complex<double> >     { typedef ZMUMPS_STRUC_C stuctType; typedef mumps_double_complex scalarType; };

  template <typename T> void mumps_c( T & );

  template <> inline void mumps_c( MUMPSType<float>::stuctType                & m ) { smumps_c( &m ); }
  template <> inline void mumps_c( MUMPSType<double>::stuctType               & m ) { dmumps_c( &m ); }
  template <> inline void mumps_c( MUMPSType<mumps_complex>::stuctType        & m ) { cmumps_c( &m ); }
  template <> inline void mumps_c( MUMPSType<mumps_double_complex>::stuctType & m ) { zmumps_c( &m ); }

  inline void scalar_assign( float                & res, float                const & a ) { res = a; }
  inline void scalar_assign( double               & res, double               const & a ) { res = a; }
  inline void scalar_assign( mumps_complex        & res, mumps_complex        const & a ) { res.r = a.r; res.i = a.i; }
  inline void scalar_assign( mumps_double_complex & res, mumps_double_complex const & a ) { res.r = a.r; res.i = a.i; }
  inline void scalar_assign( mumps_complex        & res, complex<float>       const & a ) { res.r = a.real(); res.i = a.imag(); }
  inline void scalar_assign( mumps_double_complex & res, complex<double>      const & a ) { res.r = a.real(); res.i = a.imag(); }

  template <typename T>
  class MUMPS {
  public:
    using real_type  = T;
    using stuctType  = typename MUMPSType<T>::stuctType;
    using scalarType = typename MUMPSType<T>::scalarType;

  private:

    stuctType          m_mumps_struct;
    integer            m_nrows;
    integer            m_nnz;
    Vector<scalarType> m_A_vec;
    Vector<int>        m_I_vec;
    Vector<int>        m_J_vec;

    void
    check( string_view message ) {
      if ( m_mumps_struct.info[0] >= 0 ) return;
      fmt::print( cerr,
        "{}\n"
        "info[0] = {}\n",
        message, m_mumps_struct.info[0]
      );
      exit(1);
    }

    int load();

    void
    free() {
      m_mumps_struct.job = JOB_END;
      mumps_c( m_mumps_struct );
      check( "MUMPS at free" );
    }

  public:

    MUMPS() {
      m_mumps_struct.comm_fortran = USE_COMM_WORLD; // for MPI (also required if MPI=disabled)
      m_mumps_struct.par          = 1;         // will be used for factorization AND solution
      m_mumps_struct.sym          = MATRIX_GE; // matrix is not symmetric
      m_mumps_struct.job          = JOB_INIT;  // initialize solver
      m_mumps_struct.irn          = 0; // rows
      m_mumps_struct.jcn          = 0; // columns
      m_mumps_struct.a            = 0; // values
      mumps_c(m_mumps_struct); // call MUMPS
    };

    ~MUMPS() { free(); }

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat ) {
      // prepare
      m_nows = Mat.nrows();
      m_nnz  = Mat.nnz();

      m_I_vec.resize( m_nnz );
      m_J_vec.resize( m_nnz );
      m_A_vec.resize( m_nnz );
      integer ii{0};
      for ( Mat.Begin(); Mat.End(); Mat.Next(), ++ii ) {
        m_I_vec(ii) = Mat.row()    + 1;
        m_J_vec(ii) = Mat.column() + 1;
        scalar_assign(m_A_vec(ii), Mat.value());
      }

      m_mumps_struct.n   = m_nrows;
      m_mumps_struct.nz  = m_nnz;
      m_mumps_struct.irn = &m_I_vec.front(); // row indices
      m_mumps_struct.jcn = &m_J_vec.front(); // column indices
      m_mumps_struct.a   = &m_A_vec.front(); // values
      m_mumps_struct.rhs = 0; // right hand side will be set later
      // solver params
      m_mumps_struct.icntl[0] = 0; // output errors
      m_mumps_struct.icntl[1] = 0; // output detailed messages
      m_mumps_struct.icntl[2] = 0; // output info messages
      m_mumps_struct.icntl[3] = 0; // verbose level
      m_mumps_struct.icntl[5] = 7; // automatically find out if to scale or permute
      // m_mumps_struct.icntl[6] = 7; // automatic pivot order for the factorization
      // m_mumps_struct.icntl[7] = 7; // simultaneous row and colum iterative scaling (for not symmetric matrices)
      // m_mumps_struct.icntl[7] = 1; // diagonal scaling (for symmetric matrices)

      // m_mumps_struct.icntl[0] = 6; // output errors
      // m_mumps_struct.icntl[1] = 6; // output detailed messages
      // m_mumps_struct.icntl[2] = 6; // output info messages
      // m_mumps_struct.icntl[3] = 0; // verbose level

      // analysis phase
      m_mumps_struct.job = JOB_ANALYSIS;
      mumps_c(m_mumps_struct); // call MUMPS
      check( "MUMPS error during analysis and reordering" );

      // factorization
      m_mumps_struct.job = JOB_FACTORIZATION;
      mumps_c(m_mumps_struct); // call MUMPS
      check( "MUMPS error during factorization" );
      return 0;
    }

    int
    solve(
      Vector<T> const & b,
      Vector<T>       & x,
      bool transpose = false
    ) {

      if ( x.size() < m_mumps_struct.n ) x.resize( m_mumps_struct.n );
      x = b;
      // MUMPS solution
      m_mumps_struct.job = JOB_SOLVE;
      m_mumps_struct.rhs = (typename MUMPSType<T>::scalarType *)&x.front();
      mumps_c(m_mumps_struct);
      check("MUMPS error during solution");
      m_mumps_struct.rhs = 0;
      return 0;
    }
  };

  inline
  void
  init_MUMPS( int argc, char *argv[] ) {
    #ifdef MUMPS_USE_MPI
    static int myid;
    int ierr;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    #endif
  }

  inline
  void
  end_MUMPS() {
    #ifdef MUMPS_USE_MPI
    int ierr = MPI_Finalize();
    #endif
  }

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::MUMPS;
  using ::Sparse_tool::init_MUMPS;
  using ::Sparse_tool::end_MUMPS;
}
#endif

#endif
