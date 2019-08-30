/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MUMPS    |
 |                                                                          |
 |  date         : 2011, 17 July                                            |
 |  version      : 1.0                                                      |
 |  file         : SparseTool_MUMPS.hh                                      |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email : enrico.bertolazzi@unitn.it                       |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_MUMPS_HH
#define SPARSETOOL_MUMPS_HH

#include "../sparse_tool.hh"

#include <smumps_c.h>
#include <dmumps_c.h>
#include <cmumps_c.h>
#include <zmumps_c.h>
#include <complex>

namespace SparseTool {

  using namespace std;
  
  extern int myid;

  void initMUMPS( int argc, char *argv[] );
  void endMUMPS();

  static int const JOB_END           = -2;
  static int const JOB_INIT          = -1;
  static int const JOB_ANALYSIS      =  1;
  static int const JOB_FACTORIZATION =  2;
  static int const JOB_SOLVE         =  3;
  static int const USE_COMM_WORLD    = -987654;

  static int const MATRIX_GE   = 0;
  static int const MATRIX_SPD  = 1;
  static int const MATRIX_SYMM = 2;

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

  inline void scalarAssign( float                & res, float                const & a ) { res = a; }
  inline void scalarAssign( double               & res, double               const & a ) { res = a; }
  inline void scalarAssign( mumps_complex        & res, mumps_complex        const & a ) { res.r = a.r; res.i = a.i; }
  inline void scalarAssign( mumps_double_complex & res, mumps_double_complex const & a ) { res.r = a.r; res.i = a.i; }
  inline void scalarAssign( mumps_complex        & res, complex<float>       const & a ) { res.r = a.real(); res.i = a.imag(); }
  inline void scalarAssign( mumps_double_complex & res, complex<double>      const & a ) { res.r = a.real(); res.i = a.imag(); }

  template <typename T>
  class MUMPS {
  public:  
    typedef T valueType;
    typedef typename MUMPSType<T>::stuctType  stuctType;
    typedef typename MUMPSType<T>::scalarType scalarType;

  private:
  
    stuctType mumpsStruct;

    indexType  numRows, nnz;
    Vector<scalarType> Avec;
    Vector<int>        Ivec;
    Vector<int>        Jvec;

    void
    check( char const message[] ) {
      if ( mumpsStruct.info[0] >= 0 ) return;
      cerr << message
           << "\ninfo[0] = " << mumpsStruct.info[0]
           << '\n';
      exit(1);
    }

    int load();

    void free() {
      mumpsStruct.job = JOB_END;
      mumps_c( mumpsStruct );
      check( "MUMPS at free" );
    }

  public:
  
    MUMPS() {
      mumpsStruct.comm_fortran = USE_COMM_WORLD; // for MPI (also required if MPI=disabled)
      mumpsStruct.par = 1;         // will be used for factorization AND solution
      mumpsStruct.sym = MATRIX_GE; // matrix is not symmetric 
      mumpsStruct.job = JOB_INIT;  // initialize solver
      mumpsStruct.irn = 0; // rows
      mumpsStruct.jcn = 0; // columns
      mumpsStruct.a   = 0; // values
      mumps_c(mumpsStruct); // call MUMPS
    };

    ~MUMPS() { free(); }

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat ) {
      // prepare
      numRows = Mat.numRows();
      nnz     = Mat.nnz();

      Ivec.resize( nnz );
      Jvec.resize( nnz );
      Avec.resize( nnz );
      indexType ii = 0;
      for ( Mat.Begin(); Mat.End(); Mat.Next(), ++ii ) {
        Ivec(ii) = Mat.row()    + 1;
        Jvec(ii) = Mat.column() + 1;
        scalarAssign(Avec(ii), Mat.value());
      }

      mumpsStruct.n   = numRows;
      mumpsStruct.nz  = nnz;
      mumpsStruct.irn = &Ivec.front(); // row indices
      mumpsStruct.jcn = &Jvec.front(); // column indices
      mumpsStruct.a   = &Avec.front(); // values
      mumpsStruct.rhs = 0; // right hand side will be set later
      // solver params
      mumpsStruct.icntl[0] = 0; // output errors
      mumpsStruct.icntl[1] = 0; // output detailed messages
      mumpsStruct.icntl[2] = 0; // output info messages
      mumpsStruct.icntl[3] = 0; // verbose level
      mumpsStruct.icntl[5] = 7; // automatically find out if to scale or permute
      // mumpsStruct.icntl[6] = 7; // automatic pivot order for the factorization
      // mumpsStruct.icntl[7] = 7; // simultaneous row and colum iterative scaling (for not symmetric matrices)
      // mumpsStruct.icntl[7] = 1; // diagonal scaling (for symmetric matrices)
  
      // mumpsStruct.icntl[0] = 6; // output errors
      // mumpsStruct.icntl[1] = 6; // output detailed messages
      // mumpsStruct.icntl[2] = 6; // output info messages
      // mumpsStruct.icntl[3] = 0; // verbose level

      // analysis phase
      mumpsStruct.job = JOB_ANALYSIS;
      mumps_c(mumpsStruct); // call MUMPS
      check( "MUMPS error during analysis and reordering" );

      // factorization
      mumpsStruct.job = JOB_FACTORIZATION;
      mumps_c(mumpsStruct); // call MUMPS
      check( "MUMPS error during factorization" );
      return 0;
    }

    int
    solve(
      Vector<T> const & b,
      Vector<T>       & x,
      bool transpose = false
    ) {

      if ( x.size() < mumpsStruct.n ) x . resize( mumpsStruct.n );
      x = b;
      // MUMPS solution
      mumpsStruct.job = JOB_SOLVE;
      mumpsStruct.rhs = (typename MUMPSType<T>::scalarType *)&x.front();
      mumps_c(mumpsStruct);
      check("MUMPS error during solution");
      mumpsStruct.rhs = 0;
      return 0;
    }
  };

  inline
  void
  initMUMPS( int argc, char *argv[] ) {
    #ifdef MUMPS_USE_MPI
    static int myid;
    int ierr;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    #endif
  }

  inline
  void
  endMUMPS() {
    #ifdef MUMPS_USE_MPI
    int ierr = MPI_Finalize();
    #endif
  }

}

namespace SparseToolLoad {
  using ::SparseTool::MUMPS;
  using ::SparseTool::initMUMPS;
  using ::SparseTool::endMUMPS;
}

#endif
