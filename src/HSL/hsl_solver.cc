/**
 * \file hslsolver.cc
 * Header definitions of the class HSLsolver.
 *
 * \author A. Huber and E.Bertolazzi
 * \since  13.11.2018
 */

#include "hsl_solver.hh"

namespace lapack_wrapper {

  template <typename real>
  HSLsolver<real>::HSLsolver()
  : isInitialized(false)
  , isFactorized(false)
  , last_error("")
  { }

  template <typename real>
  HSLsolver<real>::~HSLsolver() {
  }

  template class HSLsolver<float>;
  template class HSLsolver<double>;

} // namespace lapack_wrapper

//
// EOF: hsl_solver.cc
//
