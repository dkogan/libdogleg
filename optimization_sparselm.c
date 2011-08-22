#include <splm.h>
#include <stdio.h>
#include <math.h>
#include "optimization_sparselm.h"

// this is the solver using the sparse Levenberg-Marquardt implementation from
// http://www.ics.forth.gr/~lourakis/sparseLM/
//
// This implementation is licensed under the GPL, so we cannot use it in a commercial product
// without opening all code linked to it

#define MAX_ITERATIONS 1000

typedef struct
{
  optimizationFunction_splm_t* callback;
  void*                        cookie;
} solver_cookie_t;

static void callback_cost(double *p, double *hx, int nvars, int nobs, void* _solver_cookie)
{
  // These do nothing, but prevent the compiler from warning me about these unused variables
  (void)nvars;
  (void)nobs;

  /* functional relation describing measurements. Given a parameter vector p,
   * computes a prediction of the measurements \hat{x}. p is nvarsx1,
   * \hat{x} is nobsx1, maximum
   */
  solver_cookie_t* solver_cookie = (solver_cookie_t*)_solver_cookie;
  (*solver_cookie->callback)( p, hx, NULL, solver_cookie->cookie );
}

static void callback_constraint(double *p, struct splm_crsm* jac, int nvars, int nobs, void* _solver_cookie)
{
  jac->nr = nobs;
  jac->nc = nvars;

  /* function to supply the nonzero pattern of the sparse Jacobian of func and
   * evaluate it at p in CRS format. Non-zero elements are to be stored in jac
   * which has been preallocated with a capacity of Jnnz
   */
  solver_cookie_t* solver_cookie = (solver_cookie_t*)_solver_cookie;
  (*solver_cookie->callback)( p, NULL, jac, solver_cookie->cookie );
}

double optimize_sparseLM(double* p, int n,
                         int numMeasurements, int numNonzeroJacobianElements,
                         optimizationFunction_splm_t* f, void* cookie)
{
  // input options
  // I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, delta, spsolver]. Respectively the
  // scale factor for initial \mu, stopping thresholds for ||J^T e||_inf, ||dp||_2 and
  // ||e||_2. Fifth element (i.e., delta) is unused, last element specifies the sparse direct solver
  // to employ. Set to NULL for defaults to be used.
  double opts[SPLM_OPTS_SZ] __attribute__((unused)); // I'm using the defaults

  // status report when done
  // info[0]=||e||_2 at initial p.
  // info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
  // info[5]= # iterations,
  // info[6]=reason for terminating: 1 - stopped by small gradient J^T e
  //                                 2 - stopped by small dp
  //                                 3 - stopped by itmax
  //                                 4 - singular matrix. Restart from current p with increased mu 
  //                                 5 - too many failed attempts to increase damping. Restart with increased mu
  //                                 6 - stopped by small ||e||_2
  //                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. User error
  // info[7]= # function evaluations
  // info[8]= # Jacobian evaluations
  // info[9]= # linear systems solved, i.e. # attempts for reducing error
  double info[SPLM_INFO_SZ];

  solver_cookie_t solver_cookie = {.callback = f,
                                   .cookie   = cookie};

  int iterations = sparselm_dercrs(&callback_cost, &callback_constraint,
                                   p, NULL, n, 0, numMeasurements,
                                   numNonzeroJacobianElements, -1,
                                   MAX_ITERATIONS, NULL, info,
                                   &solver_cookie);

  if( iterations > 0 )
  {
    fprintf(stderr, "success! took %d iterations\n", iterations);
    return sqrt(info[1] / ((double)numMeasurements / 2.0));
  }

  fprintf(stderr, "solver failed. error: %d\n", iterations);
  return -1.0;
}
