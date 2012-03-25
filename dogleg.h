// Copyright 2011 Oblong Industries
// License: GNU LGPL <http://www.gnu.org/licenses>.

#pragma once

#include <cholmod.h>

typedef void (dogleg_callback_t)(const double*   p,
                                 double*         x,
                                 cholmod_sparse* Jt,
                                 void*           cookie);

// an operating point of the solver
typedef struct
{
  double*         p;
  double*         x;
  double          norm2_x;
  cholmod_sparse* Jt;
  double*         Jt_x;

  // the cached update vectors. It's useful to cache these so that when a step is rejected, we can
  // reuse these when we retry
  double*        updateCauchy;
  cholmod_dense* updateGN_cholmoddense;
  double         updateCauchy_lensq, updateGN_lensq; // update vector lengths

  // whether the current update vectors are correct or not
  int updateCauchy_valid, updateGN_valid;

  int didStepToEdgeOfTrustRegion;
} dogleg_operatingPoint_t;

// solver context. This has all the internal state of the solver
typedef struct
{
  cholmod_common  common;

  dogleg_callback_t* f;
  void*              cookie;

  // between steps, beforeStep contains the operating point of the last step.
  // afterStep is ONLY used while making the step. Externally, use beforeStep
  // unless you really know what you're doing
  dogleg_operatingPoint_t* beforeStep;
  dogleg_operatingPoint_t* afterStep;

  // The result of the last JtJ factorization performed. Note that JtJ is not
  // necessarily factorized at every step, so this is NOT guaranteed to contain
  // the factorization of the most recent JtJ
  cholmod_factor*          factorization;

  // Have I ever seen a singular JtJ? If so, I add this constant to the diagonal
  // from that point on. This is a simple and fast way to deal with
  // singularities. This constant starts at 0, and is increased every time a
  // singular JtJ is encountered. This is suboptimal but works for me for now.
  double                   lambda;
} dogleg_solverContext_t;


// The main optimization callback. Initial estimate of the solution passed in p,
// final optimized solution returned in p. p has Nstate variables. There are
// Nmeas measurements, the jacobian matrix has NJnnz non-zero entries.
//
// The evaluation function is given in the callback f; this function is passed
// the given cookie
//
// If we want to get the full solver state when we're done, a non-NULL
// returnContext pointer can be given. If this is done then THE USER IS
// RESPONSIBLE FOR FREEING ITS MEMORY WITH dogleg_freeContext()
double dogleg_optimize(double* p, unsigned int Nstate,
                       unsigned int Nmeas, unsigned int NJnnz,
                       dogleg_callback_t* f, void* cookie,
                       dogleg_solverContext_t** returnContext);

void dogleg_testGradient(unsigned int var, const double* p0,
                         unsigned int Nstate, unsigned int Nmeas, unsigned int NJnnz,
                         dogleg_callback_t* f, void* cookie);


// If we want to get the full solver state when we're done optimizing, we can
// pass a non-NULL returnContext pointer to dogleg_optimize(). If we do this,
// then the user MUST call dogleg_freeContext() to deallocate the pointer when
// the USER is done
void dogleg_freeContext(dogleg_solverContext_t** ctx);


////////////////////////////////////////////////////////////////
// solver parameters
////////////////////////////////////////////////////////////////
void dogleg_setMaxIterations(int n);
void dogleg_setTrustregionUpdateParameters(double downFactor, double downThreshold,
                                           double upFactor,   double upThreshold);

// lots of solver-related debug output when on
void dogleg_setDebug(int debug);


// The following parameters reflect the scaling of the problem being solved, so
// the user is strongly encouraged to tweak these. The defaults are set
// semi-arbitrarily

// The size of the trust region at start. It is cheap to reject a too-large
// trust region, so this should be something "large". Say 10x the length of an
// "expected" step size
void dogleg_setInitialTrustregion(double t);

// termination thresholds. These really depend on the scaling of the input
// problem, so the user should set these appropriately
//
// Jt_x threshold:
// The function being minimized is E = norm2(x) where x is a function of the state p.
// dE/dp = 2*Jt*x where Jt is transpose(dx/dp).
//   if( for every i  fabs(Jt_x[i]) < JT_X_THRESHOLD )
//   { we are done; }
//
// update threshold:
//   if(for every i  fabs(state_update[i]) < UPDATE_THRESHOLD)
//   { we are done; }
//
// trust region threshold:
//   if(trustregion < TRUSTREGION_THRESHOLD)
//   { we are done; }
//
// to leave a particular threshold alone, use a value <= 0 here
void dogleg_setThresholds(double Jt_x, double update, double trustregion);

