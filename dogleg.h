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

  dogleg_operatingPoint_t* beforeStep;
  dogleg_operatingPoint_t* afterStep;
  cholmod_factor*          factorization;

  // Have I ever seen a singular JtJ? If so, I add a small constant to the
  // diagonal from that point on. This is a simple and fast way to deal with
  // singularities. This is suboptimal but works for me for now.
  int               wasPositiveSemidefinite;
} dogleg_solverContext_t;



double dogleg_optimize(double* p, unsigned int Nstate,
                       unsigned int Nmeas, unsigned int NJnnz,
                       dogleg_callback_t* f, void* cookie);

void dogleg_testGradient(unsigned int var, const double* p0,
                         unsigned int Nstate, unsigned int Nmeas, unsigned int NJnnz,
                         dogleg_callback_t* f, void* cookie);




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

