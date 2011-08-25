#pragma once

#include <cholmod.h>

typedef void (dogleg_callback_t)(const double*   p,
                                 double*         x,
                                 cholmod_sparse* Jt,
                                 void*           cookie);

double dogleg_optimize(double* p, unsigned int Nstate,
                       unsigned int numMeasurements, unsigned int numNonzeroJacobianElements,
                       dogleg_callback_t* f, void* cookie);

void dogleg_testGradient(unsigned int var, const double* p0,
                         unsigned int Nstate, unsigned int Nmeas, unsigned int Jnnz,
                         dogleg_callback_t* callback, void* cookie);




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

