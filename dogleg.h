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


// these parameters likely should be messed with
void dogleg_setDebug(int debug);
void dogleg_setInitialTrustregion(double t);
void dogleg_setThresholds(double Jt_x, double update, double trustregion);

// these parameters likely should not be messed with.
void dogleg_setMaxIterations(int n);
void dogleg_setTrustregionUpdateParameters(double downFactor, double downThreshold,
                                           double upFactor,   double upThreshold);
