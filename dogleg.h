#pragma once

#include <cholmod.h>

typedef void (dogleg_callback_t)(const double*   p,
                                 double*         x,
                                 cholmod_sparse* Jt,
                                 void*           cookie);

double optimize(double* p, unsigned int Nstate,
                unsigned int numMeasurements, unsigned int numNonzeroJacobianElements,
                dogleg_callback_t* f, void* cookie);

void testGradient(unsigned int var, const double* p0,
                  unsigned int Nstate, unsigned int Nmeas, unsigned int Jnnz,
                  dogleg_callback_t* callback, void* cookie);
