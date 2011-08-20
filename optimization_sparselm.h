#pragma once

#include <cholmod.h>

typedef void (optimizationFunction_splm_t)(const double*   p,
                                           double*         x,
                                           cholmod_sparse* Jt,
                                           void*           cookie);

double optimize_sparseLM(double* p, unsigned int Nstate,
                         unsigned int numMeasurements, unsigned int numNonzeroJacobianElements,
                         optimizationFunction_splm_t* f, void* cookie);

void testGradient(unsigned int var, const double* p0,
                  unsigned int Nstate, unsigned int Nmeas, unsigned int Jnnz,
                  optimizationFunction_splm_t* callback, void* cookie);
