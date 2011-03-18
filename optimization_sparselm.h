#pragma once

typedef void (optimizationFunction_splm_t)(const double*     p,
                                           double*           x,
                                           struct splm_crsm* J,
                                           void*             cookie);

double optimize_sparseLM(double* p, int n,
                         int numMeasurements, int numNonzeroJacobianElements,
                         optimizationFunction_splm_t* f, void* cookie);

void testGradient(int var, const double* p0, void* cookie, int Nstate, int Nmeas, int Jnnz, optimizationFunction_splm_t* callback);
