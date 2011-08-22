#pragma once

typedef void (optimizationFunction_splm_t)(const double*     p,
                                           double*           x,
                                           struct splm_crsm* J,
                                           void*             cookie);

double optimize_sparseLM(double* p, int n,
                         int numMeasurements, int numNonzeroJacobianElements,
                         optimizationFunction_splm_t* f, void* cookie);
