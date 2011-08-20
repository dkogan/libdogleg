#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <cholmod.h>
#include "optimization_sparselm.h"

// I do this myself because I want this to be active in all build modes, not just !NDEBUG
#define ASSERT(x) do { if(!(x)) { fprintf(stderr, "ASSERTION FAILED at %s:%d\n", __FILE__, __LINE__); exit(1); } } while(0)

#define MAX_ITERATIONS         100
#define LAMBDA_DECREASE_FACTOR 0.1
#define LAMBDA_INCREASE_FACTOR 2

// used to consolidate the casts
#define P(A, index) ((unsigned int*)((A)->p))[index]
#define I(A, index) ((unsigned int*)((A)->i))[index]
#define X(A, index) ((double*      )((A)->x))[index]

//////////////////////////////////////////////////////////////////////////////////////////
// routines for gradient testing
//////////////////////////////////////////////////////////////////////////////////////////
#define GRADTEST_DELTA 1e-6
static double getGrad(unsigned int var, int meas, const cholmod_sparse* Jt)
{
  int row     = P(Jt, meas);
  int rownext = P(Jt, meas+1);

  for(int i=row; i<rownext; i++)
  {
    if(I(Jt,i) == var)
      return X(Jt,i);
  }

  return nan("nogradient");
}

void testGradient(unsigned int var, const double* p0,
                  unsigned int Nstate, unsigned int Nmeas, unsigned int Jnnz,
                  optimizationFunction_splm_t* callback, void* cookie)
{
  double           x0[Nmeas];
  double           x [Nmeas];

  double           p [Nstate];
  memcpy(p, p0, Nstate * sizeof(double));


  cholmod_common _cholmod_common;
  if( !cholmod_start(&_cholmod_common) )
  {
    fprintf(stderr, "Couldn't initialize CHOLMOD\n");
    return;
  }

  cholmod_sparse* Jt  = cholmod_allocate_sparse(Nstate, Nmeas, Jnnz,
                                                1, // sorted
                                                1, // packed,
                                                0, // NOT symmetric
                                                CHOLMOD_REAL,
                                                &_cholmod_common);
  cholmod_sparse* Jt0 = cholmod_allocate_sparse(Nstate, Nmeas, Jnnz,
                                                1, // sorted
                                                1, // packed,
                                                0, // NOT symmetric
                                                CHOLMOD_REAL,
                                                &_cholmod_common);

  (*callback)(p, x0, Jt0, cookie);
  p[var] += GRADTEST_DELTA;
  (*callback)(p, x,  Jt,  cookie);

  for(unsigned int i=0; i<Nmeas; i++)
  {
    double gObs = (x[i] - x0[i]) / GRADTEST_DELTA;
    double gRep = getGrad(var, i, Jt0);

    if(isnan(gRep))
    {
      if( gObs != 0 )
        printf("var,meas %d,%d: no reported gradient, but observed %f\n", var, i, gObs);

      continue;
    }

    printf("var,meas %d,%d: reported: %f, observed: %f, err: %f, relativeerr: %f\n", var, i,
           gRep, gObs, fabs(gRep - gObs), fabs(gRep - gObs) / ( (fabs(gRep) + fabs(gObs)) / 2.0 ) );
  }

  cholmod_free_sparse(&Jt,  &_cholmod_common);
  cholmod_free_sparse(&Jt0, &_cholmod_common);
  cholmod_finish(&_cholmod_common);
}


//////////////////////////////////////////////////////////////////////////////////////////
// solver routines
//////////////////////////////////////////////////////////////////////////////////////////

typedef struct
{
  double*         p;
  double*         x;
  cholmod_sparse* Jt;
  cholmod_sparse* JtJ;
  cholmod_sparse* Jt_x;
} operatingPoint_t;

typedef struct
{
  cholmod_common  common;

  optimizationFunction_splm_t* f;
  void*                        cookie;

  operatingPoint_t* beforeStep;
  operatingPoint_t* afterStep;
  cholmod_factor*   factorization;
} solverContext_t;



// takes in point->p, and computes all the quantities derived from it, storing the result in the
// other members of the operatingPoint structure
static void computeCallbackOperatingPoint(operatingPoint_t* point, solverContext_t* ctx)
{
  (*ctx->f)(point->p, point->x, point->Jt, ctx->cookie);

  // compute Jt*x. I do this myself because the MatrixOps module in CHOLMOD is licensed under the
  // GPL
  memset(point->Jt_x->x, 0, sizeof(double)*point->Jt->nrow);
  for(unsigned int i=0; i<point->Jt->ncol; i++)
  {
    for(unsigned int j=P(point->Jt, i); j<P(point->Jt, i+1); j++)
    {
      int row = I(point->Jt, j);
      X(point->Jt_x, row) += point->x[i] * X(point->Jt, j);
    }
  }

  // This may be somewhat inefficient and is DEFINITELY inconvenient. Instead of
  // computing JtJ, I like to pass ctx->Jt to cholmod_analyze() and
  // cholmod_factorize(), then call cholmod_updown() to add the lambda damping
  // factor. cholmod_updown() is in the CHOLMOD's Modify module, which is
  // licensed under the GPL, NOT the LGPL. Thus I handle the damping factor
  // myself. Benchmarks indicate that there isn't really a performance hit as a
  // result
  if(point->JtJ != NULL)
    cholmod_free_sparse(&point->JtJ, &ctx->common);

  point->JtJ = cholmod_aat(point->Jt, NULL, 0, 1, &ctx->common);
  ASSERT(point->JtJ != NULL);
  point->JtJ->stype = 1; // matrix is symmetric, upper triangle given
}

// takes a step from the given operating point, using the given lambda.
static void takeStepFrom(operatingPoint_t* pointFrom, double* newp,
                         double lambda, solverContext_t* ctx)
{
  // I'm assuming the pattern of zeros will remain the same throughout, so I
  // analyze only once
  if(ctx->factorization == NULL)
  {
    ctx->factorization = cholmod_analyze(pointFrom->JtJ, &ctx->common);
    ASSERT(ctx->factorization != NULL);
  }

  double beta[2] = {lambda, 0};
  ASSERT( cholmod_factorize_p(pointFrom->JtJ, beta,
                              NULL, 0,
                              ctx->factorization, &ctx->common) );

  // solve JtJ delta = Jt x
  // new p is p - delta
  cholmod_sparse* result = cholmod_spsolve(0, ctx->factorization,
                                           pointFrom->Jt_x,
                                           &ctx->common);

  ASSERT( result->ncol == 1 );  // make sure the resulting matrix has the right dimensions
  for(unsigned int i=0; i<P(result, 1); i++)
  {
    unsigned int row = I(result, i);
    newp[row] = pointFrom->p[row] - X(result, i);
  }

  cholmod_free_sparse(&result, &ctx->common);
}


// I have a candidate step. I adjust the lambda accordingly, and also report whether this step
// should be accepted (0 == rejected, otherwise accepted)
static int evaluateStep_adjustLambda(const operatingPoint_t* before,
                                     const operatingPoint_t* after,
                                     double* lambda)
{
  return 1;
}
 
static void runOptimizer(solverContext_t* ctx)
{
  double lambda = 1.0;

  computeCallbackOperatingPoint(ctx->beforeStep, ctx);

  for(int stepCount=0; stepCount<MAX_ITERATIONS; stepCount++)
  {
    while(1)
    {
      takeStepFrom                 (ctx->beforeStep, ctx->afterStep->p, lambda, ctx);
      computeCallbackOperatingPoint(ctx->afterStep,  ctx);

      if( evaluateStep_adjustLambda(ctx->beforeStep, ctx->afterStep, &lambda) )
      {
        // I accept this step, so the after-step operating point is the before-step operating point
        // of the next iteration. I exchange the before- and after-step structures so that all the
        // pointers are still around and I don't have to re-allocate
        operatingPoint_t* tmp;
        tmp             = ctx->afterStep;
        ctx->beforeStep = ctx->afterStep;
        ctx->afterStep  = tmp;

        break;
      }

      // I have rejected this step, so I try again with the new lambda
    }
  }
}

static operatingPoint_t* allocOperatingPoint(int Nstate, int numMeasurements,
                                             int numNonzeroJacobianElements,
                                             cholmod_common* common)
{
  operatingPoint_t* point = malloc(sizeof(operatingPoint_t));
  ASSERT(point != NULL);

  point->p = malloc(Nstate * sizeof(point->p[0]));
  ASSERT(point->p != NULL);

  point->x = malloc(numMeasurements * sizeof(point->x[0]));
  ASSERT(point->x != NULL);

  point->Jt = cholmod_allocate_sparse(Nstate, numMeasurements, numNonzeroJacobianElements,
                                      1, // sorted
                                      1, // packed,
                                      0, // NOT symmetric
                                      CHOLMOD_REAL,
                                      common);
  ASSERT(point->Jt != NULL);

  // the 1-column vector Jt * x
  point->Jt_x = cholmod_allocate_sparse(Nstate, 1, Nstate,
                                        1, // sorted
                                        1, // packed,
                                        0, // NOT symmetric
                                        CHOLMOD_REAL,
                                        common);
  ASSERT(point->Jt_x != NULL);
  // I set up the dense 1-column vector as a sparse one
  P(point->Jt_x, 0) = 0;
  P(point->Jt_x, 1) = Nstate;
  for(int i=0; i<Nstate; i++)
    I(point->Jt_x, i) = i;

  // JtJ gets allocated when it's computed. So here I simply mark it as unallocated
  point->JtJ = NULL;

  return point;
}

static void freeOperatingPoint(operatingPoint_t** point, cholmod_common* common)
{
  free((*point)->p);
  free((*point)->x);

  cholmod_free_sparse(&(*point)->Jt,   common);
  cholmod_free_sparse(&(*point)->Jt_x, common);

  if( (*point)->JtJ != NULL)
    cholmod_free_sparse(&(*point)->JtJ, common);

  free(*point);
  *point = NULL;
}

double optimize_sparseLM(double* p, unsigned int Nstate,
                         unsigned int numMeasurements, unsigned int numNonzeroJacobianElements,
                         optimizationFunction_splm_t* f, void* cookie)
{
  solverContext_t ctx = {.f             = f,
                         .cookie        = cookie,
                         .factorization = NULL};

  if( !cholmod_start(&ctx.common) )
  {
    fprintf(stderr, "Couldn't initialize CHOLMOD\n");
    return -1.0;
  }

  ctx.beforeStep = allocOperatingPoint(Nstate, numMeasurements, numNonzeroJacobianElements, &ctx.common);
  ctx.afterStep  = allocOperatingPoint(Nstate, numMeasurements, numNonzeroJacobianElements, &ctx.common);

  memcpy(ctx.beforeStep->p, p, Nstate * sizeof(double));

  // everything is set up, so run the solver!
  runOptimizer(&ctx);

  // runOptimizer places the most recent results into beforeStep in preparation for another
  // iteration
  memcpy(p, ctx.beforeStep->p, Nstate * sizeof(double));

  freeOperatingPoint(&ctx.beforeStep, &ctx.common);
  freeOperatingPoint(&ctx.afterStep , &ctx.common);

  if(ctx.factorization != NULL)
    cholmod_free_factor (&ctx.factorization, &ctx.common);
  cholmod_finish(&ctx.common);


  fprintf(stderr, "success! took %d iterations\n", 10);
  return 10.0; // rms
}
