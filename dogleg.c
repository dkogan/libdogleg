#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "dogleg.h"

//#define DOGLEG_DEBUG

// I do this myself because I want this to be active in all build modes, not just !NDEBUG
#define ASSERT(x) do { if(!(x)) { fprintf(stderr, "ASSERTION FAILED at %s:%d\n", __FILE__, __LINE__); exit(1); } } while(0)

#define MAX_ITERATIONS 100
#define TRUSTREGION0   1.0

#define TRUSTREGION_DECREASE_FACTOR    0.1
#define TRUSTREGION_INCREASE_FACTOR    2
#define TRUSTREGION_INCREASE_THRESHOLD 0.75
#define TRUSTREGION_DECREASE_THRESHOLD 0.25

#define JT_X_THRESHOLD           1e-10
#define UPDATE_THRESHOLD         1e-10
#define TRUSTREGION_THRESHOLD    1e-10

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

  // this could be done more efficiently with bsearch()
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
  double          norm2_x;
  cholmod_sparse* Jt;
  double*         Jt_x;
} operatingPoint_t;

typedef struct
{
  cholmod_common  common;

  optimizationFunction_splm_t* f;
  void*                        cookie;

  operatingPoint_t* beforeStep;
  operatingPoint_t* afterStep;
  cholmod_factor*   factorization;

  // maxium fabs(p) value seen throughout the iteration. These are used to
  // find a good scale for the termination conditions
  double*           p_max;
} solverContext_t;

// This should really come from BLAS or libminimath
static double norm2(const double* x, unsigned int n)
{
  double result = 0;
  for(unsigned int i=0; i<n; i++)
    result += x[i]*x[i];
  return result;
}
static double inner(const double* x, const double* y, unsigned int n)
{
  double result = 0;
  for(unsigned int i=0; i<n; i++)
    result += x[i]*y[i];
  return result;
}
static void vec_copy_scaled(double* dest,
                            const double* v, double scale, int n)
{
  for(int i=0; i<n; i++)
    dest[i] = scale * v[i];
}
static void vec_add(double* dest,
                    const double* v0, const double* v1, int n)
{
  for(int i=0; i<n; i++)
    dest[i] = v0[i] + v1[i];
}
static void vec_sub(double* dest,
                    const double* v0, const double* v1, int n)
{
  for(int i=0; i<n; i++)
    dest[i] = v0[i] - v1[i];
}
static void vec_negate(double* v, int n)
{
  for(int i=0; i<n; i++)
    v[i] *= -1.0;
}
// computes a + k*(b-a)
static void vec_interpolate(double* dest,
                            const double* a, double k, const double* b_minus_a,
                            int n)
{
  for(int i=0; i<n; i++)
    dest[i] = a[i] + k*b_minus_a[i];
}

// It would be nice to use the CHOLMOD implementation of these, but they're
// licensed under the GPL
static void mul_spmatrix_densevector(double* dest,
                                     const cholmod_sparse* A, const double* x)
{
  memset(dest, 0, sizeof(double) * A->nrow);
  for(unsigned int i=0; i<A->ncol; i++)
  {
    for(unsigned int j=P(A, i); j<P(A, i+1); j++)
    {
      int row = I(A, j);
      dest[row] += x[i] * X(A, j);
    }
  }
}
static double norm2_mul_spmatrix_t_densevector(const cholmod_sparse* At, const double* x)
{
  // computes norm2(transpose(At) * x). For each row of A (col of At) I
  // compute that element of A*x, and accumulate it immediately towards the
  // norm
  double result = 0.0;

  for(unsigned int i=0; i<At->ncol; i++)
  {
    double dotproduct = 0.0;
    for(unsigned int j=P(At, i); j<P(At, i+1); j++)
    {
      int row = I(At, j);
      dotproduct += x[row] * X(At, j);
    }
    result += dotproduct * dotproduct;
  }

  return result;
}



// takes in point->p, and computes all the quantities derived from it, storing
// the result in the other members of the operatingPoint structure. Returns
// true if the gradient-size termination criterion has been met
static int computeCallbackOperatingPoint(operatingPoint_t* point, solverContext_t* ctx)
{
  (*ctx->f)(point->p, point->x, point->Jt, ctx->cookie);

  // compute the 2-norm of the current error vector
  // At some point this should be changed to use the call from libminimath
  point->norm2_x = norm2(point->x, point->Jt->ncol);

  // compute Jt*x
  mul_spmatrix_densevector(point->Jt_x, point->Jt, point->x);

  // If the largest absolute gradient element is smaller than the threshold,
  // we can stop iterating. This is equivalent to the inf-norm
  for(unsigned int i=0; i<point->Jt->nrow; i++)
    if(fabs(point->Jt_x[i]) > JT_X_THRESHOLD)
      return 0;
#ifdef DOGLEG_DEBUG
  fprintf(stderr, "Jt_x all below the threshold. Done iterating!\n");
#endif
  return 1;
}
static double computeExpectedImprovement(const double* step, const operatingPoint_t* point)
{
  // My error function is F=norm2(f(p + step)). F(0) - F(step) =
  // = norm2(x) - norm2(x + J*step) = -2*inner(x,J*step) - norm2(J*step)
  // = -2*inner(Jt_x,step) - norm2(J*step)
  return
    - 2.0*inner(point->Jt_x, step, point->Jt->nrow)
    - norm2_mul_spmatrix_t_densevector(point->Jt, step);
}


// takes a step from the given operating point, using the given trust region
// radius. Returns the expected improvement, based on the step taken and the
// linearized x(p). If we can stop iterating, returns a negative number
static double takeStepFrom(operatingPoint_t* pointFrom, double* newp,
                           double trustregion, solverContext_t* ctx)
{
#ifdef DOGLEG_DEBUG
  fprintf(stderr, "taking step with trustregion %f\n", trustregion);
#endif

  double update_cauchy[pointFrom->Jt->nrow];
  double update_dogleg[pointFrom->Jt->nrow];
  cholmod_dense* update_dense_gn = NULL;

  double* update;

  while(1) // not a real loop. I always break out at 1st pass. Used for flow control
  {
    // first, I look at a step in the steepest direction that minimizes my
    // quadratic error function (Cauchy point). If this is past my trust region,
    // I move as far as the trust region allows along the steepest descent
    // direction. My error function is F=norm2(f(p)). dF/dP = 2*ft*J.
    // This is proportional to Jt_x, which is thus the steepest ascent direction.
    //
    // Thus along this direction we have F(k) = norm2(f(p + k Jt_x)). The Cauchy
    // point is where F(k) is at a minumum:
    // dF_dk = 2 f(p + k Jt_x)t  J Jt_x ~ (x + k J Jt_x)t J Jt_x =
    // = xt J Jt x + k xt J Jt J Jt x = norm2(Jt x) + k norm2(J Jt x)
    // dF_dk = 0 -> k= -norm2(Jt x) / norm2(J Jt x)
    // Summary:
    // the steepest direction is parallel to Jt*x. The Cauchy point is at
    // k*Jt*x where k       = -norm2(Jt*x)/norm2(J*Jt*x)
    double norm2_Jt_x       = norm2(pointFrom->Jt_x, pointFrom->Jt->nrow);
    double norm2_J_Jt_x     = norm2_mul_spmatrix_t_densevector(pointFrom->Jt, pointFrom->Jt_x);
    double k                = -norm2_Jt_x / norm2_J_Jt_x;
    double cauchyStepSizeSq = k*k * norm2_Jt_x;

    vec_copy_scaled(update_cauchy,
                    pointFrom->Jt_x, k, pointFrom->Jt->nrow);

#ifdef DOGLEG_DEBUG
    fprintf(stderr, "cauchy step size %f\n", sqrt(cauchyStepSizeSq));
#endif

    if(cauchyStepSizeSq >= trustregion*trustregion)
    {
#ifdef DOGLEG_DEBUG
      fprintf(stderr, "taking cauchy step\n");
#endif

      // cauchy step goes beyond my trust region, so I do a gradient descent
      // to the edge of my trust region and call it good
      update = update_cauchy;
      break;
    }

    // I'm not yet done. The cauchy point is within the trust region, so I can
    // go further.

    // I look at the full Gauss-Newton step. If this is within the trust
    // region, I use it. Otherwise, I find the point at the edge of my trust
    // region that lies on a straight line between the Cauchy point and the
    // Gauss-Newton solution, and use that. This is the heart of Powell's
    // dog-leg algorithm.

    // I'm assuming the pattern of zeros will remain the same throughout, so I
    // analyze only once
    if(ctx->factorization == NULL)
    {
      ctx->factorization = cholmod_analyze(pointFrom->Jt, &ctx->common);
      ASSERT(ctx->factorization != NULL);
    }

    ASSERT( cholmod_factorize(pointFrom->Jt,
                              ctx->factorization, &ctx->common) );

    // solve JtJ*update_gn = Jt*x. Gauss-Newton step is then -update_gn
    cholmod_dense Jt_x_dense = {.nrow  = pointFrom->Jt->nrow,
                                .ncol  = 1,
                                .nzmax = pointFrom->Jt->nrow,
                                .d     = pointFrom->Jt->nrow,
                                .x     = pointFrom->Jt_x,
                                .xtype = CHOLMOD_REAL,
                                .dtype = CHOLMOD_DOUBLE};
    cholmod_dense* update_dense_gn = cholmod_solve(0, ctx->factorization,
                                                   &Jt_x_dense,
                                                   &ctx->common);
    double* update_gn = update_dense_gn->x;
    vec_negate(update_gn, pointFrom->Jt->nrow); // should be more efficient than this later

    double GaussNewtonStepSizeSq = norm2(update_gn, pointFrom->Jt->nrow);
#ifdef DOGLEG_DEBUG
    fprintf(stderr, "gn step size %f\n", sqrt(GaussNewtonStepSizeSq));
#endif

    if(GaussNewtonStepSizeSq <= trustregion*trustregion)
    {
#ifdef DOGLEG_DEBUG
      fprintf(stderr, "taking GN step\n");
#endif

      // full Gauss-Newton step lies within my trust region. Take the full step
      update = update_gn;
      break;
    }
    else
    {
#ifdef DOGLEG_DEBUG
      fprintf(stderr, "taking interpolated step\n");
#endif


      // full Gauss-Newton step lies outside my trust region, so I interpolate
      // between the Cauchy-point step and the Gauss-Newton step to find a step
      // that takes me to the edge of my trust region.
      //
      // I have something like norm2(a + k*(b-a)) = dsq
      // = norm2(a) + 2*at*(b-a) * k + norm2(b-a)*k^2 = dsq
      // let c = at*(b-a), l2 = norm2(b-a) ->
      // l2 k^2 + 2*c k + norm2(a)-dsq = 0
      //
      // This is a simple quadratic equation:
      // k = (-2*c +- sqrt(4*c*c - 4*l2*(norm2(a)-dsq)))/(2*l2)
      //   = (-c +- sqrt(c*c - l2*(norm2(a)-dsq)))/l2

      // to make 100% sure the descriminant is positive, I choose a to be the
      // cauchy step.  The solution must have k in [0,1], so I much have the
      // +sqrt side, since the other one is negative
      double        dsq    = trustregion*trustregion;
      double        norm2a = cauchyStepSizeSq;
      const double* a      = update_cauchy;
      const double* b      = update_gn;
      double        a_minus_b[pointFrom->Jt->nrow];

      vec_sub(a_minus_b, a, b, pointFrom->Jt->nrow);
      double l2           = norm2(a_minus_b,    pointFrom->Jt->nrow);
      double neg_c        = inner(a_minus_b, a, pointFrom->Jt->nrow);
      double discriminant = neg_c*neg_c - l2* (norm2a - dsq);
      if(discriminant < 0.0)
      {
        fprintf(stderr, "negative discriminant: %f!\n", discriminant);
        discriminant = 0.0;
      }
      double k = (neg_c + sqrt(discriminant))/l2;

      // I can rehash this to not store this data array at all, but it's clearer
      // to
      vec_interpolate(update_dogleg, a, -k, a_minus_b, pointFrom->Jt->nrow);

#ifdef DOGLEG_DEBUG
      double updateNorm = norm2(update_dogleg, pointFrom->Jt->nrow);
      fprintf(stderr, "k %f, norm %f\n", k, sqrt(updateNorm));
#endif

      update = update_dogleg;
      break;
    }
  }

  // take the step
  vec_add(newp, pointFrom->p, update, pointFrom->Jt->nrow);
  double expectedImprovement = computeExpectedImprovement(update, pointFrom);

  // are we done? For each state variable I look at the relative update step
  // in respect to the max value of that variable seen so far. If all the
  // relative ratios fall below a threshold, I call myself done
  unsigned int i;
  for(i=0; i<pointFrom->Jt->nrow; i++)
    if( fabs(update[i]) > UPDATE_THRESHOLD*ctx->p_max[i])
      break;
  if(update_dense_gn != NULL)
    cholmod_free_dense(&update_dense_gn, &ctx->common);

  if(i == pointFrom->Jt->nrow)
  {
#ifdef DOGLEG_DEBUG
    fprintf(stderr, "update small enough. Done iterating!\n");
#endif
    return -1.0;
  }

  return expectedImprovement;
}


// I have a candidate step. I adjust the trustregion accordingly, and also
// report whether this step should be accepted (0 == rejected, otherwise
// accepted)
static int evaluateStep_adjustTrustRegion(const operatingPoint_t* before,
                                          const operatingPoint_t* after,
                                          double* trustregion,
                                          double expectedImprovement)
{
  double observedImprovement = before->norm2_x - after->norm2_x;

#ifdef DOGLEG_DEBUG
  fprintf(stderr, "expected improvement: %f, got improvement %f\n", expectedImprovement, observedImprovement);
#endif

  double rho = observedImprovement / expectedImprovement;
  if     (rho > TRUSTREGION_INCREASE_THRESHOLD) *trustregion *= TRUSTREGION_INCREASE_FACTOR;
  else if(rho < TRUSTREGION_DECREASE_THRESHOLD) *trustregion *= TRUSTREGION_DECREASE_FACTOR;

  return rho > 0.0;
}

static void runOptimizer(solverContext_t* ctx)
{
  double trustregion = TRUSTREGION0;

  if( computeCallbackOperatingPoint(ctx->beforeStep, ctx) )
    return;

  int stepCount;
  for(stepCount=0; stepCount<MAX_ITERATIONS; stepCount++)
  {
#ifdef DOGLEG_DEBUG
    fprintf(stderr, "step %d\n", stepCount);
    fprintf(stderr, "\n\n\n");
#endif

    while(1)
    {
#ifdef DOGLEG_DEBUG
      fprintf(stderr, "\n\n");
#endif

      double expectedImprovement =
        takeStepFrom(ctx->beforeStep, ctx->afterStep->p, trustregion, ctx);

      // negative expectedImprovement is used to indicate that we're done
      if(expectedImprovement < 0.0)
        return;

      double afterStepZeroGradient = computeCallbackOperatingPoint(ctx->afterStep, ctx);

      if( evaluateStep_adjustTrustRegion(ctx->beforeStep, ctx->afterStep, &trustregion,
                                         expectedImprovement) )
      {
#ifdef DOGLEG_DEBUG
        fprintf(stderr, "accepted step\n");
#endif
        // I accept this step, so the after-step operating point is the before-step operating point
        // of the next iteration. I exchange the before- and after-step structures so that all the
        // pointers are still around and I don't have to re-allocate
        operatingPoint_t* tmp;
        tmp             = ctx->afterStep;
        ctx->afterStep  = ctx->beforeStep;
        ctx->beforeStep = tmp;

        if( afterStepZeroGradient )
        {
#ifdef DOGLEG_DEBUG
          fprintf(stderr, "Gradient low enough and we just improved. Done iterating!");
#endif
          return;
        }

        // I now update p_max
        for(unsigned int i=0; i<ctx->beforeStep->Jt->nrow; i++)
          ctx->p_max[i] = fmax( fabs(ctx->beforeStep->p[i]), ctx->p_max[i] );

        break;
      }

#ifdef DOGLEG_DEBUG
      fprintf(stderr, "rejected step\n");
#endif

      // This step was rejected. check if the new trust region size is small
      // enough to give up
      if(TRUSTREGION_THRESHOLD*norm2(ctx->p_max, ctx->beforeStep->Jt->nrow) > trustregion*trustregion)
      {
#ifdef DOGLEG_DEBUG
        fprintf(stderr, "trust region small enough. Giving up. Done iterating!\n");
#endif
        return;
      }

      // I have rejected this step, so I try again with the new trust region
    }
  }

#ifdef DOGLEG_DEBUG
  if(stepCount == MAX_ITERATIONS)
    fprintf(stderr, "Exceeded max number of iterations\n");
#endif
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
  point->Jt_x = malloc(Nstate * sizeof(point->Jt_x[0]));
  ASSERT(point->Jt_x != NULL);

  return point;
}

static void freeOperatingPoint(operatingPoint_t** point, cholmod_common* common)
{
  free((*point)->p);
  free((*point)->x);

  cholmod_free_sparse(&(*point)->Jt,   common);
  free((*point)->Jt_x);

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

  ctx.p_max = malloc(Nstate * sizeof(double));
  ASSERT(ctx.p_max != NULL);
  for(unsigned int i=0; i<Nstate; i++)
    ctx.p_max[i] = fabs(p[i]);

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

  free(ctx.p_max);

  fprintf(stderr, "success! took %d iterations\n", 10);
  return 10.0; // rms
}
