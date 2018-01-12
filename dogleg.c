// -*- mode: C; c-basic-offset: 2 -*-
// Copyright 2011 Oblong Industries
//           2017 Dima Kogan <dima@secretsauce.net>
// License: GNU LGPL <http://www.gnu.org/licenses>.

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#include "dogleg.h"

#if (CHOLMOD_VERSION > (CHOLMOD_VER_CODE(2,2)))
#include <suitesparse/cholmod_function.h>
#endif


#define SAY_NONEWLINE(fmt, ...) fprintf(stderr, "libdogleg at %s:%d: " fmt, __FILE__, __LINE__, ## __VA_ARGS__)
#define SAY(fmt, ...)           do {  SAY_NONEWLINE(fmt, ## __VA_ARGS__); fprintf(stderr, "\n"); } while(0)

// I do this myself because I want this to be active in all build modes, not just !NDEBUG
#define ASSERT(x) do { if(!(x)) { SAY("ASSERTION FAILED: " #x "is not true"); exit(1); } } while(0)

// used to consolidate the casts
#define P(A, index) ((unsigned int*)((A)->p))[index]
#define I(A, index) ((unsigned int*)((A)->i))[index]
#define X(A, index) ((double*      )((A)->x))[index]


//////////////////////////////////////////////////////////////////////////////////////////
// parameter stuff
//////////////////////////////////////////////////////////////////////////////////////////
// These are the optimizer parameters. They have semi-arbitrary defaults. The
// user should adjust them through the API
static int DOGLEG_DEBUG   = 0;
static int MAX_ITERATIONS = 100;

// it is cheap to reject a too-large trust region, so I start with something
// "large". The solver will quickly move down to something reasonable. Only the
// user really knows if this is "large" or not, so they should change this via
// the API
static double TRUSTREGION0 = 1.0e3;

// These are probably OK to leave alone. Tweaking them can maybe result in
// slightly faster convergence
static double TRUSTREGION_DECREASE_FACTOR    = 0.1;
static double TRUSTREGION_INCREASE_FACTOR    = 2;
static double TRUSTREGION_INCREASE_THRESHOLD = 0.75;
static double TRUSTREGION_DECREASE_THRESHOLD = 0.25;

// The termination thresholds. Documented in the header
static double JT_X_THRESHOLD        = 1e-8;
static double UPDATE_THRESHOLD      = 1e-8;
static double TRUSTREGION_THRESHOLD = 1e-8;

// if I ever see a singular JtJ, I factor JtJ + LAMBDA*I from that point on
#define LAMBDA_INITIAL 1e-10

// these parameters likely should be messed with
void dogleg_setDebug(int debug)
{
  DOGLEG_DEBUG = debug;
}
void dogleg_setInitialTrustregion(double t)
{
  TRUSTREGION0 = t;
}
void dogleg_setThresholds(double Jt_x, double update, double trustregion)
{
  if(Jt_x > 0.0)        JT_X_THRESHOLD        = Jt_x;
  if(update > 0.0)      UPDATE_THRESHOLD      = update;
  if(trustregion > 0.0) TRUSTREGION_THRESHOLD = trustregion;
}

// these parameters likely should not be messed with.
void dogleg_setMaxIterations(int n)
{
  MAX_ITERATIONS = n;
}
void dogleg_setTrustregionUpdateParameters(double downFactor, double downThreshold,
                                           double upFactor,   double upThreshold)
{
  TRUSTREGION_DECREASE_FACTOR    = downFactor;
  TRUSTREGION_DECREASE_THRESHOLD = downThreshold;
  TRUSTREGION_INCREASE_FACTOR    = upFactor;
  TRUSTREGION_INCREASE_THRESHOLD = upThreshold;
}




//////////////////////////////////////////////////////////////////////////////////////////
// general boring linear algebra stuff
// should really come from BLAS or libminimath
//////////////////////////////////////////////////////////////////////////////////////////
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
__attribute__((unused))
static double inner_withstride(const double* x, const double* y, unsigned int n, unsigned int stride)
{
  double result = 0;
  for(unsigned int i=0; i<n*stride; i+=stride)
    result += x[i]*y[i];
  return result;
}
// JtJ += outer(j,j). JtJ is packed, upper-triangular (lower-triangular as far
// as LAPACK is concerned)
static void accum_outerproduct_packed( double* JtJ, const double* j, int n )
{
  int iJtJ=0;
  for(int i1=0; i1<n; i1++)
    for(int i0=i1; i0<n; i0++, iJtJ++)
      JtJ[iJtJ] += j[i0]*j[i1];
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
static void vec_negate(double* v, int n)
{
  for(int i=0; i<n; i++)
    v[i] *= -1.0;
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

// transpose(A)*x
static void mul_matrix_t_densevector(double* dest,
                                     const double* A, const double* x,
                                     int Nrows, int Ncols)
{
  memset(dest, 0, sizeof(double) * Ncols);
  for(int i=0; i<Ncols; i++)
    for(int j=0; j<Nrows; j++)
      dest[i] += A[i + j*Ncols]*x[j];
}
static double norm2_mul_matrix_vector(const double* A, const double* x, int Nrows, int Ncols)
{
  // computes norm2(A * x). For each row of A I compute that element of A*x, and
  // accumulate it immediately towards the norm
  double result = 0.0;

  for(int i=0; i<Nrows; i++)
  {
    double dotproduct = inner(x, &A[i*Ncols], Ncols);
    result += dotproduct * dotproduct;
  }

  return result;
}

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

  // This gradient is not in my sparse matrix, so its value is 0.0
  return 0.0;
}
static double getGrad_dense(unsigned int var, int meas, const double* J, int Nstate)
{
  return J[var + meas*Nstate];
}

static
void _dogleg_testGradient(unsigned int var, const double* p0,
                          unsigned int Nstate, unsigned int Nmeas, unsigned int NJnnz,
                          dogleg_callback_t* f, void* cookie)
{
  int is_sparse = NJnnz > 0;

  double* x0 = malloc(Nmeas  * sizeof(double));
  double* x  = malloc(Nmeas  * sizeof(double));
  double* p  = malloc(Nstate * sizeof(double));
  ASSERT(x0);
  ASSERT(x);
  ASSERT(p);

  memcpy(p, p0, Nstate * sizeof(double));


  cholmod_common _cholmod_common;
  cholmod_sparse* Jt;
  cholmod_sparse* Jt0;
  double* J_dense  = NULL; // setting to NULL to pacify compiler's "uninitialized" warnings
  double* J_dense0 = NULL; // setting to NULL to pacify compiler's "uninitialized" warnings


  printf("# ivar imeasurement gradient_reported gradient_observed error error_relative\n");

  if( is_sparse )
  {
    if( !cholmod_start(&_cholmod_common) )
    {
      SAY( "Couldn't initialize CHOLMOD");
      return;
    }

    Jt  = cholmod_allocate_sparse(Nstate, Nmeas, NJnnz,
                                  1, // sorted
                                  1, // packed,
                                  0, // NOT symmetric
                                  CHOLMOD_REAL,
                                  &_cholmod_common);
    Jt0 = cholmod_allocate_sparse(Nstate, Nmeas, NJnnz,
                                  1, // sorted
                                  1, // packed,
                                  0, // NOT symmetric
                                  CHOLMOD_REAL,
                                  &_cholmod_common);

    p[var] -= GRADTEST_DELTA/2.0;
    (*f)(p, x0, Jt0, cookie);
    p[var] += GRADTEST_DELTA;
    (*f)(p, x,  Jt,  cookie);
    p[var] -= GRADTEST_DELTA/2.0;
  }
  else
  {
    J_dense  = malloc( Nmeas * Nstate * sizeof(J_dense[0]) );
    J_dense0 = malloc( Nmeas * Nstate * sizeof(J_dense[0]) );

    dogleg_callback_dense_t* f_dense = (dogleg_callback_dense_t*)f;
    p[var] -= GRADTEST_DELTA/2.0;
    (*f_dense)(p, x0, J_dense0, cookie);
    p[var] += GRADTEST_DELTA;
    (*f_dense)(p, x,  J_dense,  cookie);
    p[var] -= GRADTEST_DELTA/2.0;
  }


  for(unsigned int i=0; i<Nmeas; i++)
  {
    // estimated gradients at the midpoint between x and x0
    double g_observed = (x[i] - x0[i]) / GRADTEST_DELTA;
    double g_reported;
    if( is_sparse )
      g_reported = (getGrad(var, i, Jt0) + getGrad(var, i, Jt)) / 2.0;
    else
      g_reported = (getGrad_dense(var, i, J_dense0, Nstate) + getGrad_dense(var, i, J_dense, Nstate)) / 2.0;

    double g_sum_abs = fabs(g_reported) + fabs(g_observed);
    double g_abs_err = fabs(g_reported - g_observed);

    printf( "%d %d %.6g %.6g %.6g %.6g\n", var, i,
            g_reported, g_observed, g_abs_err,

            g_sum_abs == 0.0 ? 0.0 : (g_abs_err / ( g_sum_abs / 2.0 )));
  }

  if( is_sparse )
  {
    cholmod_free_sparse(&Jt,  &_cholmod_common);
    cholmod_free_sparse(&Jt0, &_cholmod_common);
    cholmod_finish(&_cholmod_common);
  }
  else
  {
    free(J_dense);
    free(J_dense0);
  }

  free(x0);
  free(x);
  free(p);
}
void dogleg_testGradient(unsigned int var, const double* p0,
                         unsigned int Nstate, unsigned int Nmeas, unsigned int NJnnz,
                         dogleg_callback_t* f, void* cookie)
{
  if( NJnnz == 0 )
  {
    SAY( "I must have NJnnz > 0, instead I have %d", NJnnz);
    return;
  }
  return _dogleg_testGradient(var, p0, Nstate, Nmeas, NJnnz, f, cookie);
}
void dogleg_testGradient_dense(unsigned int var, const double* p0,
                               unsigned int Nstate, unsigned int Nmeas,
                               dogleg_callback_dense_t* f, void* cookie)
{
  return _dogleg_testGradient(var, p0, Nstate, Nmeas, 0, (dogleg_callback_t*)f, cookie);
}


//////////////////////////////////////////////////////////////////////////////////////////
// solver routines
//////////////////////////////////////////////////////////////////////////////////////////

static void computeCauchyUpdate(dogleg_operatingPoint_t* point,
                                const dogleg_solverContext_t* ctx)
{
  // I already have this data, so don't need to recompute
  if(point->updateCauchy_valid)
    return;
  point->updateCauchy_valid = 1;

  // I look at a step in the steepest direction that minimizes my
  // quadratic error function (Cauchy point). If this is past my trust region,
  // I move as far as the trust region allows along the steepest descent
  // direction. My error function is F=norm2(f(p)). dF/dP = 2*ft*J.
  // This is proportional to Jt_x, which is thus the steepest ascent direction.
  //
  // Thus along this direction we have F(k) = norm2(f(p + k Jt_x)). The Cauchy
  // point is where F(k) is at a minimum:
  // dF_dk = 2 f(p + k Jt_x)t  J Jt_x ~ (x + k J Jt_x)t J Jt_x =
  // = xt J Jt x + k xt J Jt J Jt x = norm2(Jt x) + k norm2(J Jt x)
  // dF_dk = 0 -> k= -norm2(Jt x) / norm2(J Jt x)
  // Summary:
  // the steepest direction is parallel to Jt*x. The Cauchy point is at
  // k*Jt*x where k = -norm2(Jt*x)/norm2(J*Jt*x)
  double norm2_Jt_x       = norm2(point->Jt_x, ctx->Nstate);
  double norm2_J_Jt_x     = ctx->is_sparse ?
    norm2_mul_spmatrix_t_densevector(point->Jt, point->Jt_x) :
    norm2_mul_matrix_vector         (point->J_dense, point->Jt_x, ctx->Nmeasurements, ctx->Nstate);
  double k                = -norm2_Jt_x / norm2_J_Jt_x;

  point->updateCauchy_lensq = k*k * norm2_Jt_x;

  vec_copy_scaled(point->updateCauchy,
                  point->Jt_x, k, ctx->Nstate);

  if( DOGLEG_DEBUG )
    SAY( "cauchy step size %.6g", sqrt(point->updateCauchy_lensq));
}

// LAPACK prototypes for a packed cholesky factorization and a linear solve
// using that factorization, respectively
int dpptrf_(char* uplo, int* n, double* ap,
            int* info, int uplo_len);
int dpptrs_(char* uplo, int* n, int* nrhs,
            double* ap, double* b, int* ldb, int* info,
            int uplo_len);


void dogleg_computeJtJfactorization(dogleg_operatingPoint_t* point, dogleg_solverContext_t* ctx)
{
  // I already have this data, so don't need to recompute
  if(point->updateGN_valid)
    return;

  if( ctx->is_sparse )
  {
    // I'm assuming the pattern of zeros will remain the same throughout, so I
    // analyze only once
    if(ctx->factorization == NULL)
    {
      ctx->factorization = cholmod_analyze(point->Jt, &ctx->common);
      ASSERT(ctx->factorization != NULL);
    }

    while(1)
    {
      if( ctx->lambda == 0.0 )
        ASSERT( cholmod_factorize(point->Jt, ctx->factorization, &ctx->common) );
      else
      {
        double beta[] = { ctx->lambda, 0 };
        ASSERT( cholmod_factorize_p(point->Jt, beta, NULL, 0,
                                    ctx->factorization, &ctx->common) );
      }

      if(ctx->factorization->minor == ctx->factorization->n)
        break;

      // singular JtJ. Raise lambda and go again
      if( ctx->lambda == 0.0) ctx->lambda = LAMBDA_INITIAL;
      else                    ctx->lambda *= 10.0;
      ASSERT( isfinite(ctx->lambda) );

      if( DOGLEG_DEBUG )
        SAY( "singular JtJ. Have rank/full rank: %zd/%d. Adding %g I from now on",
             ctx->factorization->minor, ctx->Nstate, ctx->lambda);
    }
  }
  else
  {
    if(ctx->factorization_dense == NULL)
    {
      // Need to store symmetric JtJ, so I only need one triangle of it
      ctx->factorization_dense = malloc( ctx->Nstate * (ctx->Nstate+1) / 2 *
                                         sizeof( ctx->factorization_dense[0]));
      ASSERT(ctx->factorization_dense);
    }

    while(1)
    {
      // I construct my JtJ. JtJ is packed and stored row-first. I have two
      // equivalent implementations. The one enabled here is maybe a bit faster,
      // but it's definitely clearer
#if 1
      memset(ctx->factorization_dense,
             0,
             ctx->Nstate*(ctx->Nstate+1)/2*sizeof(ctx->factorization_dense[0]));
      for(int i=0; i<ctx->Nmeasurements; i++)
        accum_outerproduct_packed( ctx->factorization_dense, &point->J_dense[ctx->Nstate*i],
                                   ctx->Nstate );
      if( ctx->lambda > 0.0 )
      {
        int iJtJ=0;
        for(int i1=0; i1<ctx->Nstate; i1++)
        {
          ctx->factorization_dense[iJtJ] += ctx->lambda;
          iJtJ                           += ctx->Nstate-i1;
        }
      }
#else
      int iJtJ = 0;
      for(int i1=0; i1<ctx->Nstate; i1++)
      {
        #error this does not work. overwritten in the following loop
        ctx->factorization_dense[iJtJ] += ctx->lambda;

        for(int i0=i1; i0<ctx->Nstate; i0++, iJtJ++)
          ctx->factorization_dense[iJtJ] = inner_withstride( &point->J_dense[i0],
                                                             &point->J_dense[i1],
                                                             ctx->Nmeasurements,
                                                             ctx->Nstate);
      }
#endif



      int info;
      dpptrf_(&(char){'L'}, &(int){ctx->Nstate}, ctx->factorization_dense,
              &info, 1);
      ASSERT(info >= 0); // we MUST either succeed or see complain of singular
      // JtJ
      if( info == 0 )
        break;

      // singular JtJ. Raise lambda and go again
      if( ctx->lambda == 0.0) ctx->lambda = LAMBDA_INITIAL;
      else                    ctx->lambda *= 10.0;

      if( DOGLEG_DEBUG )
        SAY( "singular JtJ. Adding %g I from now on", ctx->lambda);
    }
  }
}

static void computeGaussNewtonUpdate(dogleg_operatingPoint_t* point, dogleg_solverContext_t* ctx)
{
  // I already have this data, so don't need to recompute
  if(point->updateGN_valid)
    return;

  dogleg_computeJtJfactorization(point, ctx);

  // try to factorize the matrix directly. If it's singular, add a small
  // constant to the diagonal. This constant gets larger if we keep being
  // singular
  if( ctx->is_sparse )
  {
    // solve JtJ*updateGN = Jt*x. Gauss-Newton step is then -updateGN
    cholmod_dense Jt_x_dense = {.nrow  = ctx->Nstate,
                                .ncol  = 1,
                                .nzmax = ctx->Nstate,
                                .d     = ctx->Nstate,
                                .x     = point->Jt_x,
                                .xtype = CHOLMOD_REAL,
                                .dtype = CHOLMOD_DOUBLE};

    if(point->updateGN_cholmoddense != NULL)
      cholmod_free_dense(&point->updateGN_cholmoddense, &ctx->common);

    point->updateGN_cholmoddense = cholmod_solve(CHOLMOD_A,
                                                 ctx->factorization,
                                                 &Jt_x_dense,
                                                 &ctx->common);
    vec_negate(point->updateGN_cholmoddense->x,
               ctx->Nstate); // should be more efficient than this later

    point->updateGN_lensq = norm2(point->updateGN_cholmoddense->x, ctx->Nstate);
  }
  else
  {
    memcpy( point->updateGN_dense, point->Jt_x, ctx->Nstate * sizeof(point->updateGN_dense[0]));
    int info;
    dpptrs_(&(char){'L'}, &(int){ctx->Nstate}, &(int){1},
            ctx->factorization_dense,
            point->updateGN_dense, &(int){ctx->Nstate}, &info, 1);
    vec_negate(point->updateGN_dense,
               ctx->Nstate); // should be more efficient than this later

    point->updateGN_lensq = norm2(point->updateGN_dense, ctx->Nstate);
  }


  if( DOGLEG_DEBUG )
    SAY( "gn step size %.6g", sqrt(point->updateGN_lensq));

  point->updateGN_valid = 1;
}

static void computeInterpolatedUpdate(double*                  update_dogleg,
                                      dogleg_operatingPoint_t* point,
                                      double                   trustregion,
                                      const dogleg_solverContext_t* ctx)
{
  // I interpolate between the Cauchy-point step and the Gauss-Newton step
  // to find a step that takes me to the edge of my trust region.
  //
  // I have something like norm2(a + k*(b-a)) = dsq
  // = norm2(a) + 2*at*(b-a) * k + norm2(b-a)*k^2 = dsq
  // let c = at*(b-a), l2 = norm2(b-a) ->
  // l2 k^2 + 2*c k + norm2(a)-dsq = 0
  //
  // This is a simple quadratic equation:
  // k = (-2*c +- sqrt(4*c*c - 4*l2*(norm2(a)-dsq)))/(2*l2)
  //   = (-c +- sqrt(c*c - l2*(norm2(a)-dsq)))/l2

  // to make 100% sure the discriminant is positive, I choose a to be the
  // cauchy step.  The solution must have k in [0,1], so I much have the
  // +sqrt side, since the other one is negative
  double        dsq    = trustregion*trustregion;
  double        norm2a = point->updateCauchy_lensq;
  const double* a      = point->updateCauchy;
  const double* b      = ctx->is_sparse ? point->updateGN_cholmoddense->x : point->updateGN_dense;

  double l2    = 0.0;
  double neg_c = 0.0;
  for(int i=0; i<ctx->Nstate; i++)
  {
    double d = a[i] - b[i];

    l2    += d*d;
    neg_c += d*a[i];
  }

  double discriminant = neg_c*neg_c - l2* (norm2a - dsq);
  if(discriminant < 0.0)
  {
    SAY( "negative discriminant: %.6g!", discriminant);
    discriminant = 0.0;
  }
  double k = (neg_c + sqrt(discriminant))/l2;

  for(int i=0; i<ctx->Nstate; i++)
    update_dogleg[i] = a[i] + k*(b[i] - a[i]);

  if( DOGLEG_DEBUG )
  {
    double updateNorm = norm2(update_dogleg, ctx->Nstate);
    SAY( "k_cauchy_to_gn %.6g, norm %.6g", k, sqrt(updateNorm));
  }
}

// takes in point->p, and computes all the quantities derived from it, storing
// the result in the other members of the operatingPoint structure. Returns
// true if the gradient-size termination criterion has been met
static int computeCallbackOperatingPoint(dogleg_operatingPoint_t* point, dogleg_solverContext_t* ctx)
{
  if( ctx->is_sparse )
  {
    (*ctx->f)(point->p, point->x, point->Jt, ctx->cookie);

    // compute Jt*x
    mul_spmatrix_densevector(point->Jt_x, point->Jt, point->x);
  }
  else
  {
    (*ctx->f_dense)(point->p, point->x, point->J_dense, ctx->cookie);

    // compute Jt*x
    mul_matrix_t_densevector(point->Jt_x, point->J_dense, point->x,
                             ctx->Nmeasurements, ctx->Nstate);
  }

  // I just got a new operating point, so the current update vectors aren't
  // valid anymore, and should be recomputed, as needed
  point->updateCauchy_valid = 0;
  point->updateGN_valid     = 0;

  // compute the 2-norm of the current error vector
  // At some point this should be changed to use the call from libminimath
  point->norm2_x = norm2(point->x, ctx->Nmeasurements);

  // If the largest absolute gradient element is smaller than the threshold,
  // we can stop iterating. This is equivalent to the inf-norm
  for(int i=0; i<ctx->Nstate; i++)
    if(fabs(point->Jt_x[i]) > JT_X_THRESHOLD)
      return 0;
  if( DOGLEG_DEBUG )
    SAY( "Jt_x all below the threshold. Done iterating!");

  return 1;
}
static double computeExpectedImprovement(const double* step, const dogleg_operatingPoint_t* point,
                                         const dogleg_solverContext_t* ctx)
{
  // My error function is F=norm2(f(p + step)). F(0) - F(step) =
  // = norm2(x) - norm2(x + J*step) = -2*inner(x,J*step) - norm2(J*step)
  // = -2*inner(Jt_x,step) - norm2(J*step)
  if( ctx->is_sparse )
    return
      - 2.0*inner(point->Jt_x, step, ctx->Nstate)
      - norm2_mul_spmatrix_t_densevector(point->Jt, step);
  else
    return
      - 2.0*inner(point->Jt_x, step, ctx->Nstate)
      - norm2_mul_matrix_vector(point->J_dense, step, ctx->Nmeasurements, ctx->Nstate);
}


// takes a step from the given operating point, using the given trust region
// radius. Returns the expected improvement, based on the step taken and the
// linearized x(p). If we can stop iterating, returns a negative number
static double takeStepFrom(dogleg_operatingPoint_t* pointFrom, double* newp,
                           double trustregion, dogleg_solverContext_t* ctx)
{
  if( DOGLEG_DEBUG )
    SAY( "taking step with trustregion %.6g", trustregion);

  double update_array[ctx->Nstate];
  double* update;


  computeCauchyUpdate(pointFrom, ctx);

  if(pointFrom->updateCauchy_lensq >= trustregion*trustregion)
  {
    if( DOGLEG_DEBUG )
      SAY( "taking cauchy step");

    // cauchy step goes beyond my trust region, so I do a gradient descent
    // to the edge of my trust region and call it good
    vec_copy_scaled(update_array,
                    pointFrom->updateCauchy,
                    trustregion / sqrt(pointFrom->updateCauchy_lensq),
                    ctx->Nstate);
    update = update_array;
    pointFrom->didStepToEdgeOfTrustRegion = 1;
  }
  else
  {
    // I'm not yet done. The cauchy point is within the trust region, so I can
    // go further. I look at the full Gauss-Newton step. If this is within the
    // trust region, I use it. Otherwise, I find the point at the edge of my
    // trust region that lies on a straight line between the Cauchy point and
    // the Gauss-Newton solution, and use that. This is the heart of Powell's
    // dog-leg algorithm.
    computeGaussNewtonUpdate(pointFrom, ctx);
    if(pointFrom->updateGN_lensq <= trustregion*trustregion)
    {
      if( DOGLEG_DEBUG )
        SAY( "taking GN step");

      // full Gauss-Newton step lies within my trust region. Take the full step
      update = ctx->is_sparse ? pointFrom->updateGN_cholmoddense->x : pointFrom->updateGN_dense;
      pointFrom->didStepToEdgeOfTrustRegion = 0;
    }
    else
    {
      if( DOGLEG_DEBUG )
        SAY( "taking interpolated step");

      // full Gauss-Newton step lies outside my trust region, so I interpolate
      // between the Cauchy-point step and the Gauss-Newton step to find a step
      // that takes me to the edge of my trust region.
      computeInterpolatedUpdate(update_array, pointFrom, trustregion, ctx);
      update = update_array;
      pointFrom->didStepToEdgeOfTrustRegion = 1;
    }
  }



  // take the step
  vec_add(newp, pointFrom->p, update, ctx->Nstate);
  double expectedImprovement = computeExpectedImprovement(update, pointFrom, ctx);

  // are we done? For each state variable I look at the update step. If all the elements fall below
  // a threshold, I call myself done
  for(int i=0; i<ctx->Nstate; i++)
    if( fabs(update[i]) > UPDATE_THRESHOLD )
      return expectedImprovement;

  if( DOGLEG_DEBUG )
    SAY( "update small enough. Done iterating!");

  return -1.0;
}


// I have a candidate step. I adjust the trustregion accordingly, and also
// report whether this step should be accepted (0 == rejected, otherwise
// accepted)
static int evaluateStep_adjustTrustRegion(const dogleg_operatingPoint_t* before,
                                          const dogleg_operatingPoint_t* after,
                                          double* trustregion,
                                          double expectedImprovement)
{
  double observedImprovement = before->norm2_x - after->norm2_x;

  double rho = observedImprovement / expectedImprovement;
  if( DOGLEG_DEBUG )
  {
    SAY( "observed/expected improvement: %.6g/%.6g. rho = %.6g",
            observedImprovement, expectedImprovement, rho);
  }


  // adjust the trust region
  if( rho < TRUSTREGION_DECREASE_THRESHOLD )
  {
    if( DOGLEG_DEBUG )
      SAY( "rho too small. decreasing trust region");

    // Our model doesn't fit well. We should reduce the trust region size. If
    // the trust region size was affecting the attempted step, do this by a
    // constant factor. Otherwise, drop the trustregion to attempted step size
    // first
    if( !before->didStepToEdgeOfTrustRegion )
      *trustregion = sqrt(before->updateGN_lensq);

    *trustregion *= TRUSTREGION_DECREASE_FACTOR;
  }
  else if (rho > TRUSTREGION_INCREASE_THRESHOLD && before->didStepToEdgeOfTrustRegion)
  {
    if( DOGLEG_DEBUG )
      SAY( "rho large enough. increasing trust region");

    *trustregion *= TRUSTREGION_INCREASE_FACTOR;
  }

  return rho > 0.0;
}

static int runOptimizer(dogleg_solverContext_t* ctx)
{
  double trustregion = TRUSTREGION0;
  int stepCount = 0;

  if( computeCallbackOperatingPoint(ctx->beforeStep, ctx) )
    return stepCount;

  if( DOGLEG_DEBUG )
    SAY( "Initial operating point has norm2_x %.6g", ctx->beforeStep->norm2_x);


  while( stepCount<MAX_ITERATIONS )
  {
    if( DOGLEG_DEBUG )
    {
      SAY( "================= step %d", stepCount );
    }

    while(1)
    {
      if( DOGLEG_DEBUG )
        SAY("--------");

      double expectedImprovement =
        takeStepFrom(ctx->beforeStep, ctx->afterStep->p, trustregion, ctx);

      // negative expectedImprovement is used to indicate that we're done
      if(expectedImprovement < 0.0)
        return stepCount;

      int afterStepZeroGradient = computeCallbackOperatingPoint(ctx->afterStep, ctx);
      if( DOGLEG_DEBUG )
        SAY( "Evaluated operating point with norm2_x %.6g", ctx->afterStep->norm2_x);

      if( evaluateStep_adjustTrustRegion(ctx->beforeStep, ctx->afterStep, &trustregion,
                                         expectedImprovement) )
      {
        if( DOGLEG_DEBUG )
          SAY( "accepted step");

        stepCount++;

        // I accept this step, so the after-step operating point is the before-step operating point
        // of the next iteration. I exchange the before- and after-step structures so that all the
        // pointers are still around and I don't have to re-allocate
        dogleg_operatingPoint_t* tmp;
        tmp             = ctx->afterStep;
        ctx->afterStep  = ctx->beforeStep;
        ctx->beforeStep = tmp;

        if( afterStepZeroGradient )
        {
          if( DOGLEG_DEBUG )
            SAY( "Gradient low enough and we just improved. Done iterating!" );

          return stepCount;
        }

        break;
      }

      if( DOGLEG_DEBUG )
        SAY( "rejected step");

      // This step was rejected. check if the new trust region size is small
      // enough to give up
      if(trustregion < TRUSTREGION_THRESHOLD)
      {
        if( DOGLEG_DEBUG )
          SAY( "trust region small enough. Giving up. Done iterating!");

        return stepCount;
      }

      // I have rejected this step, so I try again with the new trust region
    }
  }

  if( DOGLEG_DEBUG && stepCount == MAX_ITERATIONS)
      SAY( "Exceeded max number of iterations");

  return stepCount;
}

static
dogleg_operatingPoint_t* allocOperatingPoint(int      Nstate, int numMeasurements,
                                             unsigned int NJnnz,
                                             cholmod_common* common)
{
  int is_sparse = NJnnz > 0;

  dogleg_operatingPoint_t* point = malloc(sizeof(dogleg_operatingPoint_t));
  ASSERT(point != NULL);

  point->p = malloc(Nstate * sizeof(point->p[0]));
  ASSERT(point->p != NULL);

  point->x = malloc(numMeasurements * sizeof(point->x[0]));
  ASSERT(point->x != NULL);

  if( is_sparse )
  {
    point->Jt = cholmod_allocate_sparse(Nstate, numMeasurements, NJnnz,
                                        1, // sorted
                                        1, // packed,
                                        0, // NOT symmetric
                                        CHOLMOD_REAL,
                                        common);
    ASSERT(point->Jt != NULL);
    point->updateGN_cholmoddense = NULL; // This will be allocated as it is used
  }
  else
  {
    point->J_dense = malloc( numMeasurements * Nstate * sizeof(point->J_dense[0]) );
    ASSERT(point->J_dense != NULL);

    point->updateGN_dense = malloc( Nstate * sizeof(point->updateGN_dense[0]) );
    ASSERT(point->updateGN_dense != NULL);
  }

  // the 1-column vector Jt * x
  point->Jt_x = malloc(Nstate * sizeof(point->Jt_x[0]));
  ASSERT(point->Jt_x != NULL);

  // the cached update vectors
  point->updateCauchy          = malloc(Nstate * sizeof(point->updateCauchy[0]));

  // vectors don't have any valid data yet
  point->updateCauchy_valid = 0;
  point->updateGN_valid     = 0;

  return point;
}

static void freeOperatingPoint(dogleg_operatingPoint_t** point, cholmod_common* common)
{
  free((*point)->p);
  free((*point)->x);

  int is_sparse = common != NULL;

  if( is_sparse )
  {
    cholmod_free_sparse(&(*point)->Jt,   common);

    if((*point)->updateGN_cholmoddense != NULL)
      cholmod_free_dense(&(*point)->updateGN_cholmoddense, common);
  }
  else
  {
    free( (*point)->J_dense );
    free( (*point)->updateGN_dense );
  }

  free((*point)->Jt_x);
  free((*point)->updateCauchy);

  free(*point);
  *point = NULL;
}

static int cholmod_error_callback(const char* s, ...)
{
  SAY_NONEWLINE("");

  va_list ap;
  va_start(ap, s);
  int ret = vfprintf(stderr, s, ap);
  va_end(ap);

  return ret;
}

static void set_cholmod_options(cholmod_common* cc)
{
  // I want to use LGPL parts of CHOLMOD only, so I turn off the supernodal routines. This gave me a
  // 25% performance hit in the solver for a particular set of optical calibration data.
  cc->supernodal = 0;


  // I want all output to go to STDERR, not STDOUT
#if (CHOLMOD_VERSION <= (CHOLMOD_VER_CODE(2,2)))
  cc->print_function = cholmod_error_callback;
#else
  CHOLMOD_FUNCTION_DEFAULTS ;
  CHOLMOD_FUNCTION_PRINTF(cc) = cholmod_error_callback;
#endif
}

void dogleg_freeContext(dogleg_solverContext_t** ctx)
{
  freeOperatingPoint(&(*ctx)->beforeStep, (*ctx)->is_sparse ? &(*ctx)->common : NULL);
  freeOperatingPoint(&(*ctx)->afterStep,  (*ctx)->is_sparse ? &(*ctx)->common : NULL);

  if( (*ctx)->is_sparse )
  {
    if((*ctx)->factorization != NULL)
      cholmod_free_factor (&(*ctx)->factorization, &(*ctx)->common);
    cholmod_finish(&(*ctx)->common);
  }
  else
    free((*ctx)->factorization_dense);

  free(*ctx);
  *ctx = NULL;
}

static double _dogleg_optimize(double* p, unsigned int Nstate,
                               unsigned int Nmeas, unsigned int NJnnz,
                               dogleg_callback_t* f,
                               void* cookie,
                               dogleg_solverContext_t** returnContext)
{
  int is_sparse = NJnnz > 0;


  dogleg_solverContext_t* ctx = malloc(sizeof(dogleg_solverContext_t));
  ctx->f                      = f;
  ctx->cookie                 = cookie;
  ctx->factorization          = NULL;
  ctx->lambda                 = 0.0;
  ctx->Nstate                 = Nstate;
  ctx->Nmeasurements          = Nmeas;


  if( returnContext != NULL )
    *returnContext = ctx;

  if( is_sparse )
  {
    if( !cholmod_start(&ctx->common) )
    {
      SAY( "Couldn't initialize CHOLMOD");
      return -1.0;
    }
    set_cholmod_options(&ctx->common);
    ctx->is_sparse = 1;
  }
  else
    ctx->is_sparse = 0;

  ctx->beforeStep = allocOperatingPoint(Nstate, Nmeas, NJnnz, &ctx->common);
  ctx->afterStep  = allocOperatingPoint(Nstate, Nmeas, NJnnz, &ctx->common);

  memcpy(ctx->beforeStep->p, p, Nstate * sizeof(double));

  // everything is set up, so run the solver!
  int    numsteps = runOptimizer(ctx);
  double norm2_x  = ctx->beforeStep->norm2_x;

  // runOptimizer places the most recent results into beforeStep in preparation for another
  // iteration
  memcpy(p, ctx->beforeStep->p, Nstate * sizeof(double));

  if( returnContext == NULL )
    dogleg_freeContext(&ctx);

  if( DOGLEG_DEBUG )
    SAY( "success! took %d iterations", numsteps);

  return norm2_x;
}

double dogleg_optimize(double* p, unsigned int Nstate,
                       unsigned int Nmeas, unsigned int NJnnz,
                       dogleg_callback_t* f,
                       void* cookie,
                       dogleg_solverContext_t** returnContext)
{
  if( NJnnz == 0 )
  {
    SAY( "I must have NJnnz > 0, instead I have %d", NJnnz);
    return -1.0;
  }

  return _dogleg_optimize(p, Nstate, Nmeas, NJnnz, f,
                          cookie, returnContext);
}


double dogleg_optimize_dense(double* p, unsigned int Nstate,
                             unsigned int Nmeas,
                             dogleg_callback_dense_t* f, void* cookie,
                             dogleg_solverContext_t** returnContext)
{
  return _dogleg_optimize(p, Nstate, Nmeas, 0, (dogleg_callback_t*)f,
                          cookie, returnContext);
}
