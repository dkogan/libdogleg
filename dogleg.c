// -*- mode: C; c-basic-offset: 2 -*-
// Copyright 2011 Oblong Industries
//           2017-2018 Dima Kogan <dima@secretsauce.net>
// License: GNU LGPL <http://www.gnu.org/licenses>.

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <float.h>
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
#define SAY_IF_VERBOSE(fmt,...) do { if( DOGLEG_DEBUG ) SAY(fmt, ##__VA_ARGS__); } while(0)

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


  // This is a plain text table, that can be easily parsed with "vnlog" tools
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







// Computes pinv(J) for a subset of measurements: inv(JtJ) *
// Jt[imeasurement0..imeasurement0+N-1]. Returns false if something failed.
// ASSUMES THAT THE CHOLESKY FACTORIZATION HAS ALREADY BEEN COMPUTED.
//
// This function is experimental, and subject to change
static bool pseudoinverse_J_dense(// output
                                  double* out,

                                  // inputs
                                  const dogleg_operatingPoint_t* point,
                                  const dogleg_solverContext_t* ctx,
                                  int i_meas0, int NmeasInChunk)
{
  int info;
  memcpy(out,
         &point->J_dense[i_meas0*ctx->Nstate],
         NmeasInChunk*ctx->Nstate*sizeof(double));
  dpptrs_(&(char){'L'}, &(int){ctx->Nstate}, &NmeasInChunk,
          ctx->factorization_dense,
          out,
          &(int){ctx->Nstate}, &info, 1);
  return info==0;
}

// Computes pinv(J) for a subset of measurements: inv(JtJ) *
// Jt[imeasurement0..imeasurement0+N-1]. Returns false if something failed.
// ASSUMES THAT THE CHOLESKY FACTORIZATION HAS ALREADY BEEN COMPUTED.
//
// allocates memory, returns NULL on failure. ON SUCCESS, THE CALLER IS
// RESPONSIBLE FOR FREEING THE RETURNED MEMORY
//
// This function is experimental, and subject to change
static cholmod_dense* pseudoinverse_J_sparse(// inputs
                                             const dogleg_operatingPoint_t* point,
                                             dogleg_solverContext_t* ctx,
                                             int i_meas0, int NmeasInChunk,

                                             // Pre-allocated array for the
                                             // right-hand-side. This will be used as a
                                             // workspace. Create this like so:
                                             //
                                             //   cholmod_allocate_dense( Nstate,
                                             //                           NmeasInChunk,
                                             //                           Nstate,
                                             //                           CHOLMOD_REAL,
                                             //                           &ctx->common );

                                             cholmod_dense* Jt_chunk)
{
  // I'm solving JtJ x = b where J is sparse, b is sparse, but x ends up dense.
  // cholmod doesn't have functions for this exact case. so I use the
  // dense-sparse-dense function (cholmod_solve), and densify the input. Instead
  // of sparse-sparse-sparse and the densifying the output (cholmod_spsolve).
  // This feels like it'd be more efficient

  memset( Jt_chunk->x, 0, Jt_chunk->nrow*Jt_chunk->ncol*sizeof(double) );
  int Jt_chunk_ncol_backup = Jt_chunk->ncol;
  for(int i_meas=0; i_meas<NmeasInChunk; i_meas++)
  {
    if( i_meas0 + i_meas >= (int)ctx->Nmeasurements )
    {
      // at the end, we could have one chunk with less that chunk_size
      // columns
      Jt_chunk->ncol = i_meas;
      break;
    }

    for(unsigned int i0=P(point->Jt, i_meas0+i_meas); i0<P(point->Jt, i_meas0+i_meas+1); i0++)
    {
      int irow = I(point->Jt,i0);
      double x = X(point->Jt,i0);
      ((double*)Jt_chunk->x)[irow + i_meas*Jt_chunk->nrow] = x;
    }
  }

  // solve inv(JtJ)Jt[slice]
  cholmod_dense* pinv =
    cholmod_solve(CHOLMOD_A,
                  ctx->factorization,
                  Jt_chunk,
                  &ctx->common);
  Jt_chunk->ncol = Jt_chunk_ncol_backup;
  return pinv;
}

static void accum_outlierness_factor(// output
                                     double* factor,

                                     // inputs
                                     const double* x,

                                     // A is symmetric. I store the upper triangle
                                     const double* A,

                                     // if outliers are grouped into sets, the
                                     // group size is stated here
                                     int measurementGroupSize,
                                     int Nmeasurements)
{
  // from the derivation in a big comment in dogleg_getOutliernessFactors() I
  // haven't implemented anything else yethave:
  //
  //   f = 1/N (xot (I + A - A inv(-I + A) A) xo )
  //
  // where A = Jo inv(JtJ) Jot

  // I only implemented measurementGroupSize so far
  if(measurementGroupSize <= 1)
  {
    measurementGroupSize = 1;

    double denom = 1.0 - *A;
    if( fabs(denom) < 1e-8 )
      *factor = DBL_MAX; // definitely an outlier
    else
      *factor = x[0]*x[0] / (denom * (double)Nmeasurements);
  }
  else if(measurementGroupSize == 2)
  {
    //   f = 1/N (xot (I + A - A inv(-I + A) A) xo ) =
    //     = 1/N (norm2(xo) + xot A xo - (A xo)t inv(-I + A) (A xo)) =

    double det = (A[0]-1.0)*(A[2]-1.0) - A[1]*A[1];
    if( fabs(det) < 1e-8 )
      *factor = DBL_MAX; // definitely an outlier
    else
    {
      *factor = x[0]*x[0] + x[1]*x[1];

      double Ax[] = {A[0]*x[0] + A[1]*x[1],
                     A[1]*x[0] + A[2]*x[1]};
      *factor += Ax[0]*x[0] + Ax[1]*x[1];

      double inv_A1_Ax_det[] = { (A[2]-1.0)*Ax[0] - A[1]      *Ax[1],
                                 -A[1]     *Ax[0] + (A[0]-1.0)*Ax[1] };

      *factor -= (inv_A1_Ax_det[0]*Ax[0] + inv_A1_Ax_det[1]*Ax[1]) / det;
      *factor /= (double)Nmeasurements;
    }
  }
  else
  {
    SAY("measurementGroupSize > 2 not implemented yet. Got measurementGroupSize=%d", measurementGroupSize);
    ASSERT(0);
  }
}

static bool getOutliernessFactors_dense( // output
                                        double* factors, // Ngroups factors

                                        // inputs
                                        // if outliers are grouped into sets,
                                        // the group size is stated here
                                        int measurementGroupSize,
                                        int Ngroups,
                                        const dogleg_operatingPoint_t* point,
                                        dogleg_solverContext_t* ctx )
{
  // cholmod_spsolve() and cholmod_solve()) work in chunks of 4, so I do this in
  // chunks of 4 too. I pass it rows of J, 4 at a time. Note that if I have
  // measurement groups, I don't want these to cross chunk boundaries, so I set
  // up chunk_size=lcm(N,4)
  int chunk_size = 4;
  if(measurementGroupSize <= 1)
    measurementGroupSize = 1;
  if(measurementGroupSize > 1)
  {
    // haven't implemented anything else yet. Don't have lcm() readily available
    ASSERT(measurementGroupSize == 2);
    // chunk_size = lcm(chunk_size,measurementGroupSize)
  }

  int  Nstate        = ctx->Nstate;
  int  Nmeasurements = ctx->Nmeasurements;
  bool result        = false;

  double* invJtJ_Jt = malloc(Nstate*chunk_size*sizeof(double));
  if(invJtJ_Jt == NULL)
  {
    SAY("Couldn't allocate invJtJ_Jt!");
    goto done;
  }

  int i_measurement_valid_chunk_start = -1;
  int i_measurement_valid_chunk_last  = -1;
  int i_measurement = 0;
  for( int i_group=0; i_group<Ngroups; i_group++, i_measurement+=measurementGroupSize)
  {
    if( i_measurement > i_measurement_valid_chunk_last )
    {
      bool pinvresult = pseudoinverse_J_dense(invJtJ_Jt, point, ctx,
                                              i_measurement, chunk_size);
      if(!pinvresult)
      {
        SAY("Couldn't compute pinv!");
        goto done;
      }
      i_measurement_valid_chunk_start = i_measurement;
      i_measurement_valid_chunk_last  = i_measurement+chunk_size-1;
    }

    // from the derivation in a big comment in dogleg_getOutliernessFactors() I
    // haven't implemented anything else yethave:
    //
    //   f = 1/N (xot (I + A - A inv(-I + A) A) xo )
    //
    // where A = Jo inv(JtJ) Jot
    //
    // A is symmetric. I store the upper triangle
    double A[measurementGroupSize*(measurementGroupSize+1)/2];
    int iA=0;
    for(int i=0; i<measurementGroupSize; i++)
      for(int j=i; j<measurementGroupSize; j++, iA++)
      {
        A[iA] = 0.0;

        for(int k=0; k<Nstate; k++)
          A[iA] +=
            invJtJ_Jt     [Nstate*(i_measurement+i -i_measurement_valid_chunk_start) + k] *
            point->J_dense[Nstate* i_measurement+j                                   + k];
      }
    accum_outlierness_factor(&factors[i_group],
                             &point->x[i_measurement],
                             A, measurementGroupSize, Nmeasurements);
  }

  result = true;
 done:
  free(invJtJ_Jt);
  return result;
}


static bool getOutliernessFactors_sparse( // output
                                         double* factors, // Ngroups factors

                                         // inputs
                                         // if outliers are grouped into sets,
                                         // the group size is stated here
                                         int measurementGroupSize,
                                         int Ngroups,
                                         const dogleg_operatingPoint_t* point,
                                         dogleg_solverContext_t* ctx )
{
  // cholmod_spsolve() and cholmod_solve()) work in chunks of 4, so I do this in
  // chunks of 4 too. I pass it rows of J, 4 at a time. Note that if I have
  // measurement groups, I don't want these to cross chunk boundaries, so I set
  // up chunk_size=lcm(N,4)
  int chunk_size = 4;
  if(measurementGroupSize <= 1)
    measurementGroupSize = 1;
  if(measurementGroupSize > 1)
  {
    // haven't implemented anything else yet. Don't have lcm() readily available
    ASSERT(measurementGroupSize == 2);
    // chunk_size = lcm(chunk_size,measurementGroupSize)
  }

  int  Nstate        = ctx->Nstate;
  int  Nmeasurements = ctx->Nmeasurements;
  bool result        = false;

  cholmod_dense* invJtJ_Jt = NULL;
  cholmod_dense* Jt_chunk =
    cholmod_allocate_dense( Nstate,
                            chunk_size,
                            Nstate,
                            CHOLMOD_REAL,
                            &ctx->common );
  if(!Jt_chunk)
  {
    SAY("Couldn't allocate Jt_chunk!");
    goto done;
  }

  int i_measurement_valid_chunk_start = -1;
  int i_measurement_valid_chunk_last  = -1;
  int i_measurement = 0;
  for( int i_group=0; i_group<Ngroups; i_group++, i_measurement+=measurementGroupSize)
  {
    if( i_measurement > i_measurement_valid_chunk_last )
    {
      if(invJtJ_Jt) cholmod_free_dense(&invJtJ_Jt, &ctx->common);
      invJtJ_Jt = pseudoinverse_J_sparse(point, ctx,
                                         i_measurement, chunk_size,
                                         Jt_chunk);
      if(invJtJ_Jt == NULL)
      {
        SAY("Couldn't compute pinv!");
        goto done;
      }

      i_measurement_valid_chunk_start = i_measurement;
      i_measurement_valid_chunk_last  = i_measurement+chunk_size-1;
    }

    // from the derivation in a big comment in dogleg_getOutliernessFactors() I
    // haven't implemented anything else yethave:
    //
    //   f = 1/N (xot (I + A - A inv(-I + A) A) xo )
    //
    // where A = Jo inv(JtJ) Jot
    //
    // A is symmetric. I store the upper triangle
    double A[measurementGroupSize*(measurementGroupSize+1)/2];
    int iA=0;
    for(int i=0; i<measurementGroupSize; i++)
      for(int j=i; j<measurementGroupSize; j++, iA++)
      {
        A[iA] = 0.0;

        for(unsigned int l = P(point->Jt, i_measurement+j);
            l     < P(point->Jt, i_measurement+j+1);
            l++)
        {
          int k = I(point->Jt, l);
          A[iA] +=
            ((double*)invJtJ_Jt->x)[Nstate*(i_measurement+i-i_measurement_valid_chunk_start) + k] *
            X(point->Jt, l);
        }
      }
    accum_outlierness_factor(&factors[i_group],
                             &point->x[i_measurement],
                             A, measurementGroupSize, Nmeasurements);
  }

  result = true;
 done:
  if(Jt_chunk)  cholmod_free_dense(&Jt_chunk,  &ctx->common);
  if(invJtJ_Jt) cholmod_free_dense(&invJtJ_Jt, &ctx->common);
  return result;
}

// Computes outlierness factors. This function is experimental, and subject to
// change. See comment inside function for detail.
bool dogleg_getOutliernessFactors( // output
                                  double* factors, // Ngroups factors

                                  // inputs
                                  // if outliers are grouped into sets, the group size is
                                  // stated here
                                  int measurementGroupSize,
                                  int Ngroups,
                                  dogleg_operatingPoint_t* point,
                                  dogleg_solverContext_t* ctx )
{
  // We just computed an optimum least-squares fit, and we try to determine if
  // some of the data points look like outliers.
  // The least squares problem I just solved has cost function
  //
  //   E = sum( norm2(x) )
  //
  // where x is a length-N vector of measurements. We solved it by optimizing
  // the vector of parameters p. We define an outlier as a measurement that
  // would greatly improve the MEAN cost function E/N if this measurement was
  // removed, and the problem was re-optimized.
  //
  // Each measurement contributes to the cost function in two ways:
  // - directly, with its x value
  //
  // - indirectly, since minimizing E for THIS measurement pulls p in a
  //   particular direction, and without this measurement, p is more free to
  //   move in a way that reduces errors for the other measurements
  //
  // We look at the total cost of a particular measurement, by looking at the
  // sum of these contributions. Note that this definition of an outlier has a
  // big caveat: measurements that define the problem will look like outliers
  // because taking them out will let the solver overfit to an
  // artificially-low value of the cost function. The way to address this is
  // to add another condition: a set of measurements are outliers if they have
  // a large outlierness factor AND if removing it doesn't cause a large
  // decrease in the confidence of the solution. The former is a very simple
  // computation, performed by this function. The latter is more
  // computationally costly. We only need to run that computation for
  // large-outlierness-factor measurements, so it should be relatively quick.
  //
  // Solving the full problem (BOTH inliers, outliers) I minimize E_io =
  // norm2( x_io(p_io) ). If I knew which measurements are outliers, I can
  // remove them and minimize E_i = norm2( x_i(p_i) ). In reality I don't know
  // which measurements are outliers, so I have ..._io, but not ..._i.
  //
  // Let p_i and p_io be the state vectors are their respective optima. And
  // let's assume that the system is locally linear with J = dx/dp. I can
  // split this into Ji = dxi/dp, Jo = dxo/dp. Then p_i = p_io + dp ->
  // x_i(p_i) ~ x_i(p_io) + Ji*dp -> Ei = norm2( x_i(p_io) + Ji*dp )
  // dp is the shift in p: dp = pi-pio
  //
  // p_io optimizes E_io ->
  // dEio/dp = 0 -> Jt x_io(p_io) = Jit xi + Jot xo = 0
  //     ---> Jit xi = -Jot xo
  //
  // We also know that p_i optimizes E_i ->
  // dEi/dpi = 0 = -> Jit ( x_i(p_io) + Ji dp ) = 0
  //    ---> Jit xi = -Jit Ji dp
  //    ---> dp = - inv(Jit Ji) Jit xi
  //
  // The outlierness factor is E_io/(Ni+No) - E_i/Ni. If Ni >> No ->
  // outlierness factor f = (E_io - E_i) / N. I would expect the cost function
  // to improve when data is removed, so I would expect f > 0, with LARGE f
  // indicating an outlier.
  //
  // f = 1/N (norm2( x_io(p_io) ) - norm2( x_i(p_io) + Ji*dp ) )
  //   = 1/N (norm2( x_i(p_io) ) + norm2( x_o(p_io) ) - norm2( x_i(p_io) + Ji*dp ) )
  //   = 1/N (norm2( x_o(p_io) ) - 2*inner( x_i(p_io), Ji*dp ) - norm2( Ji*dp ) )
  //   = 1/N (xot xo - 2 xit Ji dp - norm2(Ji dp)) =
  //   = 1/N (xot xo - 2 xit Ji dp - dpt Jit Ji dp =
  //   = 1/N (xot xo - 2 xit Ji dp + xit Ji dp) =
  //   = 1/N (xot xo - xit Ji dp )
  //   = 1/N (xot xo + xit Ji inv(Jit Ji) Jit xi )
  //   = 1/N (xot xo + xot Jo inv(Jit Ji) Jot xo )
  //   = 1/N (xot (I + Jo inv(Jit Ji) Jot) xo )
  //
  // I already optimized the inlier-and-outlier problem, so I have
  //
  //   JtJ = sum(outer( ji, ji ) + outer( jo, jo ))
  //       = Jit Ji + Jot Jo
  //
  // and
  //
  //   inv(Jit Ji) = inv(JtJ - Jot Jo)
  //
  // Woodbury identity:
  //
  //   inv(Jit Ji) = inv(JtJ - Jot Jo) =
  //               = inv(JtJ) - inv(JtJ) Jot inv(-I + Jo inv(JtJ) Jot) Jo inv(JtJ)
  //
  // Let M = inv(JtJ) ->
  //   inv(Jit Ji) = M - M Jot inv(-I + Jo M Jot) Jo M
  //
  // So f = 1/N (xot (I + Jo M Jot - Jo M Jot inv(-I + Jo M Jot) Jo M Jot) xo )
  //
  // Let A = Jo M Jot ->
  //    f = 1/N (xot (I + A - A inv(-I + A) A) xo )
  //
  // If I'm looking at a single outlier measurement then A is a scalar and
  //
  //    f = 1/N xo^2 (A+1 - A^2/(A-1))
  //      = 1/N xo^2 ( - 1/(A-1)) =
  //      = xo^2 / ( N* (1-A)) =
  //      = xo^2 / ( N* (1 - jt inv(JtJ) j))
  //
  // I just solved the nonlinear optimization problem, so I already have
  // inv(JtJ). And for any one measurement, the outlier factor is
  //
  //    x^2 / (1 - jt inv(JtJ) j) / Nmeasurement
  //
  // where x is the scalar measurement, j is its vector gradient and JtJ = Jt
  // * J and J is a matrix of ALL the gradients of all the measurements
  //
  // This is the the "self" error difference + the "others" error difference.
  // If I split this up I get
  //
  //   f_others = f_both - x^2/N =
  //              x^2/Nmeasurement ( 1 / (1 - jt inv(JtJ) j) - 1 ) =
  //              x^2/Nmeasurement ( (1 - (1 - jt inv(JtJ) j)) / (1 - jt inv(JtJ) j) ) =
  //              x^2/Nmeasurement ( jt inv(JtJ) j / (1 - jt inv(JtJ) j) ) =

  if(measurementGroupSize <= 1)
    measurementGroupSize = 1;

  dogleg_computeJtJfactorization( point, ctx );
  bool result;
  if(ctx->is_sparse)
    result = getOutliernessFactors_sparse(factors, measurementGroupSize, Ngroups, point, ctx);
  else
    result = getOutliernessFactors_dense(factors, measurementGroupSize, Ngroups, point, ctx);

#if 0
  if( result )
  {
    int  Nstate        = ctx->Nstate;
    int  Nmeasurements = ctx->Nmeasurements;



    static FILE* fp = NULL;
    if(fp == NULL)
      fp = fopen("/tmp/check-outlierness.py", "w");
    static int count = -1;
    count++;
    if(count > 5)
      goto done;

    fprintf(fp, "# WARNING: all J here are unscaled with SCALE_....\n");

    if(count == 0)
    {
      fprintf(fp,
              "import numpy as np\n"
              "import numpysane as nps\n"
              "np.set_printoptions(linewidth=100000)\n"
              "\n");
    }

    fprintf(fp, "x%d = np.array((", count);
    for(int j=0;j<Nmeasurements;j++)
      fprintf(fp, "%.20g,", point->x[j]);
    fprintf(fp,"))\n");

    if( ctx->is_sparse )
    {
      fprintf(fp, "J%d = np.zeros((%d,%d))\n", count, Nmeasurements, Nstate);
      for(int imeas=0;imeas<Nmeasurements;imeas++)
      {
        for(int j = P(point->Jt, imeas);
            j     < (int)P(point->Jt, imeas+1);
            j++)
        {
          int irow = I(point->Jt, j);
          fprintf(fp, "J%d[%d,%d] = %g\n", count,
                  imeas, irow, X(point->Jt, j));
        }
      }
    }
    else
    {
      fprintf(fp, "J%d = np.array((\n", count);
      for(int j=0;j<Nmeasurements;j++)
      {
        fprintf(fp, "(");
        for(int i=0;i<Nstate;i++)
          fprintf(fp, "%.20g,", point->J_dense[j*Nstate+i]);
        fprintf(fp, "),\n");
      }
      fprintf(fp,"))\n");
    }

    fprintf(fp, "Nmeasurements = %d\n", Nmeasurements);
    fprintf(fp, "Ngroups = %d\n", Ngroups);
    fprintf(fp, "measurementGroupSize = %d\n", measurementGroupSize);

    fprintf(fp, "factors_got = np.array((");
    for(int j=0;j<Ngroups;j++)
      fprintf(fp, "%.20g,", factors[j]);
    fprintf(fp,"))\n");

    fprintf(fp,
            "pinvj = np.linalg.pinv(J%1$d)\n"
            "imeas0 = [measurementGroupSize*i for i in xrange(Ngroups)]\n"
            "jslices = [J%1$d[imeas0[i]:(imeas0[i]+measurementGroupSize), :] for i in xrange(Ngroups)]\n"
            "xslices = [x%1$d[imeas0[i]:(imeas0[i]+measurementGroupSize)   ] for i in xrange(Ngroups)]\n"
            "A       = [nps.matmult(jslices[i],pinvj[:,imeas0[i]:(imeas0[i]+measurementGroupSize)]) for i in xrange(Ngroups)]\n"
            "factors_ref = np.array([nps.inner(xslices[i], nps.matmult(np.eye(measurementGroupSize) + A[i] - nps.matmult(A[i], np.linalg.solve(A[i]-np.eye(measurementGroupSize),A[i])), nps.transpose(xslices[i])).ravel()) for i in xrange(Ngroups)]) / Nmeasurements\n",
            count);

    fprintf(fp, "print 'normdiff: {}'.format(np.linalg.norm(factors_ref-factors_got))\n");

    if(measurementGroupSize <= 1)
    {
      fprintf(fp,
              "factors_ref1 = x%1$d * x%1$d / (1.0 - nps.inner(J%1$d, nps.transpose(pinvj))) / Nmeasurements\n",
              count

              );
      fprintf(fp, "print 'normdiff1: {}'.format(np.linalg.norm(factors_ref1-factors_got))\n");
      fprintf(fp, "print 'normrefref1: {}'.format(np.linalg.norm(factors_ref1-factors_ref))\n");
    }
    fflush(fp);
  }
#endif

 done:
  return result;
}



#define OUTLIER_N_STDEVS_THRESHOLD        4
#define OUTLIER_CONFIDENCE_DROP_THRESHOLD 0.05
static void markPotentialOutliers(// output, input
                                  struct dogleg_outliers_t* markedOutliers,

                                  // input
                                  const double* factors,
                                  int Ngroups)
{
    // I have the outlier factors. How outliery is too outliery? Currently I
    // look at the distribution of the outlier factors across my dataset,
    // and focus on those measurements whose outlier factors are outliers.
    //
    // First, compute the mean, variance of the factors
    double outlierness_sum_sq_diff_mean = 0.0;
    double outlierness_sum              = 0.0;

    int N_in_statistics = 0;
    {
        // I do this in 2 passes. I like my floating-point precision
        for(int i=0; i<Ngroups; i++)
        {
            if(markedOutliers[i].marked)
                // this feature has already been designated an outlier
                continue;
            if( factors[i] == DBL_MAX )
            {
                // it's an outlier, but I don't want to include it in my
                // statistics, since it'll break them
                markedOutliers[i].markedPotential = true;
                continue;
            }
            markedOutliers[i].markedPotential = false;
            outlierness_sum += factors[i];
            N_in_statistics++;
        }
        double mean = outlierness_sum / (double)N_in_statistics;

        for(int i=0; i<Ngroups; i++)
        {
            if(markedOutliers[i].marked)
                // this feature has already been designated an outlier
                continue;
            if( factors[i] == DBL_MAX )
                // it's an outlier, but I don't want to include it in my
                // statistics, since it'll break them
                continue;

            double d = factors[i] - mean;
            outlierness_sum_sq_diff_mean += d*d;
        }
    }

    // Everything with an outlier factor at least X standard deviations
    // above the mean is considered a candidate outlier. Throwing each of
    // these out I recompute the mean and standard deviation to find more
    // points X standard deviations above the mean. I keep marking potential
    // outliers in this way, updating the variance, stdev with each pass. I
    // keep going until no more potential outliers are marked. I update the
    // variance and stdev with each pass because each outlier can have a
    // strong effect on these statistics
    bool markedAnyPotential;
    do
    {
        markedAnyPotential = false;

        double mean = outlierness_sum              /(double)N_in_statistics;
        double var  = outlierness_sum_sq_diff_mean /(double)N_in_statistics;
        SAY_IF_VERBOSE("have outlierness mean,var: %g,%g", mean, var);

        for(int i=0; i<Ngroups; i++)
        {
            if( markedOutliers[i].marked ||
                markedOutliers[i].markedPotential)
                continue;

            double d = factors[i] - mean;
            if( d < 0.0 || d*d < (double)(OUTLIER_N_STDEVS_THRESHOLD*OUTLIER_N_STDEVS_THRESHOLD) * var )
                continue;

            // Outlierness factor is above X standard deviations above the mean.
            // Potential outlier
            outlierness_sum -= factors[i];

            // what happens to a variance when we remove a measurement xo?
            //
            // VN  = sum( (xi-m)^2 ) =
            //     = sum( xi^2 - 2xi m + m^2 )
            //     = sum( xi^2 )- N m^2
            //
            // VN'= sum( (xi-m')^2 )
            //    = sum( xi^2 - 2xi m' + m'^2 )
            //    = sum( xi^2 ) - xo^2 - (N-1) (m')^2
            //    = VN + N m^2 - xo^2 - (N-1) (m')^2
            //    = VN - xo^2 - (S-xo)^2/(N-1) + S^2/N
            //    = VN - xo^2  + (S^2 N - S^2 - S^2 N + 2 S xo N - xo^2 N )/(N-1)/N
            //    = VN - xo^2  - (S^2 - 2 S xo N + xo^2 N )/(N-1)/N
            //    = VN - ( xo^2 N^2 - xo^2* N + S^2 - 2 S xo N + xo^2 N )/(N-1)/N
            //    = VN - ( xo^2 N^2 + S^2 - 2 S xo N )/(N-1)/N
            //    = VN - ( N xo - S )^2/(N-1)/N
            d = outlierness_sum - factors[i]*(double)(N_in_statistics-1);
            outlierness_sum_sq_diff_mean -= d*d/(double)((N_in_statistics-1)*N_in_statistics);

            N_in_statistics--;
            markedAnyPotential                = true;
            markedOutliers[i].markedPotential = true;

            mean = outlierness_sum              / (double)N_in_statistics;
            var  = outlierness_sum_sq_diff_mean / (double)N_in_statistics;
            SAY_IF_VERBOSE("New potential outlier group: %d. new mean,var: %g,%g",
                           i, mean, var);
        }
    } while(markedAnyPotential);
}

bool dogleg_markOutliers(// output, input
                         struct dogleg_outliers_t* markedOutliers,
                         // output only
                         int* Noutliers,

                         // input
                         double (getConfidence)(int i_group_exclude),

                         // if outliers are grouped into sets, the group size is
                         // stated here
                         int measurementGroupSize,
                         int Ngroups,

                         dogleg_operatingPoint_t* point,
                         dogleg_solverContext_t* ctx)
{
    if(measurementGroupSize <= 1)
      measurementGroupSize = 1;

    // What is an outlier? Suppose I just found an optimum. I define an
    // outlier as an observation that does two things to the problem if I
    // remove that observation:
    //
    // 1. The cost function would improve significantly. Things would
    //    clearly improve because the cost function contribution of the
    //    removed point itself would get removed, but ALSO because the
    //    parameters could fit to the remaining data better without that
    //    extra observation.
    //
    // 2. The confidence of the solution does not significantly decrease.
    //    One could imagine a set of data that define the problem poorly,
    //    and produce a low cost function value for some (overfit) set of
    //    parameters. And one can imagine an extra point being added that
    //    defines the problem and increases the confidence of the solution.
    //    This extra point would suppress the overfitting, so this extra
    //    point would increase the cost function value. Condition 1 above
    //    would make this extra point look like an outlier, and this
    //    condition is meant to detect this case and to classify this point
    //    as NOT an outlier
    bool markedAny = false;

    double* factors = malloc(Ngroups * sizeof(double));
    if(factors == NULL)
    {
        SAY("Error allocating factors");
        goto done;
    }

    if(!dogleg_getOutliernessFactors(factors, measurementGroupSize, Ngroups, point, ctx))
        goto done;

    markPotentialOutliers( markedOutliers,
                           factors, Ngroups);


    // OK then. I have my list of POTENTIAL outliers. These all have
    // suspicious outlierness factors. I check to see how much confidence I
    // would lose if I were to throw out any of these measurements, and
    // accept the outlier ONLY if the confidence loss is acceptable
    double confidence0 = getConfidence(-1);
    if( confidence0 < 0.0 )
        return false;

    SAY_IF_VERBOSE("Initial confidence: %g", confidence0);

    *Noutliers = 0;
    for(int i=0; i<Ngroups; i++)
    {
        if(markedOutliers[i].marked)
        {
          (*Noutliers)++;
          continue;
        }
        if(!markedOutliers[i].markedPotential)
            continue;

        // Looking at potential new outlier
        double confidence_excluded = getConfidence(i);
        if( confidence_excluded < 0.0 )
            return false;

        double confidence_drop_relative = 1.0 - confidence_excluded / confidence0;
        if( confidence_drop_relative < OUTLIER_CONFIDENCE_DROP_THRESHOLD )
        {
            // I would lose less than X of my confidence. OK. This is an
            // outlier. Throw it away.
            markedOutliers[i].marked = true;
            markedAny                = true;
            SAY_IF_VERBOSE("Excluding group %d produces a confidence: %g. relative loss: %g... YES an outlier; confidence drops little",
                           i, confidence_excluded, confidence_drop_relative);
            (*Noutliers)++;
        }
        else
        {
            SAY_IF_VERBOSE("Excluding group %d produces a confidence: %g. relative loss: %g... NOT an outlier: confidence drops too much",
                           i, confidence_excluded, confidence_drop_relative);
        }
    }

 done:
    free(factors);
    return markedAny;
}

// This function is just for debug reporting. It is probably too slow to
// call in general: it computes the confidence for each feature to see the
// confidence change if the feature were to be removed. Normally we do this
// ONLY for potential outliers
void dogleg_reportOutliers( double (getConfidence)(int i_group_exclude),

                            // if outliers are grouped into sets, the group size
                            // is stated here
                            int measurementGroupSize,
                            int Ngroups,

                            dogleg_operatingPoint_t* point,
                            dogleg_solverContext_t* ctx)
{
    if(measurementGroupSize <= 1)
      measurementGroupSize = 1;

    double* factors = malloc(Ngroups * sizeof(double));
    if(factors == NULL)
    {
        SAY("Error allocating factors");
        goto done;
    }

    dogleg_getOutliernessFactors(factors, measurementGroupSize, Ngroups, point, ctx);

    SAY("## Outlier statistics");
    SAY("# i_feature outlier_factor confidence_drop_relative_if_removed");

    double confidence_full = getConfidence(-1);

    for(int i=0; i<Ngroups; i++)
    {
      double confidence = getConfidence(i);
      double rot_confidence_drop_relative = 1.0 - confidence / confidence_full;

      SAY("%3d %9.3g %9.3g",
          i,
          factors[i],
          rot_confidence_drop_relative);
    }

 done:
    free(factors);
}
