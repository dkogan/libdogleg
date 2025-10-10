// -*- mode: C; c-basic-offset: 2 -*-

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "dogleg.h"

// This is a trivial sample application to demonstrate libdogleg in action.
// Let's say that I have a simple non-linear model
//
// a*b * x**2 + b*c * y**2 + c * x*y + d * x + e * y * f = measurements
//
// here I'm trying to estimate the vector (a,b,c,d,e,f) to most closely fit the
// data vector measurements. This problem is clearly non-sparse, but both sparse
// and dense versions of libdogleg are demonstrated here.
//
// First I generate some noise-corrupted data, and then use libdogleg to solve
// the problem.

// My state vector (a,b,c,d,e,f) has 6 elements
#define Nstate 6

// I simulate my measurements using these as the TRUE values for the model
#define REFERENCE_A 1.0
#define REFERENCE_B 2.0
#define REFERENCE_C 3.0
#define REFERENCE_D 4.0
#define REFERENCE_E 5.0
#define REFERENCE_F 6.0

// I simulate by sampling the x-y space in a grid. This grid is defined here
#define SIMULATION_GRID_WIDTH   10
#define SIMULATION_GRID_MIN    -10
#define SIMULATION_GRID_DELTA   2.0
#define Nmeasurements (SIMULATION_GRID_WIDTH*SIMULATION_GRID_WIDTH)


static double allx                [Nmeasurements];
static double ally                [Nmeasurements];
static double allm_simulated_noisy[Nmeasurements];

static void simulate(void)
{
  for(int i=0; i<Nmeasurements; i++)
  {
    double x = allx[i];
    double y = ally[i];

    allm_simulated_noisy[i] =
      REFERENCE_A*REFERENCE_B * x*x +
      REFERENCE_B*REFERENCE_C * y*y +
      REFERENCE_C * x*y +
      REFERENCE_D * x +
      REFERENCE_E * y +
      REFERENCE_F +
      ((double)random() / (double)RAND_MAX - 0.5) * 1.0; // +- 0.5 units of uniformly-random noise
  }
}

static void generateSimulationGrid(void)
{
  int i = 0;

  for(int ix=0; ix<SIMULATION_GRID_WIDTH; ix++)
  {
    double x = SIMULATION_GRID_MIN + ix*SIMULATION_GRID_DELTA;

    for(int iy=0; iy<SIMULATION_GRID_WIDTH; iy++)
    {
      double y = SIMULATION_GRID_MIN + iy*SIMULATION_GRID_DELTA;
      allx[i] = x;
      ally[i] = y;
      i++;
    }
  }
}

static void optimizerCallback(const double*   p,
                              double*         x,
                              cholmod_sparse* Jt,
                              void*           cookie __attribute__ ((unused)) )
{

  // These are convenient so that I only apply the casts once
  int*    Jrowptr = (int*)Jt->p;
  int*    Jcolidx = (int*)Jt->i;
  double* Jval    = (double*)Jt->x;

  int iJacobian = 0;
#define STORE_JACOBIAN(col, g)                  \
        do                                      \
        {                                       \
          Jcolidx[ iJacobian ] = col;           \
          Jval   [ iJacobian ] = g;             \
          iJacobian++;                          \
        } while(0)


  double norm2_x = 0.0;

  for(int i=0; i<Nmeasurements; i++)
  {
    x[i] =
      p[0] * p[1] * allx[i]*allx[i] +
      p[1] * p[2] * ally[i]*ally[i] +
      p[2] *        allx[i]*ally[i] +
      p[3] *        allx[i] +
      p[4] *        ally[i] +
      p[5]
      - allm_simulated_noisy[i];

    norm2_x += x[i]*x[i];

    // In this sample problem, every measurement depends on every element of the
    // state vector, so I loop through all the state vectors here. In practice
    // libdogleg is meant to be applied to sparse problems, where this internal
    // loop would be MUCH shorter than Nstate long
    Jrowptr[i] = iJacobian;
    STORE_JACOBIAN( 0, p[1]*allx[i]*allx[i] );
    STORE_JACOBIAN( 1, p[0]*allx[i]*allx[i] + p[2] * ally[i]*ally[i] );
    STORE_JACOBIAN( 2, p[1] * ally[i]*ally[i] + allx[i]*ally[i] );
    STORE_JACOBIAN( 3, allx[i] );
    STORE_JACOBIAN( 4, ally[i] );
    STORE_JACOBIAN( 5, 1.0  );
  }
  Jrowptr[Nmeasurements] = iJacobian;

#undef STORE_JACOBIAN
}

static void optimizerCallback_dense(const double*   p,
                                    double*         x,
                                    double*         J,
                                    void*           cookie __attribute__ ((unused)) )
{
  int iJacobian = 0;
#define STORE_JACOBIAN(col, g) J[ iJacobian++ ] = g


  double norm2_x = 0.0;

  for(int i=0; i<Nmeasurements; i++)
  {
    x[i] =
      p[0] * p[1] * allx[i]*allx[i] +
      p[1] * p[2] * ally[i]*ally[i] +
      p[2] *        allx[i]*ally[i] +
      p[3] *        allx[i] +
      p[4] *        ally[i] +
      p[5]
      - allm_simulated_noisy[i];

    norm2_x += x[i]*x[i];

    // In this sample problem, every measurement depends on every element of the
    // state vector, so I loop through all the state vectors here. In practice
    // libdogleg is meant to be applied to sparse problems, where this internal
    // loop would be MUCH shorter than Nstate long
    STORE_JACOBIAN( 0, p[1]*allx[i]*allx[i] );
    STORE_JACOBIAN( 1, p[0]*allx[i]*allx[i] + p[2] * ally[i]*ally[i] );
    STORE_JACOBIAN( 2, p[1] * ally[i]*ally[i] + allx[i]*ally[i] );
    STORE_JACOBIAN( 3, allx[i] );
    STORE_JACOBIAN( 4, ally[i] );
    STORE_JACOBIAN( 5, 1.0  );
  }

#undef STORE_JACOBIAN
}



#define GREEN       "\x1b[32m"
#define RED         "\x1b[31m"
#define COLOR_RESET "\x1b[0m"


int main(int argc, char* argv[] )
{
  const char* usage =
    "Usage: %s [--check] [--diag vnlog|human] [--test-gradients] sparse|dense\n";

  struct option opts[] = {
    { "diag",           required_argument, NULL, 'd' },
    { "test-gradients", no_argument,       NULL, 'g' },
    { "check",          no_argument,       NULL, 'c' },
    { "help",           no_argument,       NULL, 'h' },
    {}
  };

  bool is_sparse      = false;
  bool test_gradients = false;
  bool check          = false;
  int  debug          = 0;

  int opt;
  do
  {
    // "h" means -h does something
    opt = getopt_long(argc, argv, "+h", opts, NULL);
    switch(opt)
    {
    case -1:
      break;

    case 'h':
      printf(usage, argv[0]);
      return 0;

    case 'd':
      if(0 == strcmp("vnlog", optarg))
      {
        debug |= DOGLEG_DEBUG_VNLOG;
        break;
      }
      if(0 == strcmp("human", optarg))
      {
        debug |= 1;
        break;
      }
      fprintf(stderr, "--diag must be followed by 'vnlog' or 'human'\n");
      return 1;

    case 'g':
      test_gradients = true;
      break;

    case 'c':
      check = true;
      break;

    case '?':
      fprintf(stderr, "Unknown option\n\n");
      fprintf(stderr, usage, argv[0]);
      return 1;
    }
  } while( opt != -1 );

  const int Nargs_remaining = argc-optind;
  if( Nargs_remaining != 1 )
  {
    fprintf(stderr, "Need exactly 1 non-option argument: 'sparse' or 'dense'. Got %d\n\n",Nargs_remaining);
    fprintf(stderr, usage, argv[0]);
    return 1;
  }
  if( 0 == strcmp(argv[optind], "dense") )
    is_sparse = false;
  else if( 0 == strcmp(argv[optind], "sparse") )
    is_sparse = true;
  else
  {
    fprintf(stderr, "The final argument must be 'sparse' or 'dense'\n\n");
    fprintf(stderr, usage, argv[0]);
    return 1;
  }

  if(check && test_gradients)
  {
    fprintf(stderr, "--check and --test-gradients are exclusive\n");
    return 1;
  }

  if(!check)
  {
    if(is_sparse)
      fprintf(stderr, "Using SPARSE math\n");
    else
      fprintf(stderr, "Using DENSE math\n");
  }


  srandom( 0 ); // I want determinism here

  generateSimulationGrid();
  simulate();

  dogleg_parameters2_t dogleg_parameters;
  dogleg_getDefaultParameters(&dogleg_parameters);
  dogleg_parameters.dogleg_debug = debug;

  double p[Nstate];

  // I start solving with all my state variables set to some random noise
  for(int i=0; i<Nstate; i++)
    p[i] = ((double)random() / (double)RAND_MAX - 0.1) * 1.0; // +- 0.1 units of uniformly-random noise

  if(!check)
  {
    fprintf(stderr, "starting state:\n");
    for(int i=0; i<Nstate; i++)
      fprintf(stderr, "  p[%d] = %f\n", i, p[i]);
  }

  // This demo problem is dense, so every measurement depends on every state
  // variable. Thus ever element of the jacobian is non-zero
  int Jnnz = Nmeasurements * Nstate;

  // first, let's test our gradients. This is just a verification step to make
  // sure the optimizerCallback() is written correctly. Normally, you would do
  // this as a check when developing your program, but would turn this off in
  // the final application. This will generate LOTS of output. You need to make
  // sure that the reported and observed gradients match (the relative error is
  // low)
  if(!check)
    fprintf(stderr, "have %d variables\n", Nstate);
  if( test_gradients )
  {
    for(int i=0; i<Nstate; i++)
    {
      fprintf(stderr, "checking gradients for variable %d\n", i);
      if( is_sparse )
        dogleg_testGradient(i, p, Nstate, Nmeasurements, Jnnz, &optimizerCallback, NULL);
      else
        dogleg_testGradient_dense(i, p, Nstate, Nmeasurements, &optimizerCallback_dense, NULL);
    }
    return 0;
  }

  if(!check)
    fprintf(stderr, "SOLVING:\n");

  double optimum;
  if( is_sparse )
    optimum = dogleg_optimize2(p, Nstate, Nmeasurements, Jnnz,
                               &optimizerCallback, NULL,
                               &dogleg_parameters, NULL);
  else
    optimum = dogleg_optimize_dense2(p, Nstate, Nmeasurements,
                                     &optimizerCallback_dense, NULL,
                                     &dogleg_parameters, NULL);

  if(check)
  {
    if(optimum < 0)
    {
      printf(RED "ERROR: the optimization did not converge\n" COLOR_RESET);
      return 1;
    }

    printf(GREEN "OK: the optimization converged to an optimum  of norm2(x)=%.1f\n" COLOR_RESET,
           optimum);

    bool anyfailed = false;
    const double pref[] =
      { REFERENCE_A,
        REFERENCE_B,
        REFERENCE_C,
        REFERENCE_D,
        REFERENCE_E,
        REFERENCE_F };
    for(int i=0; i<Nstate; i++)
    {
      const double err = p[i] - pref[i];
      if(fabs(err) < 5e-2)
        printf(GREEN "OK: parameter %d recovered: psolved=%.3f pref=%.3f perr=%.3f\n" COLOR_RESET,
               i, pref[i], p[i], err);
      else
      {
        printf(RED "ERROR: parameter %d was NOT recovered: psolved=%.3f pref=%.3f perr=%.3f\n" COLOR_RESET,
               i, pref[i], p[i], err);
        anyfailed = true;
      }
    }

    return anyfailed ? 1 : 0;
  }
  else
  {
    fprintf(stderr, "Done. Optimum = %f\n", optimum);
    if(optimum < 0)
      fprintf(stderr, "optimum<0: an error has occurred\n");
    fprintf(stderr, "optimal state:\n");
    for(int i=0; i<Nstate; i++)
      fprintf(stderr, "  p[%d] = %f\n", i, p[i]);
  }

  return 0;
}
