#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_rng.h>

#include "mcmc.h"

struct mcmc_configuration
{
  int n_iter;

  int n_param;
  double* parameters;

  void* joint_prior;
  void* proposal_distributions;
  void* data_probability;

  unsigned long int seed;
};

double* mcmc_allocate_results(mcmc_configuration config);
{
  int n_iter = config.n_iter;
  int n_param = config.n_param;
  
  double* results = malloc(n_iter*n_param*sizeof(double));
  if (results == NULL)
    {
      fprintf(stderr, "mcmc: error allocating memory");
      exit(1);
    }

  return results;
}
  
			   

int mcmc_run(double* data, int n_data, mcmc_configuration config, double* results)
{  
  int n_iter = config.n_iter;
  int n_param = config.n_param;

  double* params[n_param];
  double* propose_params[n_param];

  double metroplolis_ratio;
  double propose_posterior;
  double current_posterior;

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (j=0; j<n_param; j++)
    {
      results[0+j*n_iter] = params[j];
    }
  
  int ACCEPTED;
  double u;
  for (i=1; i<n_iter; i++)
    {
      ACCEPTED = 0;
      for (j=0; j<n_param; j++)
	{
	  propose_params[j] = config.proposal_distributions[j](params[j]);
	}

      propose_posterior = config.prior(propose_params)*config.data_probability(data, propose_params);
      current_posterior = config.prior(params)*config.data_probability(data, params);

      metropolis_ratio = propose_posterior/current_posterior;
      
      if (metropolis_ratio >= 1)
	{
	  ACCEPTED = 1;
	}
      else
	{
	  u = gsl_rng_uniform (r);
	  if (u <= metropolis_ratio)
	    {
	      ACCEPTED = 1;
	    }
	}
      for (j=0; j<n_param; j++)
	{
	  if (ACCEPTED)
	    {
	      results[i+j*n_param] = propose_params[j];
	    }
	  else
	    {
	      results[i+j*n_param] = params[j];
	    }
	}
    }
  gsl_rng_free (r);
}
