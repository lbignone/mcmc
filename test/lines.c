#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_rng.h>

#include "mcmc.h"

double uniform(double x, double x_min, double x_max)
{
  return 1.0/(x_max - x_min);
}

double jeffreys(double x, double x_min, double x_max)
{
  return 1.0/(x*log(x_max/x_min));
}

double joint_prior(double nu, double T)
{
  double T_min = 0.1;
  double T_max = 100;

  double nu_min = 1;
  double nu_max = 44;

  return uniform(nu, nu_min, nu_max) * jeffreys(T, T_min, T_max);
}

double data_probability(double* data, double* params)
{
  double f[
}

mcmc_configuration config;

config.n_iter = 500;
config.n_param = 2;

config.parameters = malloc(2*sizeof(double));
config.parameters[0] = 5.0;
config.parameters[1] = 30.0;





int mcmc_run(double* data, mcmc_configuration config, double* results)
