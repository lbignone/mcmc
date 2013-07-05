#ifndef MCMC_H
#define MCMC_H

#include <stdlib.h>

typedef struct mcmc_configuration;

double* mcmc_allocate_results(mcmc_configuration config);

int mcmc_run(double* data, mcmc_configuration config, double* results);

#endif MCMC_H
