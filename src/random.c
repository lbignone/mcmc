#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>

#include "mcmc.h"

unsigned long int mcmc_get_seed()
{
    unsigned long int seed;
    int err;

    /* open /dev/random */
    FILE* file = fopen("/dev/random", "rb");
    if (file == NULL)
	mcmc_error("mcmc: error opening /dev/random\n");

    /* get random seed */
    err = fread(&seed, sizeof(seed), 1, file);
    if (err != 1)
	mcmc_error("mcmc: error reading /dev/random");

    fclose(file); 

    return seed;

}

void mcmc_set_seed(mcmc_configuration* config, unsigned long int seed)
{
    gsl_rng_set (config->gslrng, seed);
    config->seed = seed;
}

void mcmc_initialize_rng(mcmc_configuration* config)
{
    const gsl_rng_type * T;
    unsigned long int seed;
    
    /* choose default random number generator */
    T = gsl_rng_default; 

    /* allocate random number generator */
    config->gslrng = gsl_rng_alloc (T);
    if (config->gslrng == NULL)
	mcmc_error("mcmc: could not allocate random number generator\n");
    
    /* get a seed from /dev/random and save it  */
    seed = mcmc_get_seed(); 
    mcmc_set_seed(config, seed);
}
