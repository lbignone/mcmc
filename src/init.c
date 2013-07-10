#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>

#include "mcmc.h"

/* Print error message to stderr and exit
 *
 * Input:
 * - message: error message
 */
void mcmc_error(const char* message)
{
    fputs(message, stderr);
    exit(1);
}

/* Allocate memory for results*
 *
 * Input config: THIS SHOULD CHANGE TO A POINTER (MAYBE)
 *
 * Return: pointer to newly allocated array
*/
double* mcmc_allocate_results (mcmc_configuration config)
{
    int n_iter = config.n_iter;
    int n_param = config.n_param;

    double* results = malloc(n_iter*n_param*sizeof(double));
    if (results == NULL)
        {
            mcmc_error("mcmc: error allocating memory for results");
        }

    return results;

}

/* Get a random seed from /dev/random
 *
 * Return: random seed
 */
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

/* Initialize random number generator
 *
 * input:
 * - config: pointer to mcmc_configuration struct
*/
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
