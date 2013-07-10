#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "mcmc.h"

/* Gaussian shape proposal.
 *
 * Input:
 * - config:    mcmc_configuration
 * - mu:        parameter value
 * - param_num: parameter number
 *
 * Returns: gaussian deviate
 */
double mcmc_normal(mcmc_configuration config, double mu, int param_num)
{
    double* p_sigma = (double*)config.proposal_distribution_args[param_num];
    double sigma = *p_sigma;
    gsl_rng *r = config.gslrng;
    double x;
    x = gsl_ran_gaussian (r, sigma);
    return mu + x;
}
