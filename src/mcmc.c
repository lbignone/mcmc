#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include <gsl/gsl_rng.h>

#include "hdf5.h"
#include "hdf5_hl.h"

#include "mcmc.h"

/* Initialize the mcmc procedure. Allocate necessary memory 
 * 
 * Input:
 * - n_param: number of parameters to fit
 * - n_iter:  number of iterations to run
 *
 * Returns: a mcmc_configuration structure
 */
mcmc_configuration mcmc_initialize (int n_param, int n_iter)
{
    mcmc_configuration config;

    config.n_iter = n_iter;
    config.n_param = n_param;

    config.parameters = calloc(n_param, sizeof(double));
    config.proposal_distributions = malloc(n_param*sizeof(double(*)(double, void*)));
    config.proposal_distribution_args = malloc(n_param*sizeof(void**));
    int FAILED = 0;
    if (config.parameters == NULL || config.proposal_distributions == NULL)
        {
            FAILED = 1;
        }
    if (config.proposal_distribution_args == NULL)
        {
            FAILED = 1;
        }

    if (FAILED)
        {
            mcmc_error("mcmc: could not allocate parameters\n");
        }

    config.results = mcmc_allocate_results(config);

    mcmc_initialize_rng(&config);

    return config;
}

/* Set proposal distribution
 *
 * For a Gaussian proposal the proposal_type should be set to NORMAL
 * and the standard deviation should be given as the only extra argument
 *
 * Input:
 * - config:        pointer to a mcmc_configuration struct
 * - param_number:  parameter number for identification
 * - proposal_type: type of pdf to use as proposal
 * - ...:           extra arguments to be pass to the proposal pdf
 */
void mcmc_set_proposal_distribution(mcmc_configuration *config,
                                    int param_num,
                                    mcmc_proposal_type proposal_type,
                                    ...)
{
    va_list arguments;
    va_start(arguments, proposal_type);

    double (*p_proposal_distribution)
        (mcmc_configuration config, double param, int param_num);

    switch(proposal_type)
        {
        case NORMAL:
            {
                double sigma = va_arg(arguments, double);
                p_proposal_distribution = &mcmc_normal;
                config->proposal_distributions[param_num]
                    = p_proposal_distribution;
                config->proposal_distribution_args[param_num]
                    = malloc(sizeof(double));
                *(double*)config->proposal_distribution_args[param_num] = sigma;
                va_end ( arguments );
            }
            break;
        }
}

/*  Set the joint prior
 *
 * The function must follow the prototype double(double*) and expect
 * an array with the parameter values as input.
 *
 * Input:
 *
 * - config:      pointer to a mcmc_configuration struct
 * - joint_prior: pointer to a function with prototype
 *                double()(double*)
*/
void mcmc_set_joint_prior(mcmc_configuration* config,
                          double(*joint_prior)(double*))
{
    double (*p_joint_prior)(double* params);
    p_joint_prior = joint_prior;
    config->joint_prior = p_joint_prior;
}

/*  Set the data_probability
 *
 * The function must follow the prototype double(double*, double*) and
 * expect an array with the data as first argument and an array with
 * the parameter values as second argument.
 *
 * Input:
 *
 * - config:           pointer to a mcmc_configuration struct
 * - data_probability: pointer to a function with prototype
 *                     double()(double*, double*)
 */
void mcmc_set_data_probability(mcmc_configuration* config,
                               double(*data_probability)(double*, double*))
{
    double (*p_data_probability)(double* data, double* params);
    p_data_probability = data_probability;
    config->data_probability = p_data_probability;
}

/* Set a random seed
 *
 * Input: 
 * - config: pointer to a mcmc_configuration struct
 * - seed:   random seed
 */
void mcmc_set_seed(mcmc_configuration* config, unsigned long int seed)
{
    gsl_rng_set (config->gslrng, seed);
    config->seed = seed;
}

/* Start mcmc loop
 *
 * Input:
 * - config: mcmc_configuration struct
 * - data:   array containing the data
 * - n_data: not really used -> TO BE REMOVED
 * 
 * Return: 0
*/
int mcmc_run (mcmc_configuration config, double* data, int n_data)
{
    int n_iter = config.n_iter;
    int n_param = config.n_param;
    double* results = config.results;

    double params[n_param];
    double proposed_params[n_param];

    double metropolis_ratio;
    double proposed_posterior;
    double current_posterior;

    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    memcpy(results, config.parameters, n_param*sizeof(double));
    memcpy(params, config.parameters, n_param*sizeof(double));

    config.current_posterior = (config.joint_prior(params)
				*config.data_probability(data, params));


    int ACCEPTED;
    double u;
    int i, j, offset;
    for (i=1; i<n_iter; i++)
        {
            ACCEPTED = 0;
            for (j=0; j<n_param; j++)
                {
                    proposed_params[j] =
                        (*config.proposal_distributions[j])
                        (config, params[j], j);
                }

            proposed_posterior = (config.joint_prior(proposed_params)
				 *config.data_probability(data, proposed_params));

            current_posterior = config.current_posterior;

            metropolis_ratio = proposed_posterior/current_posterior;

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
            if (ACCEPTED)
                {
                    memcpy(params, proposed_params, n_param*sizeof(double));
		    config.current_posterior = proposed_posterior;
                }
            offset = i*n_param;
            memcpy(results+offset, params, n_param*sizeof(double));
        }
    gsl_rng_free (r);
    return 0;
}

/* Set output file
 *
 * Input:
 * - config:    pointer to a mcmc_configuration struct
 * - file_name: name of file to be used as output
 * 
 * Return: 0
 */
int mcmc_set_file(mcmc_configuration* config, const char* file_name)
{
    hid_t file_id;
    
    file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    config->file_id = file_id;    
    return 0;
}

/* Free all memory associated with a mcmc_configuration
 *
 *
*/
void mcmc_free(mcmc_configuration* config)
{
    /* TO DO*/
}
