#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "hdf5.h"
#include "hdf5_hl.h"

#include "mcmc.h"

int mcmc_save_trace(mcmc_configuration config, int param_num, 
		    const char* param_name)
{
    double* trace;
    trace = mcmc_trace(config, param_num);

    hsize_t dims[1] = {config.n_iter};

    H5LTmake_dataset (config.file_id, param_name, 1,
		      dims, H5T_NATIVE_DOUBLE, trace);
    return 0;
}

int mcmc_set_file(mcmc_configuration* config, const char* file_name)
{
    hid_t file_id;
    
    file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    config->file_id = file_id;    
    return 0;
}

void mcmc_error(const char* message)
{
    fputs(message, stderr);
    exit(1);
}

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

void mcmc_set_joint_prior(mcmc_configuration* config,
                          double(*joint_prior)(double*))
{
    double (*p_joint_prior)(double* params);
    p_joint_prior = joint_prior;
    config->joint_prior = p_joint_prior;
}

void mcmc_set_data_probability(mcmc_configuration* config,
                               double(*data_probability)(double*, double*))
{
    double (*p_data_probability)(double* data, double* params);
    p_data_probability = data_probability;
    config->data_probability = p_data_probability;
}

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

double mcmc_normal(mcmc_configuration config, double mu, int param_num)
{
    double* p_sigma = (double*)config.proposal_distribution_args[param_num];
    double sigma = *p_sigma;
    gsl_rng *r = config.gslrng;
    double x;
    x = gsl_ran_gaussian (r, sigma);
    return mu + x;
}

double* mcmc_trace(mcmc_configuration config, int param_num)
{
    int n_iter = config.n_iter;
    int n_param = config.n_param;
    
    double* trace = malloc(n_iter*sizeof(double));    
    if (trace == NULL)
	{
	    mcmc_error("mcmc: error allocating trace\n");
	}
    int i;
    for (i=0; i<n_iter; i++)
	{
	    trace[i] = config.results[param_num+i*n_param];
	}
    return trace;
}

void mcmc_free(mcmc_configuration* config)
{
    
}

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

            current_posterior = (config.joint_prior(params)
				 *config.data_probability(data, params));

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
                }
            offset = i*n_param;
            memcpy(results+offset, params, n_param*sizeof(double));
        }
    gsl_rng_free (r);
    return 0;
}
