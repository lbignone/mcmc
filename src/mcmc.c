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
mcmc_configuration mcmc_initialize (int n_param, int n_iter, int n_times)
{
    mcmc_configuration config;

    config.n_iter = n_iter;
    config.n_times = n_times;
       
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
    config.probability = malloc(n_iter*sizeof(double));
    config.accepted = malloc(n_iter*sizeof(int));
    config.accepted[0] = 0;

    config.log_probability = 0;

    config.proposed = malloc(n_param*n_iter*sizeof(double));

    mcmc_initialize_rng(&config);

    config.save_function = &mcmc_save_traces;

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
		break;
            }
	case FIX:
	    {
		p_proposal_distribution = &mcmc_fix;
		config->proposal_distributions[param_num]
		    = p_proposal_distribution;
		break;
	    }
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
double mcmc_run (mcmc_configuration config, double* data)
{
    int n_iter = config.n_iter;
    int n_times = config.n_times;
    int n_param = config.n_param;
    double* results = config.results;
    double* proposed = config.proposed;

    double params[n_param];
    double proposed_params[n_param];

    double metropolis_ratio;
    double proposed_posterior;
    double current_posterior;
    double data_prob;
    double prior;
   
    memcpy(params, config.parameters, n_param*sizeof(double));

    data_prob = config.data_probability(data, params);
    prior = config.joint_prior(params);

    if(config.log_probability)
    {
        config.current_posterior = log(prior) + data_prob;
    }
    else
    {
        config.current_posterior = (prior*data_prob);
    }

    config.current_posterior = (config.joint_prior(params)
				*data_prob);

    int ACCEPTED;
    double u;
    int i, j, k, offset;
    int accepted_number = 0;
    for (k=0; k<n_times; k++)
    {
	for (i=0; i<n_iter; i++)
        {
            ACCEPTED = 0;
            for (j=0; j<n_param; j++)
	    {
		proposed_params[j] =
		    (*config.proposal_distributions[j])
		    (config, params[j], j);
	    }

	    data_prob = config.data_probability(data, proposed_params);
            prior = config.joint_prior(proposed_params);
            if (prior != 0.0)
            {    
                current_posterior = config.current_posterior;
                
                if(config.log_probability)
                {
                    proposed_posterior = log(prior) + data_prob;
                    metropolis_ratio = exp(proposed_posterior - current_posterior);
                }
                else
                {
                    proposed_posterior = (prior*data_prob);
                    metropolis_ratio = proposed_posterior/current_posterior;
                }

                config.accepted[i] = 0;
                if (metropolis_ratio >= 1)
                {
                    ACCEPTED = 1;
                }
                else
                {
                    u = gsl_rng_uniform (config.gslrng);
                    if (u <= metropolis_ratio)
                    {
                        ACCEPTED = 1;
                    }
                }
            }
            else
            {
                data_prob = 0.0;
            }

            if (ACCEPTED)
	    {
		memcpy(params, proposed_params, n_param*sizeof(double));
		config.current_posterior = proposed_posterior;
		accepted_number++;
		config.accepted[i] = 1;
	    }
            offset = i*n_param;x
            memcpy(results+offset, params, n_param*sizeof(double));
	    config.probability[i] = data_prob;
	    memcpy(proposed+offset, proposed_params, n_param*sizeof(double));
        }
	config.save_function(config);
    }
    double acceptance_rate = accepted_number/n_iter;
    return acceptance_rate;
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


/* Continue from existing file
 *
 * Recover the last state of a previously generated file and continue
 * from there. 
 *
 * Only recover the last parameter state, not the proposal arguments
 *
 * Input:
 * - config:    pointer to a mcmc_configuration struct
 * - file_name: name of file to be open
 * 
 * Return: 0
 */
int mcmc_continue_from_file(mcmc_configuration* config, const char* file_name)
{
    int i;
    hid_t file_id;

    herr_t err;
    hsize_t nfields, nrecords;

    file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0)
	mcmc_error("mcmc: Could not open file\n");

    err =  H5TBget_table_info (file_id, "Traces", &nfields, &nrecords);
    if (err < 0)
	mcmc_error("mcmc: Could not open Traces table\n");

    size_t *field_offset = malloc(nfields*sizeof(size_t));
    size_t *dst_size = malloc(nfields*sizeof(size_t));
    for (i=0; i < nfields; ++i)
    {
	field_offset[i] = i*sizeof(double);
	dst_size[i] = sizeof(double);
    }

    size_t type_size = nfields*sizeof(double);

    double *data = malloc(nfields*sizeof(double));

    err = H5TBread_records(file_id, "Traces", nrecords-1, 1, type_size,
    			   field_offset, dst_size, data);

    if (err < 0)
	mcmc_error("mcmc: Could not read Traces table\n");

    if (nfields != config->n_param)
	mcmc_error("mcmc: Incompatible file\n");
	
    config->file_id = file_id;
    for(i=0; i<nfields; i++)
    	config->parameters[i] = data[i];
}

int mcmc_set_meta(mcmc_configuration* config)
{
    hid_t file_id = config->file_id;

    int n_meta = 1;
    hid_t field_types[n_meta];
    size_t dst_offsets[n_meta];
    size_t dst_sizes[n_meta];

    size_t dst_size = n_meta*sizeof(double);
    field_types[0] = H5T_NATIVE_DOUBLE;	    
    dst_sizes[0] = sizeof(double);
    dst_offsets[0] = 0;
    
    const char *field_names[1] = {"probability"};

    hsize_t chunk_size = 10;
    int *fill_data = NULL;
    int compress = 0;

    herr_t status;
    status = H5TBmake_table("Meta", file_id, "Probability", n_meta, 0, 
			    dst_size, field_names, dst_offsets, field_types,
			    chunk_size, fill_data, compress, NULL);

    dst_size = n_meta*sizeof(int);
    field_types[0] = H5T_NATIVE_INT;	    
    dst_sizes[0] = sizeof(int);
    dst_offsets[0] = 0;
    
    const char *field_names_2[1] = {"accepted"};

    status = H5TBmake_table("Meta", file_id, "Accepted", n_meta, 0, 
			    dst_size, field_names_2, dst_offsets, field_types,
			    chunk_size, fill_data, compress, config->accepted);

    

}

/* Free all memory associated with a mcmc_configuration
 *
 *
*/
void mcmc_free(mcmc_configuration* config)
{
    free(config->proposal_distributions);
    free(config->proposal_distribution_args);
    free(config->results);
    
    H5Fclose (config->file_id);
    free(config);
}
