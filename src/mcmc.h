#ifndef MCMC_H
#define MCMC_H

#include <gsl/gsl_rng.h>

#include "hdf5.h"

typedef struct mcmc_configuration mcmc_configuration;
typedef enum mcmc_proposal_type mcmc_proposal_type;

struct mcmc_configuration
{
    int n_iter;
    int n_times;

    int n_param;
    double* parameters;

    double current_posterior;  

    double (*joint_prior)(double* params);
    double (**proposal_distributions)
    (mcmc_configuration config, double param, int param_num);
    double (*data_probability)(double* data, double* params);
    void** proposal_distribution_args;

    int (*save_function)(mcmc_configuration config);

    double* results;

    hid_t file_id;

    unsigned long int seed;

    gsl_rng * gslrng;

    double *probability;
    double initial_data_prob;
    double *proposed;

    int *accepted;
    
};

enum mcmc_proposal_type
{
    NORMAL,
    FIX
};


void mcmc_error(const char* message);
mcmc_configuration mcmc_initialize (int n_param, int n_iter, int n_times);
double* mcmc_allocate_results(mcmc_configuration config);
double mcmc_normal(mcmc_configuration config, double mu, int param_num);
double mcmc_fix(mcmc_configuration config, double mu, int param_num);
double* mcmc_trace(mcmc_configuration config, int param_num);
double mcmc_run(mcmc_configuration config, double* data);

int mcmc_set_table(mcmc_configuration config, const char* field_names);
int mcmc_save_traces(mcmc_configuration config);

double mcmc_uniform(double x, double x_min, double x_max);
double mcmc_jeffreys(double x, double x_min, double x_max);

#endif // MCMC_H
