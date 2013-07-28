#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <mpi.h>

#include "mcmc.h"

const gsl_rng_type * T_c;
gsl_rng * r_c;

/* total number of data point to generate */
int n = 100; 

/* Parameters of the target distribution */
double sigma = 1.0;
double mu = 5.0;

/* Uniform prior*/
double uniform(double x, double x_min, double x_max)
{
    return 1.0/(x_max - x_min);
}

double joint_prior(double* params)
{
    double mu_min = -10.0;
    double mu_max = 10.0;

    double mu_p = params[0];

    return uniform(mu_p, mu_min, mu_max);
}

/* 
   Generate n data points with the target distribution and compute the probability of obtaining those same data points with a distribution given by params 

   Only the root process computes the final result;
*/
double data_probability(double* data, double* params)
{
    double x;
    double y;

    double prop_mu = params[0];

    int rank, numtasks;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    int n_per_proc = (int) floor(n/numtasks);

    MPI_Bcast(&prop_mu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double result = 1.0;
    double final_result = 1.0;
    
    int i;
    for (i=0; i<n_per_proc; i++)
    {
	x = gsl_ran_gaussian(r_c, sigma);
	x = x + mu;
	y = x - prop_mu;
	result *= gsl_ran_gaussian_pdf(y, sigma);
    }
    MPI_Reduce(&result, &final_result, 1, MPI_DOUBLE, MPI_PROD, 
	       0, MPI_COMM_WORLD);
   
    MPI_Barrier(MPI_COMM_WORLD);
   
    return final_result;
}

int main ()
{
    
    MPI_Init(NULL, NULL);

    gsl_rng_env_setup();
    T_c = gsl_rng_default;
    r_c = gsl_rng_alloc (T_c);

    int i, rank, numtasks, n_data;
    
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    
    double *data, *results;
    struct mcmc_configuration config;

    config = mcmc_initialize(1, 1e7);

    config.parameters[0] = 0.0; 
    
    mcmc_set_proposal_distribution(&config, 0, NORMAL, 1.0);
    
    mcmc_set_joint_prior(&config, &joint_prior);
    
    mcmc_set_data_probability(&config, &data_probability);

    mcmc_run(config, data, n_data);

    /* Only processor 0 writes the results */
    if (rank == 0)
    {
	mcmc_set_file(&config, "parallel.h5");
	mcmc_save_trace(config, 0, "mu");
    }
    
    MPI_Finalize();
    return 0;
}
