#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include "hdf5.h"

#include "mcmc.h"

#define NROWS 64
#define NCOLUMNS 2

/* Define uniform prior*/
double uniform(double x, double x_min, double x_max)
{
    return 1.0/(x_max - x_min);
}

/* Define Jeffreys's prior */
double jeffreys(double x, double x_min, double x_max)
{
    return 1.0/(x*log(x_max/x_min));
}

/* Compute joint prior for T and nu */
double joint_prior(double* params)
{
    double nu = params[0];
    double T = params[1];

    double T_min = 0.1;
    double T_max = 100;

    double nu_min = 1.0;
    double nu_max = 44.0;

    int OUT_OF_RANGE = 0;
    if (nu < nu_min || nu > nu_max)
	{
	    OUT_OF_RANGE = 1;
	}
    else if (T < T_min || T > T_max)
	{
	    OUT_OF_RANGE = 1;
	}

    if (OUT_OF_RANGE)
	{
	    return 0.0;
	}

    return uniform(nu, nu_min, nu_max) * jeffreys(T, T_min, T_max);
}

/* Compute data probability given the values of T and nu*/
double data_probability(double* data, double* params)
{

    double n_rows = NROWS;
    
    double nu = params[0];
    double T = params[1];
    double partial, result;
    double nui, Ti, fi;

    double sigma = 1.0;
    
    double sigma_L = 2.0;
    int i;
    partial = 0;
    for (i=0; i<NROWS; i++ )
	{
	    nui = data[0+i*NCOLUMNS];
	    Ti = data[1+i*NCOLUMNS];
	    fi = exp( -pow(nui - nu, 2.0) / (2*pow(sigma_L, 2.0)) );
	    partial += pow(Ti - T*fi, 2.0);
	}
    result = 
	pow(2.0*M_PI, -n_rows/2.0) * pow(sigma, -n_rows)
	* exp(-partial/(2.0*pow(sigma, 2.0)));
    return result;
}

/* Read data from file*/
double* read_data(int* n_data)
{
    char ignore[100];
    FILE* file = fopen("data/spectral.dat", "r");
    if (file == NULL)
	{
	    mcmc_error("Error: could not open file\n");
	}
    double* data;

    int n_rows = NROWS;
    int n_columns = NCOLUMNS;
    *n_data = n_columns*n_rows;

    data = malloc(*n_data*sizeof(double));
    if (data == NULL)
	{
	    mcmc_error("Error: could not allocate data");
	}
    
    fgets(ignore, 100, file);
    int i;
    for (i=0; i<n_rows; i++)
	{
	    fscanf(file, "%lf %lf", &data[0+i*n_columns], &data[1+i*n_columns]);
	}
    return data;
}


int main()
{
    struct mcmc_configuration config;
    double* results;

    /* Initialize for 2 parameters and 1e5 iterations*/
    config = mcmc_initialize(2, 50, 2000);

    /* Assign starting values */
    config.parameters[0] = 30.0; 
    config.parameters[1] = 5.0;

    /* Set gaussian proposal distributions with sigma = 1 for both
       parameters */
    mcmc_set_proposal_distribution(&config, 0, NORMAL, 1.0);
    mcmc_set_proposal_distribution(&config, 1, NORMAL, 1.0);

    /* set the joint prior */
    mcmc_set_joint_prior(&config, &joint_prior);

    /* set the data probability */
    mcmc_set_data_probability(&config, &data_probability);

    /* load data */
    int n_data;
    double* data = read_data(&n_data);


    /* set output file */
    mcmc_set_file(&config, "lines.h5");

    mcmc_set_meta(&config);

    /* set table */
    const char *field_names[2] = {"nu", "T"};
    mcmc_set_table(config, field_names);
    
    /* start mcmc loop */
    mcmc_run(config, data);
    
    /* save traces to disk */
    mcmc_save_trace(config, 0, "nu");
    mcmc_save_trace(config, 1, "T");

    
}
