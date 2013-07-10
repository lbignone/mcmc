#include <stdlib.h>

#include "hdf5.h"
#include "hdf5_hl.h"

#include "mcmc.h"

/* Extract the trace of a parameter from the results.
 * 
 * Input:
 * - config:     mcmc_configuration
 * - param_num:  parameter number
 *
 * Returns: trace array. 
 */
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

/* Save the trace to disk.
 *
 * The file name is set by config.file_id so a call to mcmc_set_file()
 * should be made first.
 *
 * Input:
 * - config:     mcmc_confguration
 * - param_num:  parameter number
 * - param_name: parameter name to be given to the dataset
 *
 * Returns: 0 
*/
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
