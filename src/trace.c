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

int mcmc_set_table(mcmc_configuration config, const char* field_names[])
{
    int n_param = config.n_param;
    hid_t field_types[n_param];

    size_t dst_offsets[n_param];
    size_t dst_sizes[n_param];

    size_t dst_size = n_param*sizeof(double);

    int i;
    for (i=0; i<n_param; i++)
    {
	field_types[i] = H5T_NATIVE_DOUBLE;	    
	dst_sizes[i] = sizeof(double);
	dst_offsets[i] = i*sizeof(double);
    }

    hsize_t chunk_size = 10;
    int *fill_data = NULL;
    int compress = 0;

    
	

    herr_t status;
    status = H5TBmake_table("Traces", config.file_id, "Traces", n_param, 1, 
			    dst_size, field_names, dst_offsets, field_types,
			    chunk_size, fill_data, compress, config.parameters);

    status = H5TBmake_table("Proposed", config.file_id, "Proposed", n_param, 0, 
			    dst_size, field_names, dst_offsets, field_types,
			    chunk_size, fill_data, compress, config.parameters);

    
}

int mcmc_save_traces(mcmc_configuration config)
{
    int n_param = config.n_param;
    int n_iter = config.n_iter;
    
    size_t dst_offsets[n_param];
    size_t dst_sizes[n_param];

    size_t dst_size = n_param*sizeof(double);

    int i;
    for (i=0; i<n_param; i++)
    {
	dst_sizes[i] = sizeof(double);
	dst_offsets[i] = i*sizeof(double);
    }

    herr_t status;
    status = H5TBappend_records(config.file_id, "Traces", n_iter,
				dst_size, dst_offsets, dst_sizes, 
				config.results);

    status = H5TBappend_records(config.file_id, "Proposed", n_iter,
				dst_size, dst_offsets, dst_sizes, 
				config.proposed);

    size_t dst_offsets_meta[1];
    size_t dst_sizes_meta[1];
    dst_offsets_meta[0] = 0;
    dst_sizes_meta[0] = sizeof(double);

    status = H5TBappend_records(config.file_id, "Probability", n_iter,
				sizeof(double), dst_offsets_meta, 
				dst_sizes_meta, config.probability);
    
    dst_sizes_meta[0] = sizeof(int);
    status = H5TBappend_records(config.file_id, "Accepted", n_iter,
				sizeof(int), dst_offsets_meta, 
				dst_sizes_meta, config.accepted);
}
