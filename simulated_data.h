
#ifndef __SIMULATED_DATA_H__
#define __SIMULATED_DATA_H__

#include "required_libs.h"

char **simulated_data_data_file_names( ) ;

char *simulated_data_temperature_analysis_mask_name( ) ;

char **simulated_data_template_file_names( );

char **simulated_data_result_file_names( ) ;

void simulated_data_band_frequency_and_error_precisions( double *obs_freq, double *obs_precision ) ;


#endif
