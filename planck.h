
#ifndef __PLANCK_H__
#define __PLANCK_H__

#include "required_libs.h"

void planck_get_data_dims( int *n_freq_band, int *n_src, int *nside );

char **planck_get_data_file_names( );

char **planck_get_result_file_names( );

char *planck_get_template_file_names( );

char *planck_get_temperature_analysis_mask_name( );

void planck_get_band_frequency_and_error_precisions( double *obs_freq, double *obs_precision );

#endif
