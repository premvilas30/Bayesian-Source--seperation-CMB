
#ifndef __WMAP_H__
#define __WMAP_H__

#include "required_libs.h"

void wmap_get_7yr_data_dims( int *n_freq_band, int *n_src, int *nside );

void wmap_get_9yr_data_dims( int *n_freq_band, int *n_src, int *nside );

char **wmap_get_7yr_data_file_names( );

char **wmap_get_9yr_data_file_names( );

char **wmap_get_7yr_result_file_names( );

char **wmap_get_9yr_result_file_names( );

char *wmap_get_7yr_temperature_analysis_mask_name( );

char *wmap_get_9yr_temperature_analysis_mask_name( );

char **wmap_get_7yr_template_file_names( );

char **wmap_get_9yr_template_file_names( );

void wmap_get_7yr_band_frequency_and_error_precisions_and_offsets( double *obs_freq, double *obs_precision, double *offsets );

void wmap_get_9yr_band_frequency_and_error_precisions_and_offsets( double *obs_freq, double *obs_precision, double *offset );

void wmap_get_7yr_prior_source_dims( int n_src, int *prior_mu, double *prior_mu_ref_freq );

void wmap_get_9yr_prior_source_dims( int n_src, int *prior_mu, double *prior_mu_ref_freq );

#endif
