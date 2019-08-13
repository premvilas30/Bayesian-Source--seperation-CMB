// wmap.c :: functions for the wmap data...
// WMAP data is measured in milli Kelvin

#include "wmap.h"

void wmap_get_7yr_data_dims( int *n_freq_band, int *n_src, int *nside )
{
	*n_freq_band = 5;
	*n_src = 4;
	*nside = 512;
	return;
}

void wmap_get_9yr_data_dims( int *n_freq_band, int *n_src, int *nside )
{
	*n_freq_band = 5;
	*n_src = 4;
	*nside = 512;
	return;
}

char **wmap_get_7yr_data_file_names( )
{
	char **file;
	file = calloc( 5 , sizeof( char *) );
	
	file[0] = "../WMAP_data/wmap_band_imap_r9_7yr_K_v4.fits";
	file[1] = "../WMAP_data/wmap_band_imap_r9_7yr_Ka_v4.fits";
	file[2] = "../WMAP_data/wmap_band_imap_r9_7yr_Q_v4.fits";
	file[3] = "../WMAP_data/wmap_band_imap_r9_7yr_V_v4.fits";
	file[4] = "../WMAP_data/wmap_band_imap_r9_7yr_W_v4.fits";
	
	//trial
	//file[5] = "../TEMPLATES/haslam_NESTED.fits";
	//file[6] = "../TEMPLATES/finkbeiner_dust_NESTED.fits";
	
	return( file );
}

char **wmap_get_9yr_data_file_names( )
{
	char **file;
	file = calloc( 7 , sizeof( char *) );
	
	file[0] = "../WMAP_data/wmap_band_imap_r9_9yr_K_v5.fits";
	file[1] = "../WMAP_data/wmap_band_imap_r9_9yr_Ka_v5.fits";
	file[2] = "../WMAP_data/wmap_band_imap_r9_9yr_Q_v5.fits";
	file[3] = "../WMAP_data/wmap_band_imap_r9_9yr_V_v5.fits";
	file[4] = "../WMAP_data/wmap_band_imap_r9_9yr_W_v5.fits";
	//trial
	//file[5] = "../TEMPLATES/haslam_NESTED.fits";
	//file[6] = "../TEMPLATES/finkbeiner_dust_NESTED.fits";
	
	return( file );
}

char **wmap_get_7yr_result_file_names( )
{
	char **file; 
	file = calloc( 21 , sizeof(char *));
	
	file[0] = "../WMAP_result/cmb_7yr.fits";
	file[1] = "../WMAP_result/sync_7yr.fits";
	file[2] = "../WMAP_result/dust_7yr.fits";
	file[3] = "../WMAP_result/ffem_7yr.fits";
	file[4] = "../WMAP_result/cmb_sdev_7yr.fits";
	file[5] = "../WMAP_result/sync_sdev_7yr.fits";
	file[6] = "../WMAP_result/dust_sdev_7yr.fits";
	file[7] = "../WMAP_result/ffem_sdev_7yr.fits";
	file[8] = "../WMAP_result/sync_ind_7yr.fits";
	file[9] = "../WMAP_result/dust_ind_7yr.fits";
	file[10] = "../WMAP_result/cmb_precision_7yr.fits";
	file[11] = "../WMAP_result/sync_precision_7yr.fits";
	file[12] = "../WMAP_result/dust_precision_7yr.fits";
	file[13] = "../WMAP_result/ffem_precision_7yr.fits";
	file[14] = "../WMAP_result/K_residual_7yr.fits";
	file[15] = "../WMAP_result/Ka_residual_7yr.fits";
	file[16] = "../WMAP_result/Q_residual_7yr.fits";
	file[17] = "../WMAP_result/V_residual_7yr.fits";
	file[18] = "../WMAP_result/W_residual_7yr.fits";
	file[19] = "../WMAP_result/temp_sync_residual_7yr.fits";
	file[20] = "../WMAP_result/temp_dust_residual_7yr.fits";
	
	return( file );
}

char **wmap_get_9yr_result_file_names( )
{
	char **file; 
	file = calloc( 21 , sizeof(char *));
	
	file[0] = "../WMAP_result/cmb_9yr.fits";
	file[1] = "../WMAP_result/sync_9yr.fits";
	file[2] = "../WMAP_result/dust_9yr.fits";
	file[3] = "../WMAP_result/ffem_9yr.fits";
	file[4] = "../WMAP_result/cmb_sdev_9yr.fits";
	file[5] = "../WMAP_result/sync_sdev_9yr.fits";
	file[6] = "../WMAP_result/dust_sdev_9yr.fits";
	file[7] = "../WMAP_result/ffem_sdev_9yr.fits";
	file[8] = "../WMAP_result/sync_ind_9yr.fits";
	file[9] = "../WMAP_result/dust_ind_9yr.fits";
	file[10] = "../WMAP_result/cmb_precision_9yr.fits";
	file[11] = "../WMAP_result/sync_precision_9yr.fits";
	file[12] = "../WMAP_result/dust_precision_9yr.fits";
	file[13] = "../WMAP_result/ffem_precision_9yr.fits";
	file[14] = "../WMAP_result/K_residual_9yr.fits";
	file[15] = "../WMAP_result/Ka_residual_9yr.fits";
	file[16] = "../WMAP_result/Q_residual_9yr.fits";
	file[17] = "../WMAP_result/V_residual_9yr.fits";
	file[18] = "../WMAP_result/W_residual_9yr.fits";
	file[19] = "../WMAP_result/temp_sync_residual_9yr.fits";
	file[20] = "../WMAP_result/temp_dust_residual_9yr.fits";
	
	return( file );
}

char *wmap_get_7yr_temperature_analysis_mask_name( )
{
	char *file;
	file = "../WMAP_data/wmap_temperature_analysis_mask_r9_7yr_v4.fits" ;
	return( file );
}

char *wmap_get_9yr_temperature_analysis_mask_name( )
{
	char *file;
	file = "../WMAP_data/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits" ;
	return( file );
}

char **wmap_get_7yr_template_file_names( )
{
	int k;
	char **template_files =  calloc( 4, sizeof(char*) ) ;
	
	template_files[0] = "../TEMPLATES/wmap_ilc_9yr_v5.fits"; //cmb
	template_files[1] = "../TEMPLATES/haslam_NESTED.fits"; //sync
	template_files[2] = "../TEMPLATES/wmap_temp_dust_smth_template_v3.fits"; //dust
	template_files[3] = "../TEMPLATES/wmap_halpha_template_v5.fits"; //free-free
	
	return( template_files );
	 
}

char **wmap_get_9yr_template_file_names( )
{
	int k;
	char **template_files =  calloc( 4, sizeof(char*) ) ;
	
	template_files[0] = "../TEMPLATES/wmap_ilc_9yr_v5.fits";
	template_files[1] = "../TEMPLATES/lambda_haslam408_dsds.fits";
	template_files[2] = "../TEMPLATES/wmap_temp_dust_smth_template_v3.fits";
	template_files[3] = "../TEMPLATES/wmap_halpha_smth_template_v3.fits";
	
	return( template_files );
	 
}

void wmap_get_7yr_band_frequency_and_error_precisions_and_offsets( double *obs_freq, double *obs_precision, double *offset )
{
	obs_freq[0] = 23.; obs_freq[1] = 33.; obs_freq[2] = 41. ;
	obs_freq[3] = 61.; obs_freq[4] = 94.; obs_freq[5] = .408;
	obs_freq[6] = 94.;
	
	//precision in units of (mK)^{-2}  (milli Kelvin)
	obs_precision[0] = 1./gsl_pow_2(1.437); obs_precision[1] = 1./gsl_pow_2(1.470);
	obs_precision[2] = 1./gsl_pow_2(2.197); obs_precision[3] = 1./gsl_pow_2(3.137);
	obs_precision[4] = 1./gsl_pow_2(6.549); //obs_precision[5] = 100.;
	//obs_precision[6] = 100.;
	
	//offsets in units of mK
	//offset[0] = -.025; offset[1] = -.0054; offset[2] = -.0024; offset[3] = -.0024; offset[4] = -.0015;
	
	
}


void wmap_get_9yr_band_frequency_and_error_precisions_and_offsets( double *obs_freq, double *obs_precision, double *offset )
{
	obs_freq[0] = 23.; obs_freq[1] = 33.; obs_freq[2] = 41. ;
	obs_freq[3] = 61.; obs_freq[4] = 94.; 
	
	//templates
	//obs_freq[5] = .408;
	//obs_freq[6] = 94.;
	
	//precision in units of (mK)^{-2}  (milli Kelvin)
	obs_precision[0] = 1./gsl_pow_2(1.437); obs_precision[1] = 1./gsl_pow_2(1.470);
	obs_precision[2] = 1./gsl_pow_2(2.197); obs_precision[3] = 1./gsl_pow_2(3.137);
	obs_precision[4] = 1./gsl_pow_2(6.549); //obs_precision[5] = 100.;
	//obs_precision[6] = 100.;
	
	//offsets in units of mK
	int k;
	for( k=0; k<5; k++ ) offset[k] = 0.; 
		
}

void wmap_get_7yr_prior_source_dims( int n_src, int *prior_mu, double *prior_mu_ref_freq )
{
	//all four sources have  priors
	int k;
	for( k=0; k<n_src; k++ ) prior_mu[k] = 1;
	
	prior_mu_ref_freq[0] = 23.;
	prior_mu_ref_freq[1] = 23.;
	prior_mu_ref_freq[2] = 94.;
	prior_mu_ref_freq[3] = 23.;
	
	return; 
}

void wmap_get_9yr_prior_source_dims( int n_src, int *prior_mu, double *prior_mu_ref_freq )
{
	//all four sources have  priors
	int k;
	for( k=0; k<n_src; k++ ) prior_mu[k] = 1;
	
	prior_mu_ref_freq[0] = 23.;
	prior_mu_ref_freq[1] = 23.;
	prior_mu_ref_freq[2] = 94.;
	prior_mu_ref_freq[3] = 23.;
	
	return; 
}

