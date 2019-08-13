
#include "simulated_data.h"

char **simulated_data_data_file_names( )
{
	char **file;
	file = calloc( 5 , sizeof( char *) );
	
	file[0] = "../SIM_small_data/obs_freq_23Hz.fits";
	file[1] = "../SIM_small_data/obs_freq_33Hz.fits";
	file[2] = "../SIM_small_data/obs_freq_41Hz.fits";
	file[3] = "../SIM_small_data/obs_freq_61Hz.fits";
	file[4] = "../SIM_small_data/obs_freq_94Hz.fits";
	
	return( file );
}

char **simulated_data_template_file_names( )
{
	char **file;
	file = calloc( 4 , sizeof( char *) );
	
	file[0] = "../SIM_small_data/template_1.fits";
	file[1] = "../SIM_small_data/template_2.fits";
	file[2] = "../SIM_small_data/template_3.fits";
	file[3] = "../SIM_small_data/template_4.fits";
	
	return( file );
}

char *simulated_data_temperature_analysis_mask_name( )
{
	char *file;
	file = "../SIM_small_data/mask.fits" ;
	return( file );
}


char **simulated_data_result_file_names( )
{
	char **file; 
	file = calloc( 19 , sizeof(char *));
	
	file[0] = "../SIM_small_result/source_1.fits";
	file[1] = "../SIM_small_result/source_2.fits";
	file[2] = "../SIM_small_result/source_3.fits";
	file[3] = "../SIM_small_result/source_4.fits";
	file[4] = "../SIM_small_result/source_1_sdev.fits";
	file[5] = "../SIM_small_result/source_2_sdev.fits";
	file[6] = "../SIM_small_result/source_3_sdev.fits";
	file[7] = "../SIM_small_result/source_4_sdev.fits";
	file[8] = "../SIM_small_result/sync_ind.fits";
	file[9] = "../SIM_small_result/dust_ind.fits";
	file[10] = "../SIM_small_result/source_1_precision.fits";
	file[11] = "../SIM_small_result/source_2_precision.fits";
	file[12] = "../SIM_small_result/source_3_precision.fits";
	file[13] = "../SIM_small_result/source_4_precision.fits";
	file[14] = "../SIM_small_result/obs_23_residual.fits";
	file[15] = "../SIM_small_result/obs_33_residual.fits";
	file[16] = "../SIM_small_result/obs_41_residual.fits";
	file[17] = "../SIM_small_result/obs_61_residual.fits";
	file[18] = "../SIM_small_result/obs_94_residual.fits";
	
	return( file );
}

void simulated_data_band_frequency_and_error_precisions( double *obs_freq, double *obs_precision )
{
	obs_freq[0] = 23.; obs_freq[1] = 33.; obs_freq[2] = 41. ;
	obs_freq[3] = 61.; obs_freq[4] = 94.;
	
	int k;
	
	for( k=0; k<5; k++ ) obs_precision[k] = 1./(.001*.001);
	
	//precision in units of (mK)^{-2}  (milli Kelvin)
	/*obs_precision[0] = 1./gsl_pow_2(.2); obs_precision[1] = 1./gsl_pow_2(.2);
	obs_precision[2] = 1./gsl_pow_2(.2); obs_precision[3] = 1./gsl_pow_2(.2);
	obs_precision[4] = 1./gsl_pow_2(.2);*/
	
}



