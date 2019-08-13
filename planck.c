// planck.c :: functions for the wmap data...
// Planck data is measured in Kelvin

#include "planck.h"

void planck_get_data_dims( int *n_freq_band, int *n_src, int *nside )
{
	*n_freq_band = 7;
	*n_src = 4;
	*nside = 1024;
	return;
}

char **planck_get_data_file_names( )
{
	char **file;
	file = calloc( 7 , sizeof( char *) );
	
	file[0] = "../PLANCK_data/LFI_SkyMap_030_1024_ANALYSIS.fits";
	file[1] = "../PLANCK_data/LFI_SkyMap_044_1024_ANALYSIS.fits";
	file[2] = "../PLANCK_data/LFI_SkyMap_070_1024_ANALYSIS.fits";
	file[3] = "../PLANCK_data/HFI_SkyMap_100_1024_ANALYSIS.fits";
	file[4] = "../PLANCK_data/HFI_SkyMap_143_1024_ANALYSIS.fits";
	file[5] = "../PLANCK_data/HFI_SkyMap_217_1024_ANALYSIS.fits";
	file[6] = "../PLANCK_data/HFI_SkyMap_353_1024_ANALYSIS.fits";
	file[7] = "../PLANCK_data/HFI_SkyMap_545_1024_ANALYSIS.fits";
	file[8] = "../PLANCK_data/HFI_SkyMap_857_1024_ANALYSIS.fits";
	
	return( file );
}

char **planck_get_result_file_names( )
{
	char **file; 
	file = calloc( 21 , sizeof(char *));
	
	file[0] = "../PLANCK_result/cmb.fits";
	file[1] = "../PLANCK_result/sync.fits";
	file[2] = "../PLANCK_result/dust.fits";
	file[3] = "../PLANCK_result/ffem.fits";
	file[4] = "../PLANCK_result/cmb_sdev.fits";
	file[5] = "../PLANCK_result/sync_sdev.fits";
	file[6] = "../PLANCK_result/dust_sdev.fits";
	file[7] = "../PLANCK_result/ffem_sdev.fits";
	file[8] = "../PLANCK_result/sync_ind.fits";
	file[9] = "../PLANCK_result/dust_ind.fits";
	file[10] = "../PLANCK_result/cmb_precision.fits";
	file[11] = "../PLANCK_result/sync_precision.fits";
	file[12] = "../PLANCK_result/dust_precision.fits";
	file[13] = "../PLANCK_result/ffem_precision.fits";
	file[14] = "../PLANCK_result/SkyMap_030_residual.fits";
	file[15] = "../PLANCK_result/SkyMap_044_residual.fits";
	file[16] = "../PLANCK_result/SkyMap_070_residual.fits";
	file[17] = "../PLANCK_result/SkyMap_100_residual.fits";
	file[18] = "../PLANCK_result/SkyMap_143_residual.fits";
	file[19] = "../PLANCK_result/SkyMap_217_residual.fits";
	file[20] = "../PLANCK_result/SkyMap_353_residual.fits";
	//file[21] = "../PLANCK_result/SkyMap_545_residual.fits";	
	//file[22] = "../PLANCK_result/SkyMap_857_residual.fits";
	return( file );
}

char *planck_get_temperature_analysis_mask_name( )
{
	char *file;
	file = "../PLANCK_data/Temperature_mask_1024.fits" ;
	return( file );
}

char *planck_get_template_file_names()
{
	char *file;
	return( file );
}	

void planck_get_band_frequency_and_error_precisions( double *obs_freq, double *obs_precision )
{
	obs_freq[0] = 30.; obs_freq[1] = 44.; obs_freq[2] = 70. ;
	obs_freq[3] = 100.; obs_freq[4] = 143.; obs_freq[5] = 260.;
	obs_freq[6] = 430.; //obs_freq[7] = 630.; obs_freq[8] = 857.;
	
	//precision in units of (mK)^{-2} (Kelvin)
	obs_precision[0] = .34600379;  //629881.6;
	obs_precision[1] = .11792586;  //694444.4;
	obs_precision[2] = .04553005;  //783146.7;
	obs_precision[3] = .40643001; //12755102.0;
	obs_precision[4] = 1.63993905;  //30864197.5;
	obs_precision[5] = .85265379;  // 30864197.5;
	obs_precision[6] = .07422005;  //30864197.5;

	
	//apply a correction (to make computation comparable with WMAP)
	
	//int i;
	//for( i=0; i<7; i++ ) obs_precision[i] *= 1E-6;
	
	//obs_precision[7] = 30864197.5;
	//obs_precision[8] = 30864197.5;	
	
}



