
#ifndef __SIMULATE_H__
#define __SIMULATE_H__

#include "required_libs.h"
#include "defs.h"
#include "simulateIGMRF.h"
#include "super.h"


struct simulation
{
	int n_side;
	int n_freq; 
	int n_src;
	
	double **src;
	double **obs;
		
};

struct simulation *simulation_create( int n_side, int n_freq, int n_src ) ;

void simulation_destroy( struct simulation *sim );

void simulation_simulate_data( int n_side, int n_freq, int n_src, double *nu, double *specpar, double *obs_precisions, double *src_precisions, char **obs_files, char **src_files, char *info_file  ) ;



#endif
