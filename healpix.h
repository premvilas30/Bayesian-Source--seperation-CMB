/*
	healpix.c
	
	Functions to do operations on Healpix orderings-- to allow for index referencing
	
	Written by: Jason Wyse,
				School of Computer Science and Statistics,
				Trinity College Dublin,
				Dublin 2, 
				Ireland
			
	Last modified: Thu 21 Jul 2016 16:38:26 IST 
	
*/

#ifndef __HEALPIX_H__
#define __HEALPIX_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>

void healpix_convert_nested_to_2dlattice_ordering( int nside, int *lattice_pix_number );


#endif
