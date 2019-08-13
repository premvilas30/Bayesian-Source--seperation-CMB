//data.h : functions for data/result I/O

#ifndef __DATA_H__
#define __DATA_H__

#include "fitsio.h"
#include "hyperpar.h"
#include "block.h"
#include "patch.h"

int *data_get_mask_ones( int n_pixel );

int *data_get_mask( int id, int n_pixel, char *file_name);

void data_read_maps( int id, struct block *b, struct patch **p, struct hyperpar *h, char **file, int in_type, double *obs_intensity, double *hit_rate, int explo );

void data_read_templates( int id, struct block *b, struct patch **p, char **template_file, int explo, int prior );

void data_healpix_reorder( struct block *b, double *y, double *n_hit, int ascii, int ordering ) ;

void data_put_values( struct block *b, struct patch **p, struct hyperpar *h, int map, double *y, double *z, int explo ) ;

void data_put_spec_values(  struct block *b, struct patch **p, double *v, int synch);

void data_put_template( struct block *b, struct patch **p, int map, double *t, int prior);

void data_compute_det_C( struct block *b, struct patch **p );

void data_write_out_patch_mean( int patch_idx, int nout_map, struct patch *p, struct block *b, char *file_name );


#endif
