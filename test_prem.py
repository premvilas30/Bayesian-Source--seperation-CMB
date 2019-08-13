import subprocess

cmd = ['gcc', '-o', 'test.out', '-pg', 'chealpix.c', 'wmap.c', 'planck.c', 'simulated_data.c', 'graphs.c', 'precisions.c', 'hyperpar.c', 'block.c', 'routines.c', 'patch.c', 'update.c', 'data.c', 'super.c', 'model.c', 'model_gsl.c', 'nelder_mead_simplex.c', 'optimizer.c', 'results.c', 'grid.c', 'analysis_wmap_data_9yr.c', '-fopenmp', '-lm', '-lgsl', '-lgslcblas', '-lcholmod', '-lamd', '-lcamd', '-lcolamd', '-lsuitesparseconfig', '-lcfitsio', '-pg', 'hello']

p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

print (p.communicate())

