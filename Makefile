INC=-I/usr/include
SOURCE_TEST_PATCH = test_patch.c
SOURCE_TEST_BLOCK = test_block.c
SOURCE_TEST_SUPER = test_super.c
SOURCE_TEST_R = CMB_analyze.c
SOURCE_TEST_CMB = test_cmb.c
SOURCE_ANALYSIS_WMAP = analysis_wmap_data.c
SOURCE_ANALYSIS_WMAP_9YR = analysis_wmap_data_9yr.c
SOURCE_ANALYSIS_PLANCK = analysis_planck_data.c
SOURCE_ANALYSIS_SIMULATED = analysis_simulated_data_new.c
SOURCE_CREATE_SIMULATED = create_simulated_data.c

TARGET = test.out
FLAGS =  -pg #-Wall -Wextra
GCC = gcc
RGCC = R CMD SHLIB
LIBS =  -fopenmp -lm -lgsl -lgslcblas -lcholmod -lamd -lcamd -lcolamd -lsuitesparseconfig -lcfitsio #-g 
#MYLIBS = wmap.c planck.c simulated_data.c graphs.c precisions.c hyperpar.c block.c routines.c patch.c update.c data.c super.c model.c model_gsl.c optimizer.c results.c grid.c simulateIGMRF.c simulate.c

MYLIBS = chealpix.c wmap.c planck.c simulated_data.c graphs.c precisions.c hyperpar.c block.c routines.c patch.c update.c data.c super.c  model.c model_gsl.c nelder_mead_simplex.c optimizer.c results.c grid.c 

MYLIBS_TEST_BLOCK = 

patch_test:
	$(GCC) -o $(TARGET)  $(FLAGS) $(MYLIBS) $(SOURCE_TEST_PATCH) $(LIBS) -pg

block_test:
	$(GCC) -o $(TARGET)  $(FLAGS) $(MYLIBS) $(SOURCE_TEST_BLOCK) $(LIBS) -pg
	
super_test:
	$(GCC) -o $(TARGET)  $(FLAGS) $(MYLIBS) $(SOURCE_TEST_SUPER) $(LIBS) -pg

cmb_test:
	$(GCC) -o $(TARGET)  $(FLAGS) $(MYLIBS) $(SOURCE_TEST_CMB) $(LIBS) -pg
	
wmap_analysis:
	$(GCC) -o $(TARGET)  $(FLAGS) $(MYLIBS) $(SOURCE_ANALYSIS_WMAP) $(LIBS) -pg
	
wmap_analysis_9yr:
	$(GCC) -o $(TARGET)  $(FLAGS) $(MYLIBS) $(SOURCE_ANALYSIS_WMAP_9YR) $(LIBS) -pg

planck_anlaysis:
	$(GCC) -o $(TARGET)  $(FLAGS) $(MYLIBS) $(SOURCE_ANALYSIS_PLANCK) $(LIBS) -pg
	
r_test:
	R CMD SHLIB -o CMB_analyze.so $(FLAGS) $(MYLIBS) $(SOURCE_TEST_R) $(LIBS)
	
create_simulated:
	$(GCC) -o $(TARGET) $(FLAGS) $(MYLIBS) $(SOURCE_CREATE_SIMULATED) $(LIBS) -pg
	
simulated_analysis:
	$(GCC) -o $(TARGET) $(FLAGS) $(MYLIBS) $(SOURCE_ANALYSIS_SIMULATED) $(LIBS) -pg

R_interface:
	$(RGCC) -o Rinterface.so $(FLAGS) $(MYLIBS) interface_for_R.c $(LIBS) 

clean:
	rm $(TARGET)
	touch *

