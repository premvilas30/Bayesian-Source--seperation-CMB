#include "nelder_mead_simplex.h"

//code adapted from R source to do Nelder-Mead simplex optimisation

//Bvec is initial pars
// X is out pars
// Fmin is the value of the function at min

void optimizer_nmmin(int n, double *Bvec, double *X, double *Fmin,
	   int *fail, double abstol, double intol, 
	   double alpha, double bet, double gamm, int trace,
	   int *fncount, int maxit, void *params )
{
    char action[50];
    int C;
    int calcvert;
    double convtol, f;
    int funcount=0, H, i, j, L=0;
    int n1=0;
    double oldsize;
    double **P;
    double size, step, temp, trystep;
    char tstr[9]; // allow for 10^8 iters ...
    double VH, VL, VR;
    
    gsl_vector *x;
    x = gsl_vector_alloc( n ); 

    if (maxit <= 0) {
    for( i=0; i<n; i++ ) gsl_vector_set( x, i, Bvec[i] );
	*Fmin = gsl_objective_function( x , params );//fminfn(n, Bvec, ex);
	*fncount = 0;
	*fail = 0;
	return;
    }
    if (trace)
	printf("\n  Nelder-Mead direct search function minimizer\n");
    P = calloc( n+1, sizeof(double*));
    for( i=0; i<n+1; i++ ) P[i] = calloc(n+2,sizeof(double));
   // matrix(n, n+1);
    *fail = FALSE;
    for( i=0; i<n; i++ ) 
    {
    	gsl_vector_set( x, i, Bvec[i] );
    	//printf("\n value %d of x is %.20f", i, Bvec[i]);
    }
    f = gsl_objective_function( x , params );
    if (!isfinite(f)) {
	printf("\nfunction cannot be evaluated at initial parameters");
	*fail = TRUE;
    } else {
	if (trace) printf("\nfunction value for initial parameters = %f\n", f);
	funcount = 1;
	convtol = intol * (fabs(f) + intol);
	if (trace) printf("\n  Scaled convergence tolerance is %g\n", convtol);
	n1 = n + 1;
	C = n + 2;
	P[n1 - 1][0] = f;
	for (i = 0; i < n; i++)
	    P[i][0] = Bvec[i];

	L = 1;
	size = 0.0;

	step = 0.0;
	for (i = 0; i < n; i++) {
	    if (0.1 * fabs(Bvec[i]) > step)
		step = 0.1 * fabs(Bvec[i]);
	}
	if (step == 0.0) step = 0.1;
	if (trace) printf("\nStepsize computed as %f\n", step);
	for (j = 2; j <= n1; j++) {
	    strcpy(action, "BUILD          ");
	    for (i = 0; i < n; i++)
		P[i][j - 1] = Bvec[i];

	    trystep = step;
	    while (P[j - 2][j - 1] == Bvec[j - 2]) {
		P[j - 2][j - 1] = Bvec[j - 2] + trystep;
		trystep *= 10;
	    }
	    size += trystep;
	}
	oldsize = size;
	calcvert = TRUE;
	do {
	    if (calcvert) {
		for (j = 0; j < n1; j++) {
		    if (j + 1 != L) {
			for (i = 0; i < n; i++){
			    Bvec[i] = P[i][j];
			    gsl_vector_set( x, i, Bvec[i] );
			  }
			 f =  gsl_objective_function( x , params );
			//f = fminfn(n, Bvec, ex);
			if (!isfinite(f)) f = big;
			funcount++;
			P[n1 - 1][j] = f;
		    }
		}
		calcvert = FALSE;
	    }

	    VL = P[n1 - 1][L - 1];
	    VH = VL;
	    H = L;

	    for (j = 1; j <= n1; j++) {
		if (j != L) {
		    f = P[n1 - 1][j - 1];
		    if (f < VL) {
			L = j;
			VL = f;
		    }
		    if (f > VH) {
			H = j;
			VH = f;
		    }
		}
	    }

	    if (VH <= VL + convtol || VL <= abstol) break;

	    // avoid buffer overflow at 100001 iters. (PR#15240)
	    if (trace) {
		snprintf(tstr, 9, "%5d", funcount);
		printf("\n%s%s %f %f\n", action, tstr, VH, VL);
	    }

	    for (i = 0; i < n; i++) {
		temp = -P[i][H - 1];
		for (j = 0; j < n1; j++)
		    temp += P[i][j];
		P[i][C - 1] = temp / n;
	    }
	    for (i = 0; i < n; i++){
			Bvec[i] = (1.0 + alpha) * P[i][C - 1] - alpha * P[i][H - 1];
			gsl_vector_set( x, i, Bvec[i] );
		}
		
		f = gsl_objective_function( x , params );
	    //f = fminfn(n, Bvec, ex);
	    if (!isfinite(f)) f = big;
	    funcount++;
	    strcpy(action, "REFLECTION     ");
	    VR = f;
	    if (VR < VL) {
		P[n1 - 1][C - 1] = f;
		for (i = 0; i < n; i++) {
		    f = gamm * Bvec[i] + (1 - gamm) * P[i][C - 1];
		    P[i][C - 1] = Bvec[i];
		    Bvec[i] = f;
		    gsl_vector_set( x, i, Bvec[i] );
		}
		f = gsl_objective_function( x , params );
		//f = fminfn(n, Bvec, ex);
		if (!isfinite(f)) f = big;
		funcount++;
		if (f < VR) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = f;
		    strcpy(action, "EXTENSION      ");
		} else {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = P[i][C - 1];
		    P[n1 - 1][H - 1] = VR;
		}
	    } else {
		strcpy(action, "HI-REDUCTION   ");
		if (VR < VH) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = VR;
		    strcpy(action, "LO-REDUCTION   ");
		}

		for (i = 0; i < n; i++){
		    Bvec[i] = (1 - bet) * P[i][H - 1] + bet * P[i][C - 1];
		 	 gsl_vector_set( x, i, Bvec[i] );
		 }
		f = gsl_objective_function( x , params );
		//f = fminfn(n, Bvec, ex);
		if (!isfinite(f)) f = big;
		funcount++;

		if (f < P[n1 - 1][H - 1]) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = f;
		} else {
		    if (VR >= VH) {
			strcpy(action, "SHRINK         ");
			calcvert = TRUE;
			size = 0.0;
			for (j = 0; j < n1; j++) {
			    if (j + 1 != L) {
				for (i = 0; i < n; i++) {
				    P[i][j] = bet * (P[i][j] - P[i][L - 1])
					+ P[i][L - 1];
				    size += fabs(P[i][j] - P[i][L - 1]);
				}
			    }
			}
			if (size < oldsize) {
			    oldsize = size;
			} else {
			    if (trace)
				printf("\nPolytope size measure not decreased in shrink\n");
			    *fail = 10;
			    break;
			}
		    }
		}
	    }

	//printf("\n made it to  the end of the loop...");
	} while (funcount <= maxit);

    }

    if (trace) {
	printf("Exiting from Nelder Mead minimizer\n");
	printf("    %d function evaluations used\n", funcount);
    }
    *Fmin = P[n1 - 1][L - 1];

    for (i = 0; i < n; i++) X[i] = P[i][L - 1];

    if (funcount > maxit) *fail = 1;
    *fncount = funcount;

    for( i=0; i<n+1; i++ ) free(P[i]);
    free(P);
      
}
