/*------------------------------------------------------------------------
 *   Create dissipative boundary around the model grid
 *   The dissipative coefficients are stored in the 2-D array
 *   absorb_coeff. The interior of the model is weighted by the
 *   coefficient 1. In the absorbing frame the coefficients 
 *   are less than one. Coefficients are computed using 
 *   exponential damping (see Cerjan et al., 1985, 
 *   Geophysics, 50, 705-708)
 *
 *    T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void absorb(float *** absorb_coeff){

	/* extern variables */
	extern float	DAMPLE, DAMPRI, DAMPBA, DAMPFR, DAMPTO, DAMPBO;
	extern int	NX[3];
	extern int	AB_LE, AB_RI, AB_BA, AB_FR, AB_TO, AB_BO, FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern int	MYID;
	extern FILE	*FP;
	
	/* local variables */
	int	i, j, k, ifwle, ifwri, ifwba, ifwfr, ifwto, ifwbo, ii=0, jj=0, kk=0;
	float	amp, a, *coeff;

	ifwle = FWLE;
	ifwri = FWRI;
	ifwba = FWBA;
	ifwfr = FWFR;
	ifwto = FWTO;
	ifwbo = FWBO;

	fprintf(FP,"\n\n ----------------------- ABSORBING BOUNDARY ------------------------\n");
	fprintf(FP,"\n **Message from absorb (printed by PE %d):\n",MYID);
	fprintf(FP," Coefficients for absorbing frame are now calculated.\n");

	/* initialize array of coefficients with one */
	for (i=1;i<=NX[0];i++)
		for (j=1;j<=NX[1];j++)
			for (k=1;k<=NX[2];k++)
				absorb_coeff[i][j][k] = 1.0;


	/* compute coefficients for left grid boundary (x-direction) */
	if (AB_LE){
		fprintf(FP," Width of absorbing frame in x-direction (left  = %i gridpoints.\n",ifwle);
		fprintf(FP," Percentage of exponential damping (left)  = %5.2f\n",DAMPLE);
		amp = 1.0-DAMPLE/100.0;   /* amplitude at the edge of the numerical grid */
		coeff = vector(1,ifwle);
		a = sqrt(-log(amp)/((ifwle-1)*(ifwle-1)));
		for (i=1;i<=ifwle;i++)
			coeff[i] = exp(-(a*a*(ifwle-i)*(ifwle-i)));

		fprintf(FP," Table of coefficients in x-direction (left) \n # \t coeffle \n");
		for (i=1;i<=ifwle;i++)
			fprintf(FP," %d \t %5.3f \n", i, coeff[i]);

		for (i=1;i<=ifwle;i++){
			for (j=1;j<=NX[1];j++){
				for (k=1;k<=NX[2];k++){
					absorb_coeff[i][j][k] *= coeff[i];
				}
			}
		}
		free_vector(coeff,1,ifwle);
	}

	/* compute coefficients for right grid boundary (x-direction) */
	if (AB_RI){
		fprintf(FP," Width of absorbing frame in x-direction (right  = %i gridpoints.\n",ifwri);
		fprintf(FP," Percentage of exponential damping (right)  = %5.2f\n",DAMPRI);
		amp = 1.0-DAMPRI/100.0;   /* amplitude at the edge of the numerical grid */
		coeff = vector(1,ifwri);
		a = sqrt(-log(amp)/((ifwri-1)*(ifwri-1)));
		for (i=1;i<=ifwri;i++)
			coeff[i] = exp(-(a*a*(ifwri-i)*(ifwri-i)));

		fprintf(FP," Table of coefficients in x-direction (right) \n # \t coeffri \n");
		for (i=1;i<=ifwri;i++)
			fprintf(FP," %d \t %5.3f \n", i, coeff[i]);

		for (i=1;i<=ifwri;i++){
			ii = NX[0]-i+1;
			for (j=1;j<=NX[1];j++){
				for (k=1;k<=NX[2];k++){
					absorb_coeff[ii][j][k] *= coeff[i];
				}
			}
		}
		free_vector(coeff,1,ifwri);
	}

	/* compute coefficients for back grid boundary (y-direction) */
	if (AB_BA){
		fprintf(FP," Width of absorbing frame in y-direction (back)   = %i gridpoints.\n",ifwba);
		fprintf(FP," Percentage of exponential damping (back)    = %5.2f\n",DAMPBA);
		amp = 1.0-DAMPBA/100.0;   /* amplitude at the edge of the numerical grid */
		coeff = vector(1,ifwba);
		a = sqrt(-log(amp)/((ifwba-1)*(ifwba-1)));
		for (j=1;j<=ifwba;j++)
			coeff[j] = exp(-(a*a*(ifwba-j)*(ifwba-j)));
		fprintf(FP," Table of coefficients in y-direction (back) \n # \t coeffba \n");
		for (j=1;j<=ifwba;j++)
			fprintf(FP," %d \t %5.3f \n", j, coeff[j]);

		for (i=1;i<=NX[0];i++){
			for (j=1;j<=ifwba;j++){
				for (k=1;k<=NX[2];k++){
					absorb_coeff[i][j][k] *= coeff[j];
				}
			}
		}
		free_vector(coeff,1,ifwba);
	}

	/* compute coefficients for front grid boundary (y-direction) */
	if (AB_FR){
		fprintf(FP," Width of absorbing frame in y-direction (front)= %i gridpoints.\n",ifwfr);
		fprintf(FP," Percentage of exponential damping (front) = %5.2f\n",DAMPFR);
		amp = 1.0-DAMPFR/100.0;   /* amplitude at the edge of the numerical grid */
		coeff = vector(1,ifwfr);
		a = sqrt(-log(amp)/((ifwfr-1)*(ifwfr-1)));
		for (j=1;j<=ifwfr;j++)
			coeff[j] = exp(-(a*a*(ifwfr-j)*(ifwfr-j)));
		fprintf(FP," Table of coefficients in y-direction (front) \n # \t coefffr \n");
		for (j=1;j<=ifwfr;j++)
			fprintf(FP," %d \t %5.3f \n", j, coeff[j]);

		for (i=1;i<=NX[0];i++){
			for (j=1;j<=ifwfr;j++){
				jj = NX[1]-j+1;
				for (k=1;k<=NX[2];k++){
					absorb_coeff[i][jj][k] *= coeff[j];
				}
			}
		}
		free_vector(coeff,1,ifwfr);
	}

	/* compute coefficients for top grid boundary (z-direction) */
	if (AB_TO){
		fprintf(FP," Width of absorbing frame in z-direction (top)   = %i gridpoints.\n",ifwto);
		fprintf(FP," Percentage of exponential damping (top)    = %5.2f\n",DAMPTO);
		amp = 1.0-DAMPTO/100.0;   /* amplitude at the edge of the numerical grid */
		coeff = vector(1,ifwto);
		a = sqrt(-log(amp)/((ifwto-1)*(ifwto-1)));
		for (k=1;k<=ifwto;k++)
			coeff[k] = exp(-(a*a*(ifwto-k)*(ifwto-k)));
		fprintf(FP," Table of coefficients in z-direction (top) \n # \t coeffto \n");
		for (k=1;k<=ifwto;k++)
			fprintf(FP," %d \t %5.3f \n", k, coeff[k]);

		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
				for (k=1;k<=ifwto;k++){
					absorb_coeff[i][j][k] *= coeff[k];
				}
			}
		}
		free_vector(coeff,1,ifwto);
	}

	/* compute coefficients for bottom grid boundary (z-direction) */
	if (AB_BO){
		fprintf(FP," Width of absorbing frame in z-direction (bottom)= %i gridpoints.\n",ifwbo);
		fprintf(FP," Percentage of exponential damping (bottom) = %5.2f\n",DAMPBO);
		amp = 1.0-DAMPBO/100.0;   /* amplitude at the edge of the numerical grid */
		coeff = vector(1,ifwbo);
		a = sqrt(-log(amp)/((ifwbo-1)*(ifwbo-1)));
		for (k=1;k<=ifwbo;k++)
			coeff[k] = exp(-(a*a*(ifwbo-k)*(ifwbo-k)));
		fprintf(FP," Table of coefficients in z-direction (bottom) \n # \t coeffbo \n");
		for (k=1;k<=ifwbo;k++)
			fprintf(FP," %d \t %5.3f \n", k, coeff[k]);

		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
				for (k=1;k<=ifwbo;k++){
					kk = NX[2]-k+1;
					absorb_coeff[i][j][kk] *= coeff[k];
				}
			}
		}
		free_vector(coeff,1,ifwbo);
	}

}



