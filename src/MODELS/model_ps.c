/*
 *   Model with layer above halfspace with inclusion
 */

#include "fd.h"
#include "fd3d.h"

void model(float *** rho, float *** pi, float *** mu, float *** taus, float *** taup, float * eta, float * x, float * y, float * z){
	/*--------------------------------------------------------------------------*/
	/* extern variables */
	extern float	DT, *FL, TAU;
	extern int	NX[3], L, MYID;
	extern int	MODEL, MODEL_FORMAT;
	extern char	MODEL_FILE[STRING_SIZE];
	extern FILE	*FP;

	/* local variables */
	float	rhov, muv, piv, vp, vs;
	float	ts, tp, sumu, sumpi, ws;
	int	i, j, k, l;
	char	filename[STRING_SIZE], ext[8];


	/* model parameters */	
	/* parameters for formation 1 */
	const float vp1=2520.0, vs1=1000.0, rho1=2300.0;

	/* parameters for formation 2 */
	const float vp2=2160.0, vs2= 860.0, rho2=1980.0;

	/* depth of layer interface */
	const float Z1 = 20.0;
	const float XI = 160.0, YI = 20.0, ZI = 45.0, RI = 5.0;

	/*-----------------------------------------------------------------------*/

	fprintf(FP,"\n **Message from function model (printed by PE %d):\n",MYID);
	fprintf(FP," Generating P-wave velocity model ... \n");
	fprintf(FP," Generating S-wave velocity model ... \n");
	fprintf(FP," Generating density model ... \n");
	if (L){
		fprintf(FP," Generating QP-model ... \n");
		fprintf(FP," Generating QS-model ... \n");
	}

	/* vector for maxwellbodies */
	tp = TAU;
	ts = TAU;

	sumpi = 0.0;
	sumu  = 0.0;
	for (l=1;l<=L;l++){
		eta[l] = 2.0*PI*DT*FL[l];
		ws     = FL[l]/FL[1];
		sumpi += 1.0/(1.0+ws*ws);
		sumu  += 1.0/(1.0+ws*ws);
	}
	sumpi *= tp;
	sumu  *= ts;


	/* loop over local grid */
	for (i=1;i<=NX[0];i++){
		for (j=1;j<=NX[1];j++){
			for (k=1;k<=NX[2];k++){

				/* model generation */
				if (z[k]<Z1){
					vp=vp1; vs=vs1; rhov=rho1;	/* layer */
				}
				else{
					vp=vp2; vs=vs2; rhov=rho2;	/* halfspace */
				}

				if ((x[i]-XI)*(x[i]-XI)+(y[j]-YI)*(y[j]-YI)+(z[k]-ZI)*(z[k]-ZI) < RI*RI){
					vp=vp1; vs=vs1; rhov=rho1;	/* inclusion */
				}
				/* end model generation */

				muv = vs*vs*rhov/(1.0+sumu);
				piv = vp*vp*rhov/(1.0+sumpi);

				pi[i][j][k]   = piv;
				mu[i][j][k]   = muv;
				rho[i][j][k]  = rhov;
				taus[i][j][k] = ts;
				taup[i][j][k] = tp;
			}
		}
	}

	/* write model to disk */

	/* different data formats of output:
	   MODEL_FORMAT=1:  SU (native byte order)
	   MODEL_FORMAT=2:  ASCII
	   MODEL_FORMAT=3:  BINARY (native byte order) */

	switch(MODEL_FORMAT){
	case 1:
		sprintf(ext,".su");
		break;
	case 2:
		sprintf(ext,".asc");
		break;
	case 3:
		sprintf(ext,".bin");
		break;
	}

	if (MODEL & 4){
		sprintf(filename,"%s_mu%s",MODEL_FILE,ext);

		writemod(filename,mu,MODEL_FORMAT);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!(MYID))
			mergemod(filename,MODEL_FORMAT);
	}
	if (MODEL & 2){
		sprintf(filename,"%s_pi%s",MODEL_FILE,ext);

		writemod(filename,pi,MODEL_FORMAT);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!(MYID))
			mergemod(filename,MODEL_FORMAT);
	}
	if (MODEL & 1){
		sprintf(filename,"%s_rho%s",MODEL_FILE,ext);

		writemod(filename,rho,MODEL_FORMAT);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!(MYID))
			mergemod(filename,MODEL_FORMAT);
	}

}


