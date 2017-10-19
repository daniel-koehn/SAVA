/*------------------------------------------------------------------------
 *   Orthotropic layered model
 *
 *   D. Köhn
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "fd3d.h"

void model(float *** rho, 
           float *** c1111, float *** c1122, float *** c1133, float *** c2222, float *** c2233, float *** c3333, 
	   float *** c2323, float *** c1313, float *** c1212, 
	   float *** taus, float *** taup, float * eta, float * x, float * y, float * z){
	/*--------------------------------------------------------------------------*/
	/* extern variables */
	extern float	DT, *FL, TAU;
	extern int	NX[3], L, MYID;
	extern int	MODEL, MODEL_FORMAT;
	extern char	MODEL_FILE[STRING_SIZE];
	extern FILE	*FP;

	/* local variables */
	float	rhov, muv, piv, vp, vs;
	float	c1111t, c1122t, c1133t, c2222t, c2233t, c3333t;
	float	c2323t, c1313t, c1212t;
	float	ts, tp, sumu, sumpi, ws;
	int	i, j, k, l;
	char	filename[STRING_SIZE], ext[8];


	/* model parameters */
	/* isotropic medium 1 */
//	const float c1111iso1 = 8.0e9,	 c1122iso1 = 2.66e9,	c1133iso1 = 2.66e9; 
//	const float 			 c2222iso1 = 8.0e9,	c2233iso1 = 2.66e9;
//	const float 						c3333iso1 = 8.0e9;
//	const float c2323iso1 = 2.67e9,	 c1313iso1 = 2.67e9,	c1212iso1 = 2.67e9;

	/* orthotropic medium 1 */
//	const float c1111ortho1 = 10.8e9, c1122ortho1 = 2.2e9,	c1133ortho1 = 1.9e9;
//	const float 			  c2222ortho1 = 11.3e9, c2233ortho1 = 1.7e9;
//	const float 						c3333ortho1 = 8.5e9;
//	const float c2323ortho1 = 3.6e9,  c1313ortho1 = 3.9e9,	c1212ortho1 = 4.3e9;

	/* orthotropic medium 2 */
	const float c1111ortho2 = 4.26e9, c1122ortho2 = 2.09e9,	c1133ortho2 = 1.9e9;
	const float 			  c2222ortho2 = 11.3e9,	c2233ortho2 = 0.34e9;
	const float 						c3333ortho2 = 11.05e9;
	const float c2323ortho2 = 1.53e9, c1313ortho2 = 1.65e9,	c1212ortho2 = 2.06e9;

	/* depth of first layer */
	const float Z1 = 1.0;
	/* layer thickness */
	const float D1 = 0.5;

	/* change of material parameter */
	const float perc = 0.2;

	/* parameters for formation 1 */
	const float vp1=2000.0, vs1=1155.0, rho1=2000.0;
	const float c1=5.5e9;

	/*-----------------------------------------------------------------------*/

	fprintf(FP,"\n **Message from function model (printed by PE %d):\n",MYID);
	fprintf(FP," Generating velocity model ... \n");
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
				vp=vp1; vs=vs1;

				/* orthotropic medium */
				rhov   = rho1;

				c1111t = c1111ortho2;
				c1122t = c1122ortho2;
				c1133t = c1133ortho2;
				c2222t = c2222ortho2;
				c2233t = c2233ortho2;
				c3333t = c3333ortho2;

				c2323t = c2323ortho2;
				c1313t = c1313ortho2;
				c1212t = c1212ortho2;
   
				/* alteration of elastic parameters with depth */
				if((z[k]>Z1)&&(z[k]<Z1+D1)){
					c1111t += perc * c1111t;
				}
				if((z[k]>Z1+D1)&&(z[k]<Z1+2*D1)){
					c1122t += perc * c1122t;
				}
				if((z[k]>Z1+2*D1)&&(z[k]<Z1+3*D1)){
					c1133t += perc * c1133t;
				}
				if((z[k]>Z1+3*D1)&&(z[k]<Z1+4*D1)){
					c1212t += perc * c1212t;
				}
				if((z[k]>Z1+4*D1)&&(z[k]<Z1+5*D1)){
					c1313t += perc * c1313t;
				}
				if((z[k]>Z1+5*D1)&&(z[k]<Z1+6*D1)){
					c2222t += perc * c2222t;
				}
				if((z[k]>Z1+6*D1)&&(z[k]<Z1+7*D1)){
					c2233t += perc * c2233t;
				}
				if((z[k]>Z1+7*D1)&&(z[k]<Z1+8*D1)){
					c2323t += perc * c2323t;
				}
				if((z[k]>Z1+8*D1){
					c3333t += perc * c3333t;
				}

				/* end model generation */

				muv = vs*vs*rhov/(1.0+sumu);
				piv = vp*vp*rhov/(1.0+sumpi);

				taus[i][j][k] = ts;
				taup[i][j][k] = tp;
				rho[i][j][k]  = rhov;

				c1111[i][j][k] = c1111t;
				c1122[i][j][k] = c1122t;
				c1133[i][j][k] = c1133t;
				c2222[i][j][k] = c2222t;
				c2233[i][j][k] = c2233t;
				c3333[i][j][k] = c3333t;

				c2323[i][j][k] = c2323t;
				c1313[i][j][k] = c1313t;
				c1212[i][j][k] = c1212t;
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
		sprintf(filename,"%s_c1212%s",MODEL_FILE,ext);

		writemod(filename,c1212,MODEL_FORMAT);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!(MYID))
			mergemod(filename,MODEL_FORMAT);
	}
	if (MODEL & 2){
		sprintf(filename,"%s_c1111%s",MODEL_FILE,ext);

		writemod(filename,c1111,MODEL_FORMAT);
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


