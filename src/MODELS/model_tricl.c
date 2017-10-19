/*
 *   Model with isotropic and anisotropic (triclinic) halfspace
 */

#include "fd.h"
#include "fd3d.h"

void model(float *** rho, 
           float *** c1111, float *** c1112, float *** c1113, float *** c1122, float *** c1133, float *** c1123, float *** c1212,
           float *** c1213, float *** c1222, float *** c1223, float *** c1233, float *** c1313, float *** c1322, float *** c1323,
	   float *** c1333, float *** c2222, float *** c2223, float *** c2233, float *** c2323, float *** c2333, float *** c3333,
           float *** taus, float *** taup, float * eta, float * x, float * y, float * z){
	/*--------------------------------------------------------------------------*/
	/* extern variables */
	extern float	DT, *FL, TAU;
	extern int	NX[3], L, MYID;
	extern int	MODEL, MODEL_FORMAT;
	extern char	MODEL_FILE[STRING_SIZE];
	extern FILE	*FP;

	/* local variables */
	float	rhov;
//	float	muv, piv, vp, vs;
	float	c1111t, c1122t, c2222t, c1133t, c2233t, c3333t, c1212t;
	float	c1213t, c1313t, c1223t, c1323t, c2323t, c1112t, c1113t;
	float	c1123t, c1222t, c1322t, c2223t, c1233t, c1333t, c2333t;
	float	ts, tp, sumu, sumpi, ws;
	int	i, j, k, l;
	char	filename[STRING_SIZE], ext[8];


	/* model parameters */
	/* isotropic medium 1 */
	const float c1111iso1 = 8.0e9,	  c1122iso1 = 2.66e9,	c1133iso1 = 2.66e9; 
	const float			  c2222iso1 = 8.0e9,	c2233iso1 = 2.66e9;
	const float						c3333iso1 = 8.0e9;
	const float c2323iso1 = 2.67e9,	  c1323iso1 = 0.0,	c1223iso1 = 0.0;
	const float			  c1313iso1 = 2.67e9,	c1213iso1 = 0.0;
	const float						c1212iso1 = 2.67e9;
	const float c1123iso1 = 0.0,	  c1113iso1 = 0.0,	c1112iso1 = 0.0;	
	const float c2223iso1 = 0.0,	  c1322iso1 = 0.0,	c1222iso1 = 0.0;
	const float c2333iso1 = 0.0,	  c1333iso1 = 0.0,	c1233iso1 = 0.0;

	/* isotropic medium 2 */
//	const float c1111iso2 = 14.606e9, c1122iso2 = 10.006e9,	c1133iso2 = 10.006e9; 
//	const float			  c2222iso2 = 14.606e9,	c2233iso2 = 10.006e9;
//	const float			  			c3333iso2 = 14.606e9;
//	const float c2323iso2 = 2.300e9,  c1323iso2 = 0.0,	c1223iso2 = 0.0;
//	const float			  c1313iso2 = 2.300e9,	c1213iso2 = 0.0;
//	const float						c1212iso2 = 2.300e9;
//	const float c1123iso2 = 0.0,	  c1113iso2 = 0.0,	c1112iso2 = 0.0;	
//	const float c2223iso2 = 0.0,	  c1322iso2 = 0.0,	c1222iso2 = 0.0;
//	const float c2333iso2 = 0.0,	  c1333iso2 = 0.0,	c1233iso2 = 0.0;

	/* isotropic medium 3 */
//	const float c1111iso3 = 9.238e9,  c1122iso3 = 6.309e9,	c1133iso3 = 6.309e9; 
//	const float			  c2222iso3 = 9.238e9,	c2233iso3 = 6.309e9;
//	const float			  			c3333iso3 = 9.238e9;
//	const float c2323iso3 = 1.464e9,  c1323iso3 = 0.0,	c1223iso3 = 0.0;
//	const float			  c1313iso3 = 1.464e9,	c1213iso3 = 0.0;
//	const float						c1212iso3 = 1.464e9;
//	const float c1123iso3 = 0.0,	  c1113iso3 = 0.0,	c1112iso3 = 0.0;	
//	const float c2223iso3 = 0.0,	  c1322iso3 = 0.0,	c1222iso3 = 0.0;
//	const float c2333iso3 = 0.0,	  c1333iso3 = 0.0,	c1233iso3 = 0.0;

	/* orthotropic medium 1 */
//	const float c1111ortho1 = 10.8e9, c1122ortho1 = 2.2e9,	c1133ortho1 = 1.9e9;
//	const float			  c2222ortho1 = 11.3e9, c2233ortho1 = 1.7e9;
//	const float						c3333ortho1 = 8.5e9;
//	const float c2323ortho1 = 3.6e9,  c1323ortho1 = 0.0,	c1223ortho1 = 0.0;
//	const float			  c1313ortho1 = 3.9e9,	c1213ortho1 = 0.0;
//	const float						c1212ortho1 = 4.3e9;
//	const float c1123ortho1 = 0.0,	  c1113ortho1 = 0.0, 	c1112ortho1 = 0.0;
//	const float c2223ortho1 = 0.0,	  c1322ortho1 = 0.0,	c1222ortho1 = 0.0;
//	const float c2333ortho1 = 0.0,	  c1333ortho1 = 0.0,	c1233ortho1 = 0.0;

	/* orthotropic medium 2 */
//	const float c1111ortho2 = 4.26e9, c1122ortho2 = 2.09e9, c1133ortho2 = 1.9e9;
//	const float			  c2222ortho2 = 11.3e9,	c2233ortho2 = 0.34e9;
//	const float						c3333ortho2 = 11.05e9;
//	const float c2323ortho2 = 1.53e9, c1323ortho2 = 0.0,	c1223ortho2 = 0.0;
//	const float			  c1313ortho2 = 1.65e9,	c1213ortho2 = 0.0;
//	const float						c1212ortho2 = 2.06e9;
//	const float c1123ortho2 = 0.0,	  c1113ortho2 = 0.0,	c1112ortho2 = 0.0;
//	const float c2223ortho2 = 0.0,	  c1322ortho2 = 0.0,	c1222ortho2 = 0.0;
//	const float c2333ortho2 = 0.0,	  c1333ortho2 = 0.0,	c1233ortho2 = 0.0;

	/* triclinic medium 1 */
	const float c1111tricl1 = 4.26e9, c1122tricl1 = 2.09e9,	c1133tricl1 = 1.9e9;
	const float			  c2222tricl1 = 11.3e9,	c2233tricl1 = 0.34e9;
	const float						c3333tricl1 = 11.05e9;
	const float c2323tricl1 = 1.53e9, c1323tricl1 = -0.15e9,c1223tricl1 = -0.44e9;
	const float			  c1313tricl1 = 1.65e9,	c1213tricl1 = 0.03e9;
	const float						c1212tricl1 = 2.06e9;
	const float c1123tricl1 = 0.31e9, c1113tricl1 = -0.5e9,	c1112tricl1 = -0.06e9;
	const float c2223tricl1 = -0.24e9,c1322tricl1 = -0.47e9,c1222tricl1 = -0.36e9;
	const float c2333tricl1 = -0.54e9,c1333tricl1 = 0.44e9,	c1233tricl1 = -0.6e9;

	/* depth of layer interface */
	const float Z1 = 1000.0;

	/* parameters for formation 1 */
//	const float vp1=2000.0, vs1=1155.0;
	const float rho1=2000.0;

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
				//vp=vp1; vs=vs1;

				if (z[k] < Z1){		/* isotropic layer */
					rhov   = rho1;

					c1111t = c1111iso1;	/* (C11) */
					c1122t = c1122iso1;	/* (C12) */
					c1133t = c1133iso1;	/* (C13) */
					c2222t = c2222iso1;	/* (C22) */
					c2233t = c2233iso1;	/* (C23) */
					c3333t = c3333iso1;	/* (C33) */

					c1123t = c1123iso1;	/* (C14) */
					c1113t = c1113iso1;	/* (C15) */
					c1112t = c1112iso1;	/* (C16) */
					c2223t = c2223iso1;	/* (C24) */
					c1322t = c1322iso1;	/* (C25) */
					c1222t = c1222iso1;	/* (C26) */
					c2333t = c2333iso1;	/* (C34) */
					c1333t = c1333iso1;	/* (C35) */
					c1233t = c1233iso1;	/* (C36) */

					c2323t = c2323iso1;	/* (C44) */
					c1323t = c1323iso1;	/* (C45) */
					c1223t = c1223iso1;	/* (C46) */
					c1313t = c1313iso1;	/* (C55) */
					c1213t = c1213iso1;	/* (C56) */
					c1212t = c1212iso1;	/* (C66) */
				}
				else{			/* triclinic halfspace */
					rhov   = rho1;

					c1111t = c1111tricl1;	/* (C11) */
					c1122t = c1122tricl1;	/* (C12) */
					c1133t = c1133tricl1;	/* (C13) */
					c2222t = c2222tricl1;	/* (C22) */
					c2233t = c2233tricl1;	/* (C23) */
					c3333t = c3333tricl1;	/* (C33) */

					c1123t = c1123tricl1;	/* (C14) */
					c1113t = c1113tricl1;	/* (C15) */
					c1112t = c1112tricl1;	/* (C16) */
					c2223t = c2223tricl1;	/* (C24) */
					c1322t = c1322tricl1;	/* (C25) */
					c1222t = c1222tricl1;	/* (C26) */
					c2333t = c2333tricl1;	/* (C34) */
					c1333t = c1333tricl1;	/* (C35) */
					c1233t = c1233tricl1;	/* (C36) */

					c2323t = c2323tricl1;	/* (C44) */
					c1323t = c1323tricl1;	/* (C45) */
					c1223t = c1223tricl1;	/* (C46) */
					c1313t = c1313tricl1;	/* (C55) */
					c1213t = c1213tricl1;	/* (C56) */
					c1212t = c1212tricl1;	/* (C66) */
				}

				/* end model generation */

				//muv = vs*vs*rhov/(1.0+sumu);
				//piv = vp*vp*rhov/(1.0+sumpi);

				taus[i][j][k] = ts;
				taup[i][j][k] = tp;
				rho[i][j][k]  = rhov;

				c1111[i][j][k] = c1111t;
				c1122[i][j][k] = c1122t;
				c1133[i][j][k] = c1133t;
				c2222[i][j][k] = c2222t;
				c2233[i][j][k] = c2233t;
				c3333[i][j][k] = c3333t;

				c1123[i][j][k] = c1123t;
				c1113[i][j][k] = c1113t;
				c1112[i][j][k] = c1112t;
				c2223[i][j][k] = c2223t;
				c1322[i][j][k] = c1322t;
				c1222[i][j][k] = c1222t;
				c2333[i][j][k] = c2333t;
				c1333[i][j][k] = c1333t;
				c1233[i][j][k] = c1233t;

				c2323[i][j][k] = c2323t;
				c1323[i][j][k] = c1323t;
				c1223[i][j][k] = c1223t;
				c1313[i][j][k] = c1313t;
				c1213[i][j][k] = c1213t;
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


