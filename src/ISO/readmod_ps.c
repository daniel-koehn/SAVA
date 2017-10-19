/*------------------------------------------------------------------------
 *   Read elastic model properties (vp,vs,density) from files
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "fd3d.h"

void readmod(float *** rho, float *** pi, float *** mu, float *** taus, float *** taup, float * eta){

	/* extern variables */
	extern float	DT, *FL, TAU;
	extern int	NX[3], NXG[3], POS[3], L, MYID;
	extern int	READMOD;
	extern int	MODEL, MODEL_FORMAT;
	extern char	MFILE[STRING_SIZE], MODEL_FILE[STRING_SIZE];
	extern FILE	*FP;

	/* local variables */
	float	rhov, muv, piv, vp, vs;
	float	qp, qs;
	float	ts, tp, sumu, sumpi, ws;
	int	i, j, k, l, ii, jj, kk;
	char	filename[STRING_SIZE], ext[8];
	FILE	*fp_vs, *fp_vp, *fp_rho;
	FILE	*fp_qp=NULL, *fp_qs=NULL;

	/*-----------------------------------------------------------------------*/

	fprintf(FP,"\n **Message from function readmod (printed by PE %d):\n",MYID);

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


	fprintf(FP," Reading model information from model-files...\n");

	sprintf(filename,"%s.vp",MFILE);
	fprintf(FP,"\t P-wave velocities:\n\t %s \n\n",filename);
	fp_vp=fopen(filename,"r");
	if (fp_vp==NULL) error(" Could not open model file for P velocities ! ");

	sprintf(filename,"%s.vs",MFILE);
	fprintf(FP,"\t Shear wave velocities:\n\t %s \n\n",filename);
	fp_vs=fopen(filename,"r");
	if (fp_vs==NULL) error(" Could not open model file for shear velocities ! ");

	sprintf(filename,"%s.rho",MFILE);
	fprintf(FP,"\t Density:\n\t %s \n\n",filename);
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) error(" Could not open model file for densities ! ");

	if (L){
		sprintf(filename,"%s.qp",MFILE);
		fprintf(FP,"\t Qp:\n\t %s \n\n",filename);
		fp_qp=fopen(filename,"r");
		if (fp_qp==NULL) error(" Could not open model file for Qp-values ! ");

		sprintf(filename,"%s.qs",MFILE);
		fprintf(FP,"\t Qs:\n\t %s \n\n",filename);
		fp_qs=fopen(filename,"r");
		if (fp_qs==NULL) error(" Could not open model file for Qs-values ! ");
	}

	/* loop over global grid */
	for (i=1;i<=NXG[0];i++){
		if (!(i%100))
			fprintf(FP,"\n Reading model: step %d of %d \n",i,NXG[0]);
		for (j=1;j<=NXG[1];j++){
			if (READMOD==2){
				/* jump over 240 byte SU header */
				fseek(fp_vp,240,SEEK_CUR);
				fseek(fp_vs,240,SEEK_CUR);
				fseek(fp_rho,240,SEEK_CUR);
				if (L){
					fseek(fp_qp,240,SEEK_CUR);
					fseek(fp_qs,240,SEEK_CUR);
				}
			}
			for (k=1;k<=NXG[2];k++){
				fread(&vp,   sizeof(float), 1, fp_vp);
				fread(&vs,   sizeof(float), 1, fp_vs);
				fread(&rhov, sizeof(float), 1, fp_rho);
				if (L){
					fread(&qp, sizeof(float), 1, fp_qp);
					fread(&qs, sizeof(float), 1, fp_qs);
				}

				muv=vs*vs*rhov/(1.0+sumu);
				piv=vp*vp*rhov/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint 
				is saving model parameters in its local arrays */
				if ((POS[0]==((i-1)/NX[0])) && (POS[1]==((j-1)/NX[1])) && (POS[2]==((k-1)/NX[2]))){
					ii = i-POS[0]*NX[0];
					jj = j-POS[1]*NX[1];
					kk = k-POS[2]*NX[2];

					if (L){
						taus[ii][jj][kk]=2.0/qs;
						taup[ii][jj][kk]=2.0/qp;
					}
					pi[ii][jj][kk]  = piv;
					mu[ii][jj][kk]  = muv;
					rho[ii][jj][kk] = rhov;
				}
			}
		}
	}

	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_rho);
	if (L){
		fclose(fp_qp);
		fclose(fp_qs);
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




