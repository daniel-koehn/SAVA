/*------------------------------------------------------------------------
 *   Check parameters in input-file
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void check_par(void){

	/* extern variables */
	extern int	NXG[3], MODEL, SNAP;
//	extern int	SRCSIGNAL, MODEL_FORMAT, SNAP_FORMAT, L;
	extern float	TIME, DT;
//	extern float	*FL, TAU;
//	extern float	XREC1[3], XREC2[3];
	extern int	SEISMO, NDT;
//	extern int	SEIS_FORMAT, READMOD, READREC;
	extern int	DXREC[3];
//	extern int	LOG;
	extern int	ABTYPE_LE, ABTYPE_RI, ABTYPE_BA, ABTYPE_FR, ABTYPE_TO, ABTYPE_BO;
	extern int	FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern int	PERIODIC[3];
//	extern int	OUT_TIMING;
	extern float	N_SIGMA_LE, SIGMA_LE, N_KAPPA_LE, KAPPA_LE, N_ALPHA_LE, ALPHA_LE;
	extern float	N_SIGMA_RI, SIGMA_RI, N_KAPPA_RI, KAPPA_RI, N_ALPHA_RI, ALPHA_RI;
	extern float	N_SIGMA_BA, SIGMA_BA, N_KAPPA_BA, KAPPA_BA, N_ALPHA_BA, ALPHA_BA;
	extern float	N_SIGMA_FR, SIGMA_FR, N_KAPPA_FR, KAPPA_FR, N_ALPHA_FR, ALPHA_FR;
	extern float	N_SIGMA_TO, SIGMA_TO, N_KAPPA_TO, KAPPA_TO, N_ALPHA_TO, ALPHA_TO;
	extern float	N_SIGMA_BO, SIGMA_BO, N_KAPPA_BO, KAPPA_BO, N_ALPHA_BO, ALPHA_BO;
	extern float	DAMPLE, DAMPRI, DAMPBA, DAMPFR, DAMPTO, DAMPBO;
	extern float	TSNAP2;
//	extern float	TSNAP1, TSNAPINC, REFMOD[3], REFSRC[3], REFREC[3];
	extern float	XSNAPMIN, XSNAPMAX, YSNAPMIN, YSNAPMAX, ZSNAPMIN, ZSNAPMAX;
//	extern char	DX_FILE[STRING_SIZE], DY_FILE[STRING_SIZE], DZ_FILE[STRING_SIZE];
//	extern char	MFILE[STRING_SIZE], MODEL_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
//	extern char	SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
//	extern char	SEIS_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern int	NPROCX[3], MODEL_IDX[3], IDX[3], MYID;
//	extern int	CHECKPTREAD, CHECKPTWRITE; 

	/* local variables */
	float	a;


	/* check parameters */
	fprintf(stdout,"\n **Message from check_par (printed by PE %d):\n",MYID);
	fprintf(stdout," Check parameters defined in input file ... \n");

	/*if ((XSNAPMIN < 0)||(XSNAPMAX < 0)){
		fprintf(stdout," XSNAPMIN or XSNAPMAX is less than 0. \n   SNAP is set to 0 (no output of snapshots)!\n");
		SNAP = 0;
	}
	if ((YSNAPMIN < 0)||(YSNAPMAX < 0)){
		fprintf(stdout," YSNAPMIN or YSNAPMAX is less than 0. \n   SNAP is set to 0 (no output of snapshots)!\n");
		SNAP = 0;
	}
	if ((ZSNAPMIN < 0)||(ZSNAPMAX < 0)){
		fprintf(stdout," ZSNAPMIN or ZSNAPMAX is less than 0. \n   SNAP is set to 0 (no output of snapshots)!\n");
		SNAP = 0;
	}*/
	if (XSNAPMIN > XSNAPMAX){
		fprintf(stdout," XSNAPMIN is greater than XSNAPMAX. \n   Values of XSNAPMIN and XSNAPMAX are exchanged!\n");
		a        = XSNAPMIN;
		XSNAPMIN = XSNAPMAX;
		XSNAPMAX = a;
	}
	if (YSNAPMIN > YSNAPMAX){
		fprintf(stdout," YSNAPMIN is greater than YSNAPMAX. \n   Values of YSNAPMIN and YSNAPMAX are exchanged!\n");
		a        = YSNAPMIN;
		YSNAPMIN = YSNAPMAX;
		YSNAPMAX = a;
	}
	if (ZSNAPMIN > ZSNAPMAX){
		fprintf(stdout," ZSNAPMIN is greater than ZSNAPMAX. \n   Values of ZSNAPMIN and ZSNAPMAX are exchanged!\n");
		a        = ZSNAPMIN;
		ZSNAPMIN = ZSNAPMAX;
		ZSNAPMAX = a;
	}

	if (XSNAPMIN == XSNAPMAX){
		fprintf(stdout," XSNAPMIN equals XSNAPMAX. \n   SNAP is set to 0 (no output of snapshots)!\n");
		SNAP = 0;
	}
	if (YSNAPMIN == YSNAPMAX){
		fprintf(stdout," YSNAPMIN equals YSNAPMAX. \n   SNAP is set to 0 (no output of snapshots)!\n");
		SNAP = 0;
	}
	if (ZSNAPMIN == ZSNAPMAX){
		fprintf(stdout," ZSNAPMIN equals ZSNAPMAX. \n   SNAP is set to 0 (no output of snapshots)!\n");
		SNAP = 0;
	}

	if (TSNAP2 > TIME){
		fprintf(stdout," TSNAP2 (output time of last snapshot) is greater than TIME (time of wave propagation). \n   TSNAP2 is set to TIME!\n");
		TSNAP2 = TIME;
	}

	if ((MODEL)&&(NXG[0]%(NPROCX[0]*MODEL_IDX[0])))
		error(" Number of grid points in x-direction (NX) is not divisible by NPROCX[0]*MODEL_IDX[0]! \n   This results in problems when merging model output!\n");
	if ((MODEL)&&(NXG[1]%(NPROCX[1]*MODEL_IDX[1])))
		error(" Number of grid points in y-direction (NY) is not divisible by NPROCX[1]*MODEL_IDX[1]! \n   This results in problems when merging model output!\n");
	if ((MODEL)&&(NXG[2]%(NPROCX[2]*MODEL_IDX[2])))
		error(" Number of grid points in z-direction (NZ) is not divisible by NPROCX[2]*MODEL_IDX[2]! \n   This results in problems when merging model output!\n");	

	if ((SNAP)&&(NXG[0]%(NPROCX[0]*IDX[0])))
		error(" Number of grid points in x-direction (NX) is not divisible by NPROCX[0]*IDX[0]! \n   This results in problems when merging snapshots!\n");
	if ((SNAP)&&(NXG[1]%(NPROCX[1]*IDX[1])))
		error(" Number of grid points in y-direction (NY) is not divisible by NPROCX[1]*IDX[1]! \n   This results in problems when merging snapshots!\n");
	if ((SNAP)&&(NXG[2]%(NPROCX[2]*IDX[2])))
		error(" Number of grid points in z-direction (NZ) is not divisible by NPROCX[2]*IDX[2]! \n   This results in problems when merging snapshots!\n");	

	/* boundary conditions */
	if ((ABTYPE_LE < 0) || (ABTYPE_LE > 2)){
		if (PERIODIC[0]){
			fprintf(stdout," ABTYPE_LE has improper value. \n   ABTYPE_LE is set to 0 (free surface)!\n");
		}
		else{
			fprintf(stdout," ABTYPE_LE has improper value. \n   ABTYPE_LE is set to 0 (periodic BC)!\n");
		}
		ABTYPE_LE = 0;
	}
	if ((ABTYPE_RI < 0) || (ABTYPE_RI > 2)){
		if (PERIODIC[0]){
			fprintf(stdout," ABTYPE_RI has improper value. \n   ABTYPE_RI is set to 0 (free surface)!\n");
		}
		else{
			fprintf(stdout," ABTYPE_RI has improper value. \n   ABTYPE_RI is set to 0 (periodic BC)!\n");
		}
		ABTYPE_RI = 0;
	}
	if ((ABTYPE_BA < 0) || (ABTYPE_BA > 2)){
		if (PERIODIC[1]){
			fprintf(stdout," ABTYPE_BA has improper value. \n   ABTYPE_BA is set to 0 (free surface)!\n");
		}
		else{
			fprintf(stdout," ABTYPE_BA has improper value. \n   ABTYPE_BA is set to 0 (periodic BC)!\n");
		}
		ABTYPE_BA = 0;
	}
	if ((ABTYPE_FR < 0) || (ABTYPE_FR > 2)){
		if (PERIODIC[1]){
			fprintf(stdout," ABTYPE_FR has improper value. \n   ABTYPE_FR is set to 0 (free surface)!\n");
		}
		else{
			fprintf(stdout," ABTYPE_FR has improper value. \n   ABTYPE_FR is set to 0 (periodic BC)!\n");
		}
		ABTYPE_FR = 0;
	}
	if ((ABTYPE_TO < 0) || (ABTYPE_TO > 2)){
		if (PERIODIC[2]){
			fprintf(stdout," ABTYPE_TO has improper value. \n   ABTYPE_TO is set to 0 (free surface)!\n");
		}
		else{
			fprintf(stdout," ABTYPE_TO has improper value. \n   ABTYPE_TO is set to 0 (periodic BC)!\n");
		}
		ABTYPE_TO = 0;
	}
	if ((ABTYPE_BO < 0) || (ABTYPE_BO > 2)){
		if (PERIODIC[2]){
			fprintf(stdout," ABTYPE_BO has improper value. \n   ABTYPE_BO is set to 0 (free surface)!\n");
		}
		else{
			fprintf(stdout," ABTYPE_BO has improper value. \n   ABTYPE_BO is set to 0 (periodic BC)!\n");
		}
		ABTYPE_BO = 0;
	}

	if (ABTYPE_LE == 1){
		if (N_SIGMA_LE < 0.0){
			fprintf(stdout," N_SIGMA_LE is less than 0. \n   N_SIGMA_LE is set to 1!\n");
			N_SIGMA_LE = 1.0;
		}
		if (SIGMA_LE < 0.0){
			fprintf(stdout," SIGMA_LE is less than 0. \n   SIGMA_LE is set to 0!\n");
			SIGMA_LE = 0.0;
		}
		if (N_KAPPA_LE < 0.0){
			fprintf(stdout," N_KAPPA_LE is less than 0. \n   N_KAPPA_LE is set to 1!\n");
			N_KAPPA_LE = 1.0;
		}
		if (KAPPA_LE < 1.0){
			fprintf(stdout," KAPPA_LE is less than 1. \n   KAPPA_LE is set to 1!\n");
			N_KAPPA_LE = 1.0;
		}
		if (N_ALPHA_LE < 0.0){
			fprintf(stdout," N_ALPHA_LE is less than 0. \n   N_ALPHA_LE is set to 1!\n");
			N_ALPHA_LE = 1.0;
		}
		if (ALPHA_LE < 0.0){
			fprintf(stdout," ALPHA_LE is less than 0. \n   ALPHA_LE is set to 0!\n");
			N_ALPHA_LE = 0.0;
		}
	}
	if (ABTYPE_RI == 1){
		if (N_SIGMA_RI < 0.0){
			fprintf(stdout," N_SIGMA_RI is less than 0. \n   N_SIGMA_RI is set to 1!\n");
			N_SIGMA_RI = 1.0;
		}
		if (SIGMA_RI < 0.0){
			fprintf(stdout," SIGMA_RI is less than 0. \n   SIGMA_RI is set to 0!\n");
			SIGMA_RI = 0.0;
		}
		if (N_KAPPA_RI < 0.0){
			fprintf(stdout," N_KAPPA_RI is less than 0. \n   N_KAPPA_RI is set to 1!\n");
			N_KAPPA_RI = 1.0;
		}
		if (KAPPA_RI < 1.0){
			fprintf(stdout," KAPPA_RI is less than 1. \n   KAPPA_RI is set to 1!\n");
			N_KAPPA_RI = 1.0;
		}
		if (N_ALPHA_RI < 0.0){
			fprintf(stdout," N_ALPHA_RI is less than 0. \n   N_ALPHA_RI is set to 1!\n");
			N_ALPHA_RI = 1.0;
		}
		if (ALPHA_RI < 0.0){
			fprintf(stdout," ALPHA_RI is less than 0. \n   ALPHA_RI is set to 0!\n");
			N_ALPHA_RI = 0.0;
		}
	}
	if (ABTYPE_BA == 1){
		if (N_SIGMA_BA < 0.0){
			fprintf(stdout," N_SIGMA_BA is less than 0. \n   N_SIGMA_BA is set to 1!\n");
			N_SIGMA_BA = 1.0;
		}
		if (SIGMA_BA < 0.0){
			fprintf(stdout," SIGMA_BA is less than 0. \n   SIGMA_BA is set to 0!\n");
			SIGMA_BA = 0.0;
		}
		if (N_KAPPA_BA < 0.0){
			fprintf(stdout," N_KAPPA_BA is less than 0. \n   N_KAPPA_BA is set to 1!\n");
			N_KAPPA_BA = 1.0;
		}
		if (KAPPA_BA < 1.0){
			fprintf(stdout," KAPPA_BA is less than 1. \n   KAPPA_BA is set to 1!\n");
			N_KAPPA_BA = 1.0;
		}
		if (N_ALPHA_BA < 0.0){
			fprintf(stdout," N_ALPHA_BA is less than 0. \n   N_ALPHA_BA is set to 1!\n");
			N_ALPHA_BA = 1.0;
		}
		if (ALPHA_BA < 0.0){
			fprintf(stdout," ALPHA_BA is less than 0. \n   ALPHA_BA is set to 0!\n");
			N_ALPHA_BA = 0.0;
		}
	}
	if (ABTYPE_FR == 1){
		if (N_SIGMA_FR < 0.0){
			fprintf(stdout," N_SIGMA_FR is less than 0. \n   N_SIGMA_FR is set to 1!\n");
			N_SIGMA_FR = 1.0;
		}
		if (SIGMA_FR < 0.0){
			fprintf(stdout," SIGMA_FR is less than 0. \n   SIGMA_FR is set to 0!\n");
			SIGMA_FR = 0.0;
		}
		if (N_KAPPA_FR < 0.0){
			fprintf(stdout," N_KAPPA_FR is less than 0. \n   N_KAPPA_FR is set to 1!\n");
			N_KAPPA_FR = 1.0;
		}
		if (KAPPA_FR < 1.0){
			fprintf(stdout," KAPPA_FR is less than 1. \n   KAPPA_FR is set to 1!\n");
			N_KAPPA_FR = 1.0;
		}
		if (N_ALPHA_FR < 0.0){
			fprintf(stdout," N_ALPHA_FR is less than 0. \n   N_ALPHA_FR is set to 1!\n");
			N_ALPHA_FR = 1.0;
		}
		if (ALPHA_FR < 0.0){
			fprintf(stdout," ALPHA_FR is less than 0. \n   ALPHA_FR is set to 0!\n");
			N_ALPHA_FR = 0.0;
		}
	}
	if (ABTYPE_TO == 1){
		if (N_SIGMA_TO < 0.0){
			fprintf(stdout," N_SIGMA_TO is less than 0. \n   N_SIGMA_TO is set to 1!\n");
			N_SIGMA_TO = 1.0;
		}
		if (SIGMA_TO < 0.0){
			fprintf(stdout," SIGMA_TO is less than 0. \n   SIGMA_TO is set to 0!\n");
			N_SIGMA_TO = 0.0;
		}
		if (N_KAPPA_TO < 0.0){
			fprintf(stdout," N_KAPPA_TO is less than 0. \n   N_KAPPA_TO is set to 1!\n");
			N_KAPPA_TO = 1.0;
		}
		if (KAPPA_TO < 1.0){
			fprintf(stdout," KAPPA_TO is less than 1. \n   KAPPA_TO is set to 1!\n");
			KAPPA_TO = 1.0;
		}
		if (N_ALPHA_TO < 0.0){
			fprintf(stdout," N_ALPHA_TO is less than 0. \n   N_ALPHA_TO is set to 1!\n");
			N_ALPHA_TO = 1.0;
		}
		if (ALPHA_TO < 0.0){
			fprintf(stdout," ALPHA_TO is less than 0. \n   ALPHA_TO is set to 0!\n");
			N_ALPHA_TO = 0.0;
		}
	}
	if (ABTYPE_BO == 1){
		if (N_SIGMA_BO < 0.0){
			fprintf(stdout," N_SIGMA_BO is less than 0. \n   N_SIGMA_BO is set to 1!\n");
			N_SIGMA_BO = 1.0;
		}
		if (SIGMA_BO < 0.0){
			fprintf(stdout," SIGMA_BO is less than 0. \n   SIGMA_BO is set to 0!\n");
			N_SIGMA_BO = 0.0;
		}
		if (N_KAPPA_BO < 0.0){
			fprintf(stdout," N_KAPPA_BO is less than 0. \n   N_KAPPA_BO is set to 1!\n");
			N_KAPPA_BO = 1.0;
		}
		if (KAPPA_BO < 1.0){
			fprintf(stdout," KAPPA_BO is less than 1. \n   KAPPA_BO is set to 1!\n");
			KAPPA_BO = 1.0;
		}
		if (N_ALPHA_BO < 0.0){
			fprintf(stdout," N_ALPHA_BO is less than 0. \n   N_ALPHA_BO is set to 1!\n");
			N_ALPHA_BO = 1.0;
		}
		if (ALPHA_BO < 0.0){
			fprintf(stdout," ALPHA_BO is less than 0. \n   ALPHA_BO is set to 0!\n");
			N_ALPHA_BO = 0.0;
		}
	}
	if (FWLE > NXG[0]/NPROCX[0]){
		fprintf(stdout," Absorbing boundary in x-direction (left boundary) is greater than number of grid cells	per CPU in x-direction. \n   FWLE is set to NX!\n");
		FWLE = NXG[0]/NPROCX[0];
	}
	if (FWRI > NXG[0]/NPROCX[0]){
		fprintf(stdout," Absorbing boundary in x-direction (right boundary) is greater than number of grid cells per CPU in x-direction. \n   FWRI is set to NX!\n");
		FWRI = NXG[0]/NPROCX[0];
	}
	if (FWBA > NXG[1]/NPROCX[1]){
		fprintf(stdout," Absorbing boundary in y-direction (back boundary) is greater than number of grid cells	per CPU in y-direction. \n   FWBA is set to NY!\n");
		FWBA = NXG[1]/NPROCX[1];
	}
	if (FWFR > NXG[1]/NPROCX[1]){
		fprintf(stdout," Absorbing boundary in y-direction (front boundary) is greater than number of grid cells per CPU in y-direction. \n   FWFR is set to NY!\n");
		FWFR = NXG[1]/NPROCX[1];
	}
	if (FWTO > NXG[2]/NPROCX[2]){
		fprintf(stdout," Absorbing boundary in z-direction (top boundary) is greater than number of grid cells per CPU in z-direction. \n   FWTO is set to NZ!\n");
		FWTO = NXG[2]/NPROCX[2];
	}
	if (FWBO > NXG[2]/NPROCX[2]){
		fprintf(stdout," Absorbing boundary in z-direction (bottom boundary) is greater than number of grid cells per CPU in z-direction. \n   FWBO is set to NZ!\n");
		FWBO = NXG[2]/NPROCX[2];
	}

	if (((ABTYPE_LE == 1) && (ABTYPE_RI == 2)) || ((ABTYPE_LE == 2) && (ABTYPE_RI == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Left and right boundaries are set to PML (ABTYPE_LE=1 and ABTYPE_RI=1).\n");
		ABTYPE_LE = 1;
		ABTYPE_RI = 1;
	}
	if (((ABTYPE_LE == 1) && (ABTYPE_BA == 2)) || ((ABTYPE_LE == 2) && (ABTYPE_BA == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Left and back boundaries are set to PML (ABTYPE_LE=1 and ABTYPE_BA=1).\n");
		ABTYPE_LE = 1;
		ABTYPE_BA = 1;
	}
	if (((ABTYPE_LE == 1) && (ABTYPE_FR == 2)) || ((ABTYPE_LE == 2) && (ABTYPE_FR == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Left and front boundaries are set to PML (ABTYPE_LE=1 and ABTYPE_FR=1).\n");
		ABTYPE_LE = 1;
		ABTYPE_FR = 1;
	}
	if (((ABTYPE_LE == 1) && (ABTYPE_TO == 2)) || ((ABTYPE_LE == 2) && (ABTYPE_TO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Left and top boundaries are set to PML (ABTYPE_LE=1 and ABTYPE_TO=1).\n");
		ABTYPE_LE = 1;
		ABTYPE_TO = 1;
	}
	if (((ABTYPE_LE == 1) && (ABTYPE_BO == 2)) || ((ABTYPE_LE == 2) && (ABTYPE_BO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Left and bottom boundaries are set to PML (ABTYPE_LE=1 and ABTYPE_BO=1).\n");
		ABTYPE_LE = 1;
		ABTYPE_BO = 1;
	}
	if (((ABTYPE_RI == 1) && (ABTYPE_BA == 2)) || ((ABTYPE_RI == 2) && (ABTYPE_BA == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Right and back boundaries are set to PML (ABTYPE_RI=1 and ABTYPE_BA=1).\n");
		ABTYPE_RI = 1;
		ABTYPE_BA = 1;
	}
	if (((ABTYPE_RI == 1) && (ABTYPE_FR == 2)) || ((ABTYPE_RI == 2) && (ABTYPE_FR == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Right and front boundaries are set to PML (ABTYPE_RI=1 and ABTYPE_FR=1).\n");
		ABTYPE_RI = 1;
		ABTYPE_FR = 1;
	}
	if (((ABTYPE_RI == 1) && (ABTYPE_TO == 2)) || ((ABTYPE_RI == 2) && (ABTYPE_TO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Right and top boundaries are set to PML (ABTYPE_RI=1 and ABTYPE_TO=1).\n");
		ABTYPE_RI = 1;
		ABTYPE_TO = 1;
	}
	if (((ABTYPE_RI == 1) && (ABTYPE_BO == 2)) || ((ABTYPE_RI == 2) && (ABTYPE_BO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Right and bottom boundaries are set to PML (ABTYPE_RI=1 and ABTYPE_BO=1).\n");
		ABTYPE_RI = 1;
		ABTYPE_BO = 1;
	}
	if (((ABTYPE_BA == 1) && (ABTYPE_FR == 2)) || ((ABTYPE_BA == 2) && (ABTYPE_FR == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Back and front boundaries are set to PML (ABTYPE_BA=1 and ABTYPE_FR=1).\n");
		ABTYPE_BA = 1;
		ABTYPE_FR = 1;
	}
	if (((ABTYPE_BA == 1) && (ABTYPE_TO == 2)) || ((ABTYPE_BA == 2) && (ABTYPE_TO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Back and top boundaries are set to PML (ABTYPE_BA=1 and ABTYPE_TO=1).\n");
		ABTYPE_BA = 1;
		ABTYPE_TO = 1;
	}
	if (((ABTYPE_BA == 1) && (ABTYPE_BO == 2)) || ((ABTYPE_BA == 2) && (ABTYPE_BO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Back and bottom boundaries are set to PML (ABTYPE_BA=1 and ABTYPE_BO=1).\n");
		ABTYPE_BA = 1;
		ABTYPE_BO = 1;
	}
	if (((ABTYPE_FR == 1) && (ABTYPE_TO == 2)) || ((ABTYPE_FR == 2) && (ABTYPE_TO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Front and top boundaries are set to PML (ABTYPE_FR=1 and ABTYPE_TO=1).\n");
		ABTYPE_FR = 1;
		ABTYPE_TO = 1;
	}
	if (((ABTYPE_FR == 1) && (ABTYPE_BO == 2)) || ((ABTYPE_FR == 2) && (ABTYPE_BO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Front and bottom boundaries are set to PML (ABTYPE_FR=1 and ABTYPE_BO=1).\n");
		ABTYPE_FR = 1;
		ABTYPE_BO = 1;
	}
	if (((ABTYPE_TO == 1) && (ABTYPE_BO == 2)) || ((ABTYPE_TO == 2) && (ABTYPE_BO == 1))){
		fprintf(stdout," It is not possible to combine PMLs and absorbing frames .\n");
		fprintf(stdout," Top and bottom boundaries are set to PML (ABTYPE_TO=1 and ABTYPE_BO=1).\n");
		ABTYPE_TO = 1;
		ABTYPE_BO = 1;
	}

	if (ABTYPE_LE == 2){
		if (FWLE<25){
			fprintf(stdout," Width (FWLE) of absorbing frame in x-direction (left boundary) should be at least 25 gridpoints.\n");
			fprintf(stdout," You have specified a width of %d gridpoints.\n",FWLE);
			fprintf(stdout," FWLE is set to 25 gridpoints to avoid artificial reflections from grid boundaries ! \n");
			FWLE = 25;
		}
		if ((DAMPLE > 100.0) || (DAMPLE < 0.0)){
			fprintf(stdout," DAMPLE has an improper value. \n   DAMPLE is set to 4 percent !\n");
			DAMPLE = 4.0;
		}
	}
	if (ABTYPE_RI == 2){
		if (FWRI<25){
			fprintf(stdout," Width (FWRI) of absorbing frame in x-direction (right boundary) should be at least 25 gridpoints.\n");
			fprintf(stdout," You have specified a width of %d gridpoints.\n",FWRI);
			fprintf(stdout," FWRI is set to 25 gridpoints to avoid artificial reflections from grid boundaries ! \n");
			FWRI = 25;
		}
		if ((DAMPRI > 100.0) || (DAMPRI < 0.0)){
			fprintf(stdout," DAMPRI has an improper value. \n   DAMPRI is set to 4 percent !\n");
			DAMPRI = 4.0;
		}
	}
	if (ABTYPE_BA == 2){
		if (FWBA<25){
			fprintf(stdout," Width (FWBA) of absorbing frame in y-direction (back boundary) should be at least 25 gridpoints.\n");
			fprintf(stdout," You have specified a width of %d gridpoints.\n",FWBA);
			fprintf(stdout," FWBA is set to 25 gridpoints to avoid artificial reflections from grid boundaries ! \n");
			FWBA = 25;
		}
		if ((DAMPBA > 100.0) || (DAMPBA < 0.0)){
			fprintf(stdout," DAMPBA has an improper value. \n   DAMPBA is set to 4 percent !\n");
			DAMPBA = 4.0;
		}
	}
	if (ABTYPE_FR == 2){
		if (FWFR<25){
			fprintf(stdout," Width (FWFR) of absorbing frame in y-direction (front boundary) should be at least 25 gridpoints.\n");
			fprintf(stdout," You have specified a width of %d gridpoints.\n",FWFR);
			fprintf(stdout," FWFR is set to 25 gridpoints to avoid artificial reflections from grid boundaries ! \n");
			FWFR = 25;
		}
		if ((DAMPFR > 100.0) || (DAMPFR < 0.0)){
			fprintf(stdout," DAMPFR has an improper value. \n   DAMPFR is set to 4 percent !\n");
			DAMPFR = 4.0;
		}
	}
	if (ABTYPE_TO == 2){
		if (FWTO<25){
			fprintf(stdout," Width (FWTO) of absorbing frame in z-direction (top boundary) should be at least 25 gridpoints.\n");
			fprintf(stdout," You have specified a width of %d gridpoints.\n",FWTO);
			fprintf(stdout," FWTO is set to 25 gridpoints to avoid artificial reflections from grid boundaries ! \n");
			FWTO = 25;
		}
		if ((DAMPTO > 100.0) || (DAMPTO < 0.0)){
			fprintf(stdout," DAMPTO has an improper value. \n   DAMPTO is set to 4 percent !\n");
			DAMPTO = 4.0;
		}
	}
	if (ABTYPE_BO == 2){
		if (FWBO<25){
			fprintf(stdout," Width (FWBO) of absorbing frame in z-direction (bottom boundary) should be at least 25 gridpoints.\n");
			fprintf(stdout," You have specified a width of %d gridpoints.\n",FWBO);
			fprintf(stdout," FWBO is set to 25 gridpoints to avoid artificial reflections from grid boundaries ! \n");
			FWBO = 25;
		}
		if ((DAMPBO> 100.0) || (DAMPBO < 0.0)){
			fprintf(stdout," DAMPBO has an improper value. \n   DAMPBO is set to 4 percent !\n");
			DAMPBO = 4.0;
		}
	}
	
	/* seismograms */
	if (((TIME/(DT*NDT)) > 32767)&&(SEISMO)){
		while ((TIME/(DT*NDT)) > 32767)
			NDT <<= 1; /* multiplication by two */
		fprintf(stdout," Too many samples per trace in output seismogramms. \n   NDT is set to %i !\n",NDT);
	}

	/* receiver positions */
	if (DXREC[0] <= 0){
		DXREC[0] = 1;
	}
	if (DXREC[1] <= 0){
		DXREC[1] = 1;
	}
	if (DXREC[2] <= 0){
		DXREC[2] = 1;
	}

	fprintf(stdout,"\n Parameter check finished. \n\n");

}
