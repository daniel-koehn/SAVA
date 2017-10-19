/*------------------------------------------------------------------------
 *   Exchange FD-Parameters between PEs
 *
 *   T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(void){

	/* extern variables */
	extern int	NXG[3], SRCSIGNAL, MODEL, MODEL_FORMAT, SNAP, SNAP_FORMAT, L;
	extern float	TIME, DT, TS, *FL, TAU;
	extern float	XREC1[3], XREC2[3];
	extern int	SEISMO, NDT, DXREC[3], SEIS_FORMAT, READMOD, READREC;
	extern int	LOG;
	extern int	ABTYPE_LE, ABTYPE_RI, ABTYPE_BA, ABTYPE_FR, ABTYPE_TO, ABTYPE_BO, FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern int	PERIODIC[3];
	extern int	OUT_TIMING;
	extern int	FFID;
	extern float	TSNAP1, TSNAP2, TSNAPINC, REFMOD[3], REFSRC[3], REFREC[3];
	extern float	XSNAPMIN, XSNAPMAX, YSNAPMIN, YSNAPMAX, ZSNAPMIN, ZSNAPMAX;
	extern float	N_SIGMA_LE, SIGMA_LE, N_KAPPA_LE, KAPPA_LE, N_ALPHA_LE, ALPHA_LE;
	extern float	N_SIGMA_RI, SIGMA_RI, N_KAPPA_RI, KAPPA_RI, N_ALPHA_RI, ALPHA_RI;
	extern float	N_SIGMA_BA, SIGMA_BA, N_KAPPA_BA, KAPPA_BA, N_ALPHA_BA, ALPHA_BA;
	extern float	N_SIGMA_FR, SIGMA_FR, N_KAPPA_FR, KAPPA_FR, N_ALPHA_FR, ALPHA_FR;
	extern float	N_SIGMA_TO, SIGMA_TO, N_KAPPA_TO, KAPPA_TO, N_ALPHA_TO, ALPHA_TO;
	extern float	N_SIGMA_BO, SIGMA_BO, N_KAPPA_BO, KAPPA_BO, N_ALPHA_BO, ALPHA_BO;
	extern float	DAMPLE, DAMPRI, DAMPBA, DAMPFR, DAMPTO, DAMPBO;
	extern char	DX_FILE[STRING_SIZE], DY_FILE[STRING_SIZE], DZ_FILE[STRING_SIZE];
	extern char	MFILE[STRING_SIZE], MODEL_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
	extern char	SNAP_FILE[STRING_SIZE], SRC_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char	SEIS_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern int	NPROC, NPROCX[3], MYID, MODEL_IDX[3], IDX[3], CHECKPTREAD, CHECKPTWRITE; 

	/* local variables */
	int	idum[NPAR];
	float	fdum[NPAR];


	if (!(MYID)){ 
		fdum[1]  = TIME;
		fdum[2]  = DT;
		fdum[3]  = TS;
		fdum[4]  = TAU;
		fdum[5]  = TSNAP1;
		fdum[6]  = TSNAP2;
		fdum[7]  = TSNAPINC;
		fdum[8]  = XSNAPMIN;
		fdum[9]  = XSNAPMAX;
		fdum[10] = YSNAPMIN;
		fdum[11] = YSNAPMAX;
		fdum[12] = ZSNAPMIN;
		fdum[13] = ZSNAPMAX;
		fdum[14] = DAMPLE;
		fdum[15] = DAMPRI;
		fdum[16] = DAMPBA;
		fdum[17] = DAMPFR;
		fdum[18] = DAMPTO;
		fdum[19] = DAMPBO;
		fdum[20] = N_SIGMA_LE;
		fdum[21] = SIGMA_LE;
		fdum[22] = N_KAPPA_LE;
		fdum[23] = KAPPA_LE;
		fdum[24] = N_ALPHA_LE;
		fdum[25] = ALPHA_LE;
		fdum[26] = N_SIGMA_RI;
		fdum[27] = SIGMA_RI;
		fdum[28] = N_KAPPA_RI;
		fdum[29] = KAPPA_RI;
		fdum[30] = N_ALPHA_RI;
		fdum[31] = ALPHA_RI;
		fdum[32] = N_SIGMA_BA;
		fdum[33] = SIGMA_BA;
		fdum[34] = N_KAPPA_BA;
		fdum[35] = KAPPA_BA;
		fdum[36] = N_ALPHA_BA;
		fdum[37] = ALPHA_BA;
		fdum[38] = N_SIGMA_FR;
		fdum[39] = SIGMA_FR;
		fdum[40] = N_KAPPA_FR;
		fdum[41] = KAPPA_FR;
		fdum[42] = N_ALPHA_FR;
		fdum[43] = ALPHA_FR;
		fdum[44] = N_SIGMA_TO;
		fdum[45] = SIGMA_TO;
		fdum[46] = N_KAPPA_TO;
		fdum[47] = KAPPA_TO;
		fdum[48] = N_ALPHA_TO;
		fdum[49] = ALPHA_TO;
		fdum[50] = N_SIGMA_BO;
		fdum[51] = SIGMA_BO;
		fdum[52] = N_KAPPA_BO;
		fdum[53] = KAPPA_BO;
		fdum[54] = N_ALPHA_BO;
		fdum[55] = ALPHA_BO;
                                                                                                               
		idum[1]  = LOG;
		idum[2]  = NPROCX[0]*NPROCX[1]*NPROCX[2];
		idum[3]  = ABTYPE_LE;
		idum[4]  = ABTYPE_RI;
		idum[5]  = ABTYPE_BA;
		idum[6]  = ABTYPE_FR;
		idum[7]  = ABTYPE_TO;
		idum[8]  = ABTYPE_BO;
		idum[9]  = SRCSIGNAL;
		idum[10] = READMOD;
		idum[11] = MODEL;
		idum[12] = MODEL_FORMAT;
		idum[13] = L;
		idum[14] = FWLE;
		idum[15] = FWRI;
		idum[16] = FWBA;
		idum[17] = FWFR;
		idum[18] = FWTO;
		idum[19] = FWBO;
		idum[20] = SNAP;
		idum[21] = SNAP_FORMAT;
		idum[22] = SEISMO;
		idum[23] = READREC;
		idum[24] = NDT;
		idum[25] = SEIS_FORMAT;
		idum[26] = OUT_TIMING;
		idum[27] = FFID;
		idum[28] = CHECKPTREAD;
		idum[29] = CHECKPTWRITE;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&idum,NPAR,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdum,NPAR,MPI_FLOAT,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	TIME		= fdum[1];
	DT		= fdum[2];
	TS		= fdum[3];
	TAU		= fdum[4];
	TSNAP1		= fdum[5];
	TSNAP2		= fdum[6];
	TSNAPINC	= fdum[7];
	XSNAPMIN	= fdum[8];
	XSNAPMAX	= fdum[9];
	YSNAPMIN	= fdum[10];
	YSNAPMAX	= fdum[11];
	ZSNAPMIN	= fdum[12];
	ZSNAPMAX	= fdum[13];
	DAMPLE		= fdum[14];
	DAMPRI		= fdum[15];
	DAMPBA		= fdum[16];
	DAMPFR		= fdum[17];
	DAMPTO		= fdum[18];
	DAMPBO		= fdum[19];
	N_SIGMA_LE	= fdum[20];
	SIGMA_LE	= fdum[21];
	N_KAPPA_LE	= fdum[22];
	KAPPA_LE	= fdum[23];
	N_ALPHA_LE	= fdum[24];
	ALPHA_LE	= fdum[25];
	N_SIGMA_RI	= fdum[26];
	SIGMA_RI	= fdum[27];
	N_KAPPA_RI	= fdum[28];
	KAPPA_RI	= fdum[29];
	N_ALPHA_RI	= fdum[30];
	ALPHA_RI	= fdum[31];
	N_SIGMA_BA	= fdum[32];
	SIGMA_BA	= fdum[33];
	N_KAPPA_BA	= fdum[34];
	KAPPA_BA	= fdum[35];
	N_ALPHA_BA	= fdum[36];
	ALPHA_BA	= fdum[37];
	N_SIGMA_FR	= fdum[38];
	SIGMA_FR	= fdum[39];
	N_KAPPA_FR	= fdum[40];
	KAPPA_FR	= fdum[41];
	N_ALPHA_FR	= fdum[42];
	ALPHA_FR	= fdum[43];
	N_SIGMA_TO	= fdum[44];
	SIGMA_TO	= fdum[45];
	N_KAPPA_TO	= fdum[46];
	KAPPA_TO	= fdum[47];
	N_ALPHA_TO	= fdum[48];
	ALPHA_TO	= fdum[49];
	N_SIGMA_BO	= fdum[50];
	SIGMA_BO	= fdum[51];
	N_KAPPA_BO	= fdum[52];
	KAPPA_BO	= fdum[53];
	N_ALPHA_BO	= fdum[54];
	ALPHA_BO	= fdum[55];

	LOG		= idum[1];
	NPROC		= idum[2];
	ABTYPE_LE	= idum[3];
	ABTYPE_RI	= idum[4];
	ABTYPE_BA	= idum[5];
	ABTYPE_FR	= idum[6];
	ABTYPE_TO	= idum[7];
	ABTYPE_BO	= idum[8];
	SRCSIGNAL	= idum[9];
	READMOD		= idum[10];
	MODEL		= idum[11];
	MODEL_FORMAT	= idum[12];
	L		= idum[13];
	FWLE		= idum[14];
	FWRI		= idum[15];
	FWBA		= idum[16];
	FWFR		= idum[17];
	FWTO		= idum[18];
	FWBO		= idum[19];
	SNAP		= idum[20];
	SNAP_FORMAT	= idum[21];
	SEISMO		= idum[22];
	READREC		= idum[23];
	NDT		= idum[24];
	SEIS_FORMAT	= idum[25];
	OUT_TIMING	= idum[26];
	FFID		= idum[27];
	CHECKPTREAD	= idum[28];
	CHECKPTWRITE	= idum[29];

	if (MYID) 
		FL = vector(1,L);

	MPI_Bcast(&NPROCX,   3,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&NXG,      3,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&PERIODIC, 3,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&MODEL_IDX,3,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&IDX,      3,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&DXREC,    3,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(&XREC1,    3,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&XREC2,    3,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&REFMOD,   3,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&REFSRC,   3,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&REFREC,   3,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&FL[1],    L,MPI_FLOAT,0,MPI_COMM_WORLD);

	MPI_Bcast(&DX_FILE,    STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&DY_FILE,    STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&DZ_FILE,    STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SRC_FILE,   STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MFILE,      STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MODEL_FILE, STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SNAP_FILE,  STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&REC_FILE,   STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE,  STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&LOG_FILE,   STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SIGNAL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&CHECKPTFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

}
