/*------------------------------------------------------------------------
 *   globvar.h - global variables of viscoelastic FD programs   
 *
 *  ----------------------------------------------------------------------*/

/* definition of global variables used in the finite difference program FD2D */
/* Througout this program  uppercase letters are used for the names of the global variables */

float	XS, YS, ZS, TIME, DT, TS;
float	TSNAP1, TSNAP2, TSNAPINC, *FL, TAU;
float	XSNAPMIN, XSNAPMAX, YSNAPMIN, YSNAPMAX, ZSNAPMIN, ZSNAPMAX;
float	XREC1[3], XREC2[3];
float	REFMOD[3]={0.0, 0.0, 0.0}, REFSRC[3]={0.0, 0.0, 0.0}, REFREC[3]={0.0, 0.0, 0.0};
float	N_SIGMA_LE, SIGMA_LE, N_KAPPA_LE, KAPPA_LE, N_ALPHA_LE, ALPHA_LE;
float	N_SIGMA_RI, SIGMA_RI, N_KAPPA_RI, KAPPA_RI, N_ALPHA_RI, ALPHA_RI;
float	N_SIGMA_BA, SIGMA_BA, N_KAPPA_BA, KAPPA_BA, N_ALPHA_BA, ALPHA_BA;
float	N_SIGMA_FR, SIGMA_FR, N_KAPPA_FR, KAPPA_FR, N_ALPHA_FR, ALPHA_FR;
float	N_SIGMA_TO, SIGMA_TO, N_KAPPA_TO, KAPPA_TO, N_ALPHA_TO, ALPHA_TO;
float	N_SIGMA_BO, SIGMA_BO, N_KAPPA_BO, KAPPA_BO, N_ALPHA_BO, ALPHA_BO;
float	DAMPLE, DAMPRI, DAMPFR, DAMPBA, DAMPTO, DAMPBO;
int	ABTYPE_LE, ABTYPE_RI, ABTYPE_BA, ABTYPE_FR, ABTYPE_TO, ABTYPE_BO;
int	FWLE=0, FWRI=0, FWBA=0, FWFR=0, FWTO=0, FWBO=0;
int	PML=0, PML_LE=0, PML_RI=0, PML_BA=0, PML_FR=0, PML_TO=0, PML_BO=0;
int	AB=0,  AB_LE=0,  AB_RI=0,  AB_BA=0,  AB_FR=0,  AB_TO=0,  AB_BO=0;
int	PERIODIC[3]={0, 0, 0};
int	SEISMO, NDT, DXREC[3], SEIS_FORMAT, READMOD, READREC;
int	OUT_DIV_CURL=0, OUT_ACCEL=0;
int	OUT_TIMING;
int	NX[3], NT, MODEL, MODEL_FORMAT, SNAP, SNAP_FORMAT, LOG, SRCSIGNAL;
int	L, NXG[3], MODEL_IDX[3], IDX[3], CHECKPTREAD, CHECKPTWRITE;
int	IXSNAPMIN[3], IXSNAPMAX[3];
int	FFID;
char	DX_FILE[STRING_SIZE], DY_FILE[STRING_SIZE], DZ_FILE[STRING_SIZE];
char	SNAP_FILE[STRING_SIZE], SRC_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
char	MFILE[STRING_SIZE], MODEL_FILE[STRING_SIZE], REC_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
char	LOG_FILE[STRING_SIZE];
char	SEIS_FILE[STRING_SIZE];
char	AXIS[]="xyz";
FILE	*FP;

/* Mpi-variables */
int		NP, NPROC, NPROCX[3], MYID;
int		POS[3], INDEX[6];
const int	TAG1 =1,  TAG2 =2,  TAG3 =3,  TAG4 =4,  TAG5 =5,  TAG6 =6;
const int	TAG7 =7,  TAG8 =8,  TAG9 =9,  TAG10=10, TAG11=11, TAG12=12;
const int	TAG13=13, TAG14=14, TAG15=15, TAG16=16, TAG17=17, TAG18=18;
MPI_Comm	COMM_CART;

/* Acquisition geometry */
int NTR, NTR_LOC, NSRC, NSRC_LOC;
