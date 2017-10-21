/*------------------------------------------------------------------------
 *  Parallel 3-D Viscoacoustic Finite Difference Modelling 
 *  in Cartesian coordinates
 *  
 *  In case of questions or if you require further information
 *  you are welcome to contact the authors:                       
 *
 *  ----------------------------------------------
 *  Dr. Olaf Hellwig
 *  TU Bergakademie Freiberg
 *  Institut fuer Geophysik
 *  Gustav-Zeuner-Str. 12
 *  D-09596 Freiberg, Germany
 *  Phone: +49 (0)3731 392233
 *  Fax: +49 (0)3731 392636
 *  e-mail: olaf.hellwig[at]geophysik.tu-freiberg.de
 *  http://www.geophysik.tu-freiberg.de
 *  ----------------------------------------------
 *
 *  ------------------------------------------------------------------------
 *  Dr. Daniel Koehn 
 *  Christian-Albrechts Universität Kiel
 *  Institute of Geoscience, Department of Geophysics
 *  Otto-Hahn-Platz 1
 *  D-24098 Kiel, Germany 
 *  Phone: +49 431 880 4566
 *  e-mail:dkoehn[at]geophysik.uni-kiel.de,
 *  http://www.geophysik.uni-kiel.de/~dkoehn
 *  ------------------------------------------------------------------------
 */

#include "fd.h"		/* general include file for FD programs */
#include "fd3d.h"	/* general include file for program FD3D */
#include "fd3dp.h"	/* include file for acoustic version */
#include "fd3diso.h"	/* include file for isotropic version */
#include "AC_struct.h"	/* data structures for AC forward modelling */

void FD_AC(){

/* global variables */
extern float	XS, YS, ZS, TIME, DT, TS, TSNAP1, TSNAP2, TSNAPINC, *FL, TAU, XSNAPMIN, XSNAPMAX, YSNAPMIN, YSNAPMAX, ZSNAPMIN, ZSNAPMAX, XREC1[3], XREC2[3];
extern float	REFMOD[3], REFSRC[3], REFREC[3];
extern float	N_SIGMA_LE, SIGMA_LE, N_KAPPA_LE, KAPPA_LE, N_ALPHA_LE, ALPHA_LE;
extern float	N_SIGMA_RI, SIGMA_RI, N_KAPPA_RI, KAPPA_RI, N_ALPHA_RI, ALPHA_RI;
extern float	N_SIGMA_BA, SIGMA_BA, N_KAPPA_BA, KAPPA_BA, N_ALPHA_BA, ALPHA_BA;
extern float	N_SIGMA_FR, SIGMA_FR, N_KAPPA_FR, KAPPA_FR, N_ALPHA_FR, ALPHA_FR;
extern float	N_SIGMA_TO, SIGMA_TO, N_KAPPA_TO, KAPPA_TO, N_ALPHA_TO, ALPHA_TO;
extern float	N_SIGMA_BO, SIGMA_BO, N_KAPPA_BO, KAPPA_BO, N_ALPHA_BO, ALPHA_BO;
extern float	DAMPLE, DAMPRI, DAMPFR, DAMPBA, DAMPTO, DAMPBO;
extern int	ABTYPE_LE, ABTYPE_RI, ABTYPE_BA, ABTYPE_FR, ABTYPE_TO, ABTYPE_BO;
extern int	FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
extern int	PML, PML_LE, PML_RI, PML_BA, PML_FR, PML_TO, PML_BO;
extern int	AB,  AB_LE,  AB_RI,  AB_BA,  AB_FR,  AB_TO,  AB_BO;
extern int	PERIODIC[3], SEISMO, NDT, DXREC[3], SEIS_FORMAT, READMOD, READREC;
extern int	OUT_DIV_CURL, OUT_ACCEL, OUT_TIMING, NX[3], NT, MODEL, MODEL_FORMAT; 
extern int 	SNAP, SNAP_FORMAT, LOG, SRCSIGNAL, L, NXG[3], MODEL_IDX[3], IDX[3]; 
extern int 	CHECKPTREAD, CHECKPTWRITE, IXSNAPMIN[3], IXSNAPMAX[3], FFID;
extern char	DX_FILE[STRING_SIZE], DY_FILE[STRING_SIZE], DZ_FILE[STRING_SIZE];
extern char	SNAP_FILE[STRING_SIZE], SRC_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
extern char	MFILE[STRING_SIZE], MODEL_FILE[STRING_SIZE], REC_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
extern char	LOG_FILE[STRING_SIZE];
extern char	SEIS_FILE[STRING_SIZE];
extern char	AXIS[];
extern FILE	*FP;

/* Mpi-variables */
extern int		NP, NPROC, NPROCX[3], MYID;
extern int		POS[3], INDEX[6];
extern const int	TAG1,  TAG2,  TAG3,  TAG4,  TAG5,  TAG6;
extern const int	TAG7,  TAG8,  TAG9,  TAG10, TAG11, TAG12;
extern const int	TAG13, TAG14, TAG15, TAG16, TAG17, TAG18;
extern MPI_Comm		COMM_CART;

/* Acquisition geometry */
extern int 		NTR, NTR_LOC, NSRC, NSRC_LOC;
NTR = 0, NTR_LOC = 0, NSRC = 0, NSRC_LOC = 0;

/* variables in FD_AC */
int	ns, nseismograms=0, ndyn, nt, infoout;
int	lsnap, nsnap=0, lsamp=0;

float	memdyn, memmodel, memsrc, memseismograms, membuffer, mempml, memtotal;
float	fac1, fac2, fac3;

char	seisfile[STRING_SIZE];
FILE    *fp;

times.time1 = MPI_Wtime(); 

/* open log-file (each PE is using different file) and print info */
note(stdout);

/* domain decomposition */
initproc();

MPI_Barrier(MPI_COMM_WORLD);

NT    = iround(TIME/DT);	/* number of timesteps */
ns    = iround(NT/NDT);		/* number of samples per trace */
lsnap = iround(TSNAP1/DT);	/* first snapshot at this timestep */


/* output of parameters to log-file or stdout */
if (!(MYID)) 
	write_par(FP);

/* set variables for absorbing boundary condition at left boundary */
if (!(POS[0])){
	switch (ABTYPE_LE){
		case 1 : PML_LE = 1; break;	/* PML */
		case 2 : AB_LE  = 1; break;	/* wavefield damping */
	}
}
/* set variables for absorbing boundary condition at right boundary*/
if (POS[0]==NPROCX[0]-1){
	switch (ABTYPE_RI){
		case 1 : PML_RI = 1; break;	/* PML */
		case 2 : AB_RI  = 1; break;	/* wavefield damping */
	}
}
/* set variables for absorbing boundary condition at back boundary */
if (!(POS[1])){
	switch (ABTYPE_BA){
		case 1 : PML_BA = 1; break;	/* PML */
		case 2 : AB_BA  = 1; break;	/* wavefield damping */
	}
}
/* set variables for absorbing boundary condition at front boundary */
if (POS[1]==NPROCX[1]-1){
	switch (ABTYPE_FR){
		case 1 : PML_FR = 1; break;	/* PML */
		case 2 : AB_FR  = 1; break;	/* wavefield damping */
	}
}
/* set variables for absorbing boundary condition at top boundary */
if (!(POS[2])){
	switch (ABTYPE_TO){
		case 1 : PML_TO = 1; break;	/* PML */
		case 2 : AB_TO  = 1; break;	/* wavefield damping */
	}
}
/* set variables for absorbing boundary condition at bottom boundary */
if (POS[2]==NPROCX[2]-1){
	switch (ABTYPE_BO){
		case 1 : PML_BO = 1; break;	/* PML */
		case 2 : AB_BO  = 1; break;	/* wavefield damping */
	}
}
PML = PML_LE + PML_RI + PML_BA + PML_FR + PML_TO + PML_BO;
AB  = AB_LE  + AB_RI  + AB_BA  + AB_FR  + AB_TO  + AB_BO;

/* allocate memory for x- and z-grid point positions */
geom.x    = vector(0,NX[0]+1);
geom.y    = vector(0,NX[1]+1);
geom.z    = vector(0,NX[2]+1);
geom.xp   = vector(0,NX[0]+1);
geom.yp   = vector(0,NX[1]+1);
geom.zp   = vector(0,NX[2]+1);
geom.xg   = vector(0,NXG[0]+1);
geom.yg   = vector(0,NXG[1]+1);
geom.zg   = vector(0,NXG[2]+1);
geom.xpg  = vector(0,NXG[0]+1);
geom.ypg  = vector(0,NXG[1]+1);
geom.zpg  = vector(0,NXG[2]+1);

/* allocate memory for FD-coefficients for spatial derivatives */
geom.dx   = vector(1,NX[0]);
geom.dxp  = vector(1,NX[0]);
geom.dy   = vector(1,NX[1]);
geom.dyp  = vector(1,NX[1]);
geom.dz   = vector(1,NX[2]);
geom.dzp  = vector(1,NX[2]);

/* read files with grid-spacing */
read_grid(FP,DX_FILE,0,geom.x,geom.xp,geom.xg,geom.xpg);
read_grid(FP,DY_FILE,1,geom.y,geom.yp,geom.yg,geom.ypg);
read_grid(FP,DZ_FILE,2,geom.z,geom.zp,geom.zg,geom.zpg);

/* read source positions and source parameters from SOURCE_FILE and assign positions to local grids*/
acq.srcpos     = sources(geom.xg, geom.yg, geom.zg, geom.xpg, geom.ypg, geom.zpg);
acq.srcpos_loc = splitsrc(acq.srcpos);


/* read receiver positions and assign positions to local grids */
if (SEISMO){
	acq.recpos      = receiver(FP, geom.xg, geom.yg, geom.zg, geom.xpg, geom.ypg, geom.zpg);
        acq.recswitch   = ivector(1,NTR);
	acq.recpos_loc  = splitrec(acq.recpos, acq.recswitch);
	seis.section_fulldata = matrix(1,NTR,1,ns);
}


/* estimate memory requirement of the variables in megabytes */
/* determine number of seismograms */
if (SEISMO){
	if (SEISMO & 8)		/* particle acceleration */
		nseismograms += 3;
	if (SEISMO & 4)		/* div */    
		nseismograms += 1;
	if (SEISMO & 2)		/* pressure */
		nseismograms += 1;
	if (SEISMO & 1)		/* particle velocities */     
		nseismograms += 3;
}

/* compute indices for snapshot limits */
IXSNAPMIN[0] = irnd_to_grid_min(XSNAPMIN,geom.xp,1,NX[0],IDX[0]);
//if (IXSNAPMIN[0] > NX[0])
//	SNAP = 0;
IXSNAPMAX[0] = irnd_to_grid_max(XSNAPMAX,geom.xp,1,NX[0],IDX[0]);
//if (IXSNAPMAX[0] < 1)
//	SNAP = 0;
IXSNAPMIN[1] = irnd_to_grid_min(YSNAPMIN,geom.yp,1,NX[1],IDX[1]);
//if (IXSNAPMIN[1] > NX[1])
//	SNAP = 0;
IXSNAPMAX[1] = irnd_to_grid_max(YSNAPMAX,geom.yp,1,NX[1],IDX[1]);
//if (IXSNAPMAX[1] < 1)
//	SNAP = 0;
IXSNAPMIN[2] = irnd_to_grid_min(ZSNAPMIN,geom.zp,1,NX[2],IDX[2]);
//if (IXSNAPMIN[2] > NX[2])
//	SNAP = 0;
IXSNAPMAX[2] = irnd_to_grid_max(ZSNAPMAX,geom.zp,1,NX[2],IDX[2]);
//if (IXSNAPMAX[2] < 1)
//	SNAP = 0;


/* determine number of dynamic variable fields */
ndyn = 4;
if (L){
	ndyn += 1;
}
if ((SEISMO & 8) || (SNAP & 8)){
	/* particle acceleration */
	ndyn += 3;
	OUT_ACCEL = 1;
}
if ((SEISMO & 4) || (SNAP & 4)){
	/* div */    
	ndyn += 1;
	OUT_DIV_CURL = 1;
}

/* estimate memory requirement for dynamic, static and buffer arrays */
fac1 = (NX[0]+2)*(NX[1]+2)*(NX[2]+2);
fac2 = pow(2.0,-20.0);
fac3 = sizeof(float)*fac2;

memdyn         = (float)(ndyn*fac1)*fac3;							/* dynamic arrays */
memmodel       = (float)((5+(!(!(L)))+(!(!(AB))))*fac1+(4*(NX[0]+NX[1]+NX[2])+2*(NXG[0]+NXG[1]+NXG[2])+24)+2*4*7)*fac3;	/* static arrays (material parameters+axes+timing arrays) */
memsrc         = (float)(NSRC_LOC*(NT+7))*fac3;							/* sources */
memseismograms = (float)(NTR_LOC*(nseismograms*ns+5))*fac3;					/* seismograms */
membuffer      = (float)(4*(NX[0]*NX[1]+NX[1]*NX[2]+NX[0]*NX[2]))*fac3;				/* buffer arrays for exchange */
mempml         = (float)(2*((3+NX[0]*NX[1])*(PML_TO*FWTO+PML_BO*FWBO)+(3+NX[1]*NX[2])*(PML_LE*FWLE+PML_RI*FWRI)+(3+NX[0]*NX[2])*(PML_BA*FWBA+PML_FR*FWFR)))*fac3;	/* PML arrays */
memtotal       = memdyn+memmodel+memsrc+memseismograms+membuffer+mempml;


fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
fprintf(FP," Size of local grids: NX =%d \t NY =%d \t NZ =%d\n",  NX[0], NX[1], NX[2]);
fprintf(FP," Size of global grid: NXG=%d \t NYG=%d \t NZG=%d\n\n",NXG[0],NXG[1],NXG[2]);
fprintf(FP," Each process is now trying to allocate memory for:\n");
fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
fprintf(FP," Sources: \t\t\t %6.2f MB\n", memsrc);
fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
if (PML){
	fprintf(FP," PML variables: \t\t %6.2f MB\n", mempml);
}
fprintf(FP," ------------------------------------------------ \n");
fprintf(FP," Max. total memory required per CPU:  %6.2f MB.\n\n", memtotal);


fprintf(FP,"\n **Message from main (printed by PE %d): \n",MYID);
fprintf(FP," Start memory allocation ...\n");

/* allocation for timing arrays used for performance analysis */
times.time_update = dvector(1,7);
times.time_sum    = dvector(1,7);
times.time_avg    = dvector(1,7);
times.time_std    = dvector(1,7);


/* memory allocation for dynamic (wavefield) arrays */
wave.p =        f3tensor(0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);
wave.v = vector3d_tensor(0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);

if (L){ /* memory variables for VISCO-elastic simulation */
	wave.rpp = f3tensor(1,NX[0],1,NX[1],1,NX[2]);
}

/* memory allocation for static (model) arrays */
mat.rho  = f3tensor(0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);
mat.pi   = f3tensor(0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);
mat.taup = f3tensor(0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);

/* memory allocation for averaged material properties */
mat.rhoijpkp = f3tensor(1,NX[0],1,NX[1],1,NX[2]);
mat.rhoipjkp = f3tensor(1,NX[0],1,NX[1],1,NX[2]);
mat.rhoipjpk = f3tensor(1,NX[0],1,NX[1],1,NX[2]);

mat.eta = vector(1,L);

/* memory allocation for absorbing boundaries */
if (AB) pmls.absorb_coeff = f3tensor(1,NX[0],1,NX[1],1,NX[2]);

/* memory allocation for PML arrays */
if (PML_LE){
	/* PML-profile */
	pmls.pmlle	= pml_vector(1,FWLE);

	/* left boundary */
	pmls.Pp_x_l	= f3tensor(1,FWLE,1,NX[1],1,NX[2]);
	pmls.Pvx_x_l	= f3tensor(1,FWLE,1,NX[1],1,NX[2]);
}
if (PML_RI){
	/* PML-profile */
	pmls.pmlri	= pml_vector(1,FWRI);

	/* right boundary */
	pmls.Pp_x_r	= f3tensor(NX[0]-FWRI+1,NX[0],1,NX[1],1,NX[2]);
	pmls.Pvx_x_r	= f3tensor(NX[0]-FWRI+1,NX[0],1,NX[1],1,NX[2]);
}
if (PML_BA){
	/* PML-profile */
	pmls.pmlba	= pml_vector(1,FWBA);

	/* back boundary */
	pmls.Pp_y_b	= f3tensor(1,NX[0],1,FWBA,1,NX[2]);
	pmls.Pvy_y_b	= f3tensor(1,NX[0],1,FWBA,1,NX[2]);
}
if (PML_FR){
	/* PML-profile */
	pmls.pmlfr	= pml_vector(1,FWFR);

	/* front boundary */
	pmls.Pp_y_f	= f3tensor(1,NX[0],NX[1]-FWFR+1,NX[1],1,NX[2]);
	pmls.Pvy_y_f	= f3tensor(1,NX[0],NX[1]-FWFR+1,NX[1],1,NX[2]);
}
if (PML_TO) {
	/* PML-profile */
	pmls.pmlto	= pml_vector(1,FWTO);

	/* top boundary */
	pmls.Pp_z_t	= f3tensor(1,NX[0],1,NX[1],1,FWTO);
	pmls.Pvz_z_t	= f3tensor(1,NX[0],1,NX[1],1,FWTO);
}
if (PML_BO){
	/* PML-profile */
	pmls.pmlbo	= pml_vector(1,FWBO);

	/* bottom boundary */
	pmls.Pp_z_b	= f3tensor(1,NX[0],1,NX[1],NX[2]-FWBO+1,NX[2]);
	pmls.Pvz_z_b	= f3tensor(1,NX[0],1,NX[1],NX[2]-FWBO+1,NX[2]);
}


/* memory allocation for buffer arrays in which the wavefield
	   information which is exchanged between neighbouring PEs is stored */
mpi.sbuf_lef_to_rig = matrix(1,NX[1],1,NX[2]);
mpi.sbuf_rig_to_lef = matrix(1,NX[1],1,NX[2]);
mpi.sbuf_bac_to_fro = matrix(1,NX[0],1,NX[2]);
mpi.sbuf_fro_to_bac = matrix(1,NX[0],1,NX[2]);
mpi.sbuf_top_to_bot = matrix(1,NX[0],1,NX[1]);
mpi.sbuf_bot_to_top = matrix(1,NX[0],1,NX[1]);
mpi.rbuf_lef_to_rig = matrix(1,NX[1],1,NX[2]);
mpi.rbuf_rig_to_lef = matrix(1,NX[1],1,NX[2]);
mpi.rbuf_bac_to_fro = matrix(1,NX[0],1,NX[2]);
mpi.rbuf_fro_to_bac = matrix(1,NX[0],1,NX[2]);
mpi.rbuf_top_to_bot = matrix(1,NX[0],1,NX[1]);
mpi.rbuf_bot_to_top = matrix(1,NX[0],1,NX[1]);


/* memory allocation for dynamic variables (particle acceleration, div and curl)*/
if (OUT_ACCEL){
	/* particle acceleration */
	wave.a = vector3d_tensor(1,NX[0],1,NX[1],1,NX[2]);
}
if (OUT_DIV_CURL){
	/* div */
	wave.diverg = f3tensor(1,NX[0],1,NX[1],1,NX[2]);
}


/* memory allocation for seismogram arrays */
if (SEISMO && (NTR_LOC>0)){      

	if (SEISMO & 8){
		/* particle acceleration */
		seis.sectionax = matrix(1,NTR_LOC,1,ns);
		seis.sectionay = matrix(1,NTR_LOC,1,ns);
		seis.sectionaz = matrix(1,NTR_LOC,1,ns);
	}
	if (SEISMO & 4){
		/* div */
		seis.sectiondiv = matrix(1,NTR_LOC,1,ns);
	}
	if (SEISMO & 2){
		/* pressure */
		seis.sectionp = matrix(1,NTR_LOC,1,ns);
	}
	if (SEISMO & 1){
		/* particle velocities */     
		seis.sectionvx = matrix(1,NTR_LOC,1,ns);
		seis.sectionvy = matrix(1,NTR_LOC,1,ns);
		seis.sectionvz = matrix(1,NTR_LOC,1,ns);
	}
}

fprintf(FP," ... memory allocation for PE %d was successful.\n\n", MYID);


/* calculate wavelet for each source point */
if (NSRC_LOC) acq.signals = wavelet(acq.srcpos_loc);


/* output source signal e.g. for cross-correlation or comparison with analytical solutions */
if ((SRCSIGNAL<6)&&(NSRC_LOC)){
	char  source_signal_file[STRING_SIZE];
	sprintf(source_signal_file,"%s_source_signal.%d.su", MFILE, MYID);
	fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
	output_source_signal(fopen(source_signal_file,"w"), acq.signals, acq.srcpos_loc, NT, 1, geom.xg, geom.yg, geom.zg, geom.xpg, geom.ypg, geom.zpg);
}


MPI_Barrier(MPI_COMM_WORLD);


/* create model grids */
if (READMOD)
	readmod(mat.rho,mat.pi,mat.taup,mat.eta);
else
	model(mat.rho,mat.pi,mat.taup,mat.eta,geom.xp,geom.yp,geom.zp);


/* check if the FD run will be stable and free of numerical dispersion */
checkfd(FP,mat.rho,mat.pi,mat.taup,mat.eta,geom.x,geom.y,geom.z);

/* calculate 3-D array for exponential wavefield(!) damping in the absorbing frame */
if (AB) absorb(pmls.absorb_coeff);

/* calculate 1-D profile for PML damping in the absorbing frame */
if (PML) pml_profile(pmls.pmlle, pmls.pmlri, pmls.pmlba, pmls.pmlfr, pmls.pmlto, pmls.pmlbo, geom.xg, geom.yg, geom.zg, geom.xpg, geom.ypg, geom.zpg);

/* copy material parameters from bottom / front / right boundary to 
top / back / left boundary for averaging at periodic boundaries */
fprintf(FP,"\n **Message from matcopy (written by PE %d):\n",MYID);
fprintf(FP," Copy material properties from right to left / front to back / bottom to top boundary ... \n");

matcopy(mat.rho);

fprintf(FP," Copy material properties finished \n\n");

/* averaging of material properties */
fprintf(FP,"\n**Message from av_mat (written by PE %d):\n",MYID);

av_mat(mat.rho,mat.rhoijpkp,2,'a',geom.x,geom.xp,geom.y,geom.yp,geom.z,geom.zp);
av_mat(mat.rho,mat.rhoipjkp,4,'a',geom.x,geom.xp,geom.y,geom.yp,geom.z,geom.zp);
av_mat(mat.rho,mat.rhoipjpk,6,'a',geom.x,geom.xp,geom.y,geom.yp,geom.z,geom.zp);

fprintf(FP," Averaging of material properties finished \n\n");


/* compute FD-coefficients for spatial derivatives and wavefield averaging */
fd_coeff(1, NX[0], 1, NX[1], 1, NX[2], geom.x, geom.y, geom.z, geom.xp, geom.yp, geom.zp, geom.dx, geom.dxp, geom.dy, geom.dyp, geom.dz, geom.dzp, 
	mat.pi, mat.rhoijpkp, mat.rhoipjkp, mat.rhoipjpk);

/* initialize persistent MPI communication */
/* exchange_v: exchange one layer left-to-right, one layer back-to-front and one layer top-to-bottom */
MPI_Recv_init(&mpi.rbuf_lef_to_rig[1][1],NX[1]*NX[2],MPI_FLOAT,INDEX[1],TAG1,COMM_CART,&mpi.request_v[0]);
MPI_Recv_init(&mpi.rbuf_bac_to_fro[1][1],NX[0]*NX[2],MPI_FLOAT,INDEX[3],TAG3,COMM_CART,&mpi.request_v[1]);
MPI_Recv_init(&mpi.rbuf_top_to_bot[1][1],NX[0]*NX[1],MPI_FLOAT,INDEX[5],TAG5,COMM_CART,&mpi.request_v[2]);
MPI_Send_init(&mpi.sbuf_lef_to_rig[1][1],NX[1]*NX[2],MPI_FLOAT,INDEX[0],TAG1,COMM_CART,&mpi.request_v[3]);
MPI_Send_init(&mpi.sbuf_bac_to_fro[1][1],NX[0]*NX[2],MPI_FLOAT,INDEX[2],TAG3,COMM_CART,&mpi.request_v[4]);
MPI_Send_init(&mpi.sbuf_top_to_bot[1][1],NX[0]*NX[1],MPI_FLOAT,INDEX[4],TAG5,COMM_CART,&mpi.request_v[5]);
	
/* exchange_s: exchange one layer right-to-left, one layer front-to-back and one layer bottom-to-top */
MPI_Recv_init(&mpi.rbuf_rig_to_lef[1][1],NX[1]*NX[2],MPI_FLOAT,INDEX[0],TAG2,COMM_CART,&mpi.request_s[0]);
MPI_Recv_init(&mpi.rbuf_fro_to_bac[1][1],NX[0]*NX[2],MPI_FLOAT,INDEX[2],TAG4,COMM_CART,&mpi.request_s[1]);
MPI_Recv_init(&mpi.rbuf_bot_to_top[1][1],NX[0]*NX[1],MPI_FLOAT,INDEX[4],TAG6,COMM_CART,&mpi.request_s[2]);
MPI_Send_init(&mpi.sbuf_rig_to_lef[1][1],NX[1]*NX[2],MPI_FLOAT,INDEX[1],TAG2,COMM_CART,&mpi.request_s[3]);
MPI_Send_init(&mpi.sbuf_fro_to_bac[1][1],NX[0]*NX[2],MPI_FLOAT,INDEX[3],TAG4,COMM_CART,&mpi.request_s[4]);
MPI_Send_init(&mpi.sbuf_bot_to_top[1][1],NX[0]*NX[1],MPI_FLOAT,INDEX[5],TAG6,COMM_CART,&mpi.request_s[5]);

/* Synchronizing PEs before starting loop of time steps */
MPI_Barrier(MPI_COMM_WORLD);

/* Acoustic forward modelling */
forward_shot_AC(&wave, &pmls, &mat, &geom, &mpi, &seis, &acq, &times, ns);

/* output timing information (real times for update and exchange) */
timing(times.time_avg, times.time_std, 1);

/* free request handles for persistent MPI communication */
MPI_Request_free(&mpi.request_v[0]);
MPI_Request_free(&mpi.request_v[1]);
MPI_Request_free(&mpi.request_v[2]);
MPI_Request_free(&mpi.request_v[3]);
MPI_Request_free(&mpi.request_v[4]);
MPI_Request_free(&mpi.request_v[5]);

MPI_Request_free(&mpi.request_s[0]);
MPI_Request_free(&mpi.request_s[1]);
MPI_Request_free(&mpi.request_s[2]);
MPI_Request_free(&mpi.request_s[3]);
MPI_Request_free(&mpi.request_s[4]);
MPI_Request_free(&mpi.request_s[5]);


fprintf(FP,"\n **Message from main (printed by PE %d): \n",MYID);
fprintf(FP," Deallocation of memory ...\n");

/* deallocation of memory */
free_vector(geom.x,  0,NX[0]+1);
free_vector(geom.xp, 0,NX[0]+1);
free_vector(geom.y,  0,NX[1]+1);
free_vector(geom.yp, 0,NX[1]+1);
free_vector(geom.z,  0,NX[2]+1);
free_vector(geom.zp, 0,NX[2]+1);
free_vector(geom.xg, 0,NXG[0]+1);
free_vector(geom.xpg,0,NXG[0]+1);
free_vector(geom.yg, 0,NXG[1]+1);
free_vector(geom.ypg,0,NXG[1]+1);
free_vector(geom.zg, 0,NXG[2]+1);
free_vector(geom.zpg,0,NXG[2]+1);

/* free memory for FD-coefficients for spatial derivatives */
free_vector(geom.dx, 1,NX[0]);
free_vector(geom.dxp,1,NX[0]);
free_vector(geom.dy, 1,NX[1]);
free_vector(geom.dyp,1,NX[1]);
free_vector(geom.dz, 1,NX[2]);
free_vector(geom.dzp,1,NX[2]);

/* free timing arrays */
free_dvector(times.time_update,1,7);
free_dvector(times.time_sum,   1,7);
free_dvector(times.time_avg,   1,7);
free_dvector(times.time_std,   1,7);

/* free pressure arrays */
free_f3tensor(wave.p,0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);

/* free particle velocity arrays */
free_vector3d_tensor(wave.v,0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);

/* free stress memory arrays */
if (L){
	free_f3tensor(wave.rpp,1,NX[0],1,NX[1],1,NX[2]);
}

/* free arrays for material parameters */
free_f3tensor(mat.rho, 0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);
free_f3tensor(mat.pi,  0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);
free_f3tensor(mat.taup,0,NX[0]+1,0,NX[1]+1,0,NX[2]+1);

free_f3tensor(mat.rhoijpkp,1,NX[0],1,NX[1],1,NX[2]);
free_f3tensor(mat.rhoipjkp,1,NX[0],1,NX[1],1,NX[2]);
free_f3tensor(mat.rhoipjpk,1,NX[0],1,NX[1],1,NX[2]);

free_vector(mat.eta,1,L);

if (AB) free_f3tensor(pmls.absorb_coeff,1,NX[0],1,NX[1],1,NX[2]);

/* free variables for PML */
if (PML_LE){
	/* PML-profile */
	free_pml_vector(pmls.pmlle,1,FWLE);

	/* left boundary */
	free_f3tensor(pmls.Pp_x_l, 1,FWLE,1,NX[1],1,NX[2]);
	free_f3tensor(pmls.Pvx_x_l,1,FWLE,1,NX[1],1,NX[2]);
}
if (PML_RI){
	/* PML-profile */
	free_pml_vector(pmls.pmlri,1,FWRI);

	/* right boundary */
	free_f3tensor(pmls.Pp_x_r, NX[0]-FWRI+1,NX[0],1,NX[1],1,NX[2]);
	free_f3tensor(pmls.Pvx_x_r,NX[0]-FWRI+1,NX[0],1,NX[1],1,NX[2]);
}
if (PML_BA){
	/* PML-profile */
	free_pml_vector(pmls.pmlba,1,FWBA);

	/* back boundary */
	free_f3tensor(pmls.Pp_y_b,1, NX[0],1,FWBA,1,NX[2]);
	free_f3tensor(pmls.Pvy_y_b,1,NX[0],1,FWBA,1,NX[2]);
}
if (PML_FR){
	/* PML-profile */
	free_pml_vector(pmls.pmlfr,1,FWFR);

	/* front boundary */
	free_f3tensor(pmls.Pp_y_f,1, NX[0],NX[1]-FWFR+1,NX[1],1,NX[2]);
	free_f3tensor(pmls.Pvy_y_f,1,NX[0],NX[1]-FWFR+1,NX[1],1,NX[2]);
}
if (PML_TO){
	/* PML-profile */
	free_pml_vector(pmls.pmlto,1,FWTO);

	/* top boundary */
	free_f3tensor(pmls.Pp_z_t,1, NX[0],1,NX[1],1,FWTO);
	free_f3tensor(pmls.Pvz_z_t,1,NX[0],1,NX[1],1,FWTO);
}
if (PML_BO){
	/* PML-profile */
	free_pml_vector(pmls.pmlbo,1,FWBO);

	/* bottom boundary */
	free_f3tensor(pmls.Pp_z_b, 1,NX[0],1,NX[1],NX[2]-FWBO+1,NX[2]);
	free_f3tensor(pmls.Pvz_z_b,1,NX[0],1,NX[1],NX[2]-FWBO+1,NX[2]);
}

/* free memory for buffer arrays */
free_matrix(mpi.sbuf_lef_to_rig,1,NX[1],1,NX[2]);
free_matrix(mpi.sbuf_rig_to_lef,1,NX[1],1,NX[2]);
free_matrix(mpi.sbuf_bac_to_fro,1,NX[0],1,NX[2]);
free_matrix(mpi.sbuf_fro_to_bac,1,NX[0],1,NX[2]);
free_matrix(mpi.sbuf_top_to_bot,1,NX[0],1,NX[1]);
free_matrix(mpi.sbuf_bot_to_top,1,NX[0],1,NX[1]);
free_matrix(mpi.rbuf_lef_to_rig,1,NX[1],1,NX[2]);
free_matrix(mpi.rbuf_rig_to_lef,1,NX[1],1,NX[2]);
free_matrix(mpi.rbuf_bac_to_fro,1,NX[0],1,NX[2]);
free_matrix(mpi.rbuf_fro_to_bac,1,NX[0],1,NX[2]);
free_matrix(mpi.rbuf_top_to_bot,1,NX[0],1,NX[1]);
free_matrix(mpi.rbuf_bot_to_top,1,NX[0],1,NX[1]);

/* free memory for dynamic variables (particle acceleration, div and curl)*/
if (OUT_ACCEL){
	/* particle acceleration */
	free_vector3d_tensor(wave.a,1,NX[0],1,NX[1],1,NX[2]);
}
if (OUT_DIV_CURL){
	/* div */
	free_f3tensor(wave.diverg,1,NX[0],1,NX[1],1,NX[2]);
}

if (SEISMO && (NTR_LOC>0)){
	free_imatrix(acq.recpos_loc,1,NTR,1,7);

	if (SEISMO & 8){
		/* particle acceleration */
		free_matrix(seis.sectionax,1,NTR_LOC,1,ns);
		free_matrix(seis.sectionay,1,NTR_LOC,1,ns);
		free_matrix(seis.sectionaz,1,NTR_LOC,1,ns);
	}
	if (SEISMO & 4){
		/* div */
		free_matrix(seis.sectiondiv,1,NTR_LOC,1,ns);
	}
	if (SEISMO & 2){
		/* pressure */
		free_matrix(seis.sectionp,1,NTR_LOC,1,ns);
	}
	if (SEISMO & 1){
		/* particle velocities */
		free_matrix(seis.sectionvx,1,NTR_LOC,1,ns);
		free_matrix(seis.sectionvy,1,NTR_LOC,1,ns);
		free_matrix(seis.sectionvz,1,NTR_LOC,1,ns);
	}
}

/* free memory for section_fulldata, recswitch and global receiver positions */
if(SEISMO){
   free_matrix(seis.section_fulldata,1,NTR,1,ns);
   free_ivector(acq.recswitch,1,NTR);
   free_imatrix(acq.recpos,1,NTR,1,7);
}

/* free memory for local source positions and source signals */
if (NSRC_LOC){
	free_matrix(acq.signals,   1,NSRC_LOC,1,NT);
	free_matrix(acq.srcpos_loc,1,NSRC_LOC,1,9);
}

/* free memory for global source positions */
free_matrix(acq.srcpos,1,NSRC,1,8);

fprintf(FP," Deallocation of memory finished \n\n");


MPI_Barrier(MPI_COMM_WORLD);

fprintf(FP,"\n **Message from main (printed by PE %d): \n",MYID);
times.time4 = MPI_Wtime();
fprintf(FP," Total real time of program: %4.2f seconds.\n\n",times.time4-times.time1);


fclose(FP);

}
