/*------------------------------------------------------------------------
 *   Write FD-Parameters to stdout
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void write_par(FILE *fp){

	/*extern variables */
	extern int	NXG[3], NT, SRCSIGNAL;
	extern int	SNAP, SNAP_FORMAT, L;
	extern float	TIME, DT, *FL, TAU;
	extern float	DAMPLE, DAMPRI, DAMPBA, DAMPFR, DAMPTO, DAMPBO;
	extern float	XREC1[3], XREC2[3];
	extern int	SEISMO, NDT, DXREC[3], SEIS_FORMAT;
	extern int	READMOD, READREC, IDX[3];
	extern int	ABTYPE_LE, ABTYPE_RI, ABTYPE_BA, ABTYPE_FR, ABTYPE_TO, ABTYPE_BO;
	extern int	FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern int	PERIODIC[3];
	extern float	TSNAP1, TSNAP2, TSNAPINC, REFMOD[3], REFSRC[3], REFREC[3];
	extern char	DX_FILE[STRING_SIZE], DY_FILE[STRING_SIZE], DZ_FILE[STRING_SIZE];
	extern char	SNAP_FILE[STRING_SIZE], SRC_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char	SEIS_FILE[STRING_SIZE];
	extern char	SIGNAL_FILE[STRING_SIZE];
	extern char	MFILE[STRING_SIZE];
	extern int	NP, NPROCX[3], MYID; 

	/* local variables */
	int	l;
	char	s1[STRING_SIZE], s2[STRING_SIZE];


	fprintf(fp,"\n **Message from write_par (printed by PE %d):\n",MYID);
	fprintf(fp,"\n");
	fprintf(fp,"------------------------- Processors ------------------------\n");
	fprintf(fp," Number of PEs in x-direction (NPROCX): %d\n",NPROCX[0]);
	fprintf(fp," Number of PEs in y-direction (NPROCY): %d\n",NPROCX[1]);
	fprintf(fp," Number of PEs in vertical (z-) direction (NPROCZ): %d\n",NPROCX[2]);
	fprintf(fp," Total number of PEs in use: %d\n",NP);
	fprintf(fp,"\n");
	fprintf(fp," ----------------------- Discretization  ---------------------\n");
	fprintf(fp," Total number of gridpoints in x-direction (NX): %i\n", NXG[0]);
	fprintf(fp," Total number of gridpoints in y-direction (NY): %i\n", NXG[1]);
	fprintf(fp," Total number of gridpoints in z-direction (NZ): %i\n", NXG[2]);
	fprintf(fp," File with grid-spacing in x-direction (DX_FILE): %s\n", DX_FILE);
	fprintf(fp," File with grid-spacing in y-direction (DY_FILE): %s\n", DY_FILE);
	fprintf(fp," File with grid-spacing in z-direction (DZ_FILE): %s\n", DZ_FILE);
	fprintf(fp," reference point for coordinate system:\n");
	fprintf(fp," x=%f \ty=%f \tz=%f\n",REFMOD[0], REFMOD[1], REFMOD[2]);
	fprintf(fp," Time of wave propagation (T): %e seconds\n",TIME);
	fprintf(fp," Timestep (DT): %e seconds\n", DT);
	fprintf(fp," Number of timesteps: %i \n",NT);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- SOURCE -----------------------------\n");
	fprintf(fp," reference point for source coordinate system:\n");
	fprintf(fp," x=%f \ty=%f \tz=%f\n",REFSRC[0], REFSRC[1], REFSRC[2]);
	fprintf(fp," reading source positions, time delay, centre frequency \n");
	fprintf(fp," and initial amplitude from ASCII-file \n");
	fprintf(fp,"\t%s\n\n",SRC_FILE);

	fprintf(fp," wavelet of source:");

	switch (SRCSIGNAL){
	case 1 :
		fprintf(fp," Ricker\n");
		break;
	case 2 :
		fprintf(fp," Fuchs-Mueller\n");
		break;
	case 3 :
		fprintf(fp," sinus raised to the power of 3.0 \n");
		
		break;
	case 4 :
		fprintf(fp," reading from ASCII-file \n\t %s\n",SIGNAL_FILE);
		break;
	case 5 :
		fprintf(fp," reading from BIN-file \n\t %s\n",SIGNAL_FILE);
		break;
	default :
		error(" Incorrect specification of source wavelet ! ");
	}

	if (SEISMO){
		fprintf(fp,"\n");
		fprintf(fp," ------------------------- RECEIVER  ------- -------------------\n");
		fprintf(fp," reference point for receiver coordinate system:\n");
		fprintf(fp," x=%f \ty=%f \tz=%f\n",REFREC[0], REFREC[1], REFREC[2]);
		fprintf(fp,"\n");
		if (READREC){
			fprintf(fp," reading receiver positions from file \n");
			fprintf(fp,"\t%s\n\n",REC_FILE);
		} 
		else{
			fprintf(fp," position of first receiver (XREC1,YREC1,ZREC1) = (%e, %e, %e) m\n",XREC1[0],XREC1[1],XREC1[2]);
			fprintf(fp," position of last receiver  (XREC2,YREC2,ZREC2) = (%e, %e, %e) m\n",XREC2[0],XREC2[1],XREC2[2]);
			fprintf(fp," horizontal distance between gridpoints (DXREC) = %i gridpoints\n",DXREC[0]);
			fprintf(fp," horizontal distance between gridpoints (DYREC) = %i gridpoints\n",DXREC[1]);
			fprintf(fp," vertical distance between gridpoints (DZREC)   = %i gridpoints\n",DXREC[2]);
		}
	}

	fprintf(fp,"\n");
	fprintf(fp," ------------------------- MODEL BOUNDARIES --------------------\n");
	if (PERIODIC[0]){
		sprintf(s1,"Periodic boundary condition");
		sprintf(s2,"periodic boundary");
	}
	else{
		sprintf(s1,"Free surface boundary condition");
		sprintf(s2,"free surface boundary");
	}
	switch (ABTYPE_LE){
	case 0 :
		fprintf(fp," %s will be applied at the left boundary.\n\n",s1);
		break;
	case 1 :
		fprintf(fp," PML boundary condition \n");
		fprintf(fp," together with %s will be applied at the left boundary.\n",s2);
		fprintf(fp," Width of PML zone is %d grid points.\n\n",FWLE);
		break;
	case 2 :
		fprintf(fp," Exponential damping of the wavefield within the absorbing frame \n");
		fprintf(fp," together with %s will be applied at the left boundary. \n",s2);
		fprintf(fp," Width of absorbing frame is %d gridpoints.\n",FWLE);
		fprintf(fp," Percentage of amplitude decay: %f .\n\n",DAMPLE);
		break;
	default :
		warning(" Wrong integer value for ABTYPE_LE specified ! \n");
		warning(" No absorbing boundary condition will be applied at the left boundary. \n");
		ABTYPE_LE=0;
		break;
	}
	switch (ABTYPE_RI){
	case 0 :
		fprintf(fp," %s will be applied at the right boundary.\n\n",s1);
		break;
	case 1 :
		fprintf(fp," PML boundary condition \n");
		fprintf(fp," together with %s will be applied at the right boundary.\n",s2);
		fprintf(fp," Width of PML zone is %d grid points.\n\n",FWRI);
		break;
	case 2 :
		fprintf(fp," Exponential damping of the wavefield within the absorbing frame \n");
		fprintf(fp," together with %s will be applied at the right boundary. \n",s2);
		fprintf(fp," Width of absorbing frame is %d gridpoints.\n",FWRI);
		fprintf(fp," Percentage of amplitude decay: %f .\n\n",DAMPRI);
		break;
	default :
		warning(" Wrong integer value for ABTYPE_RI specified ! \n");
		warning(" No absorbing boundary condition will be applied at the right boundary. \n");
		ABTYPE_RI=0;
		break;
	}

	if (PERIODIC[1]){
		sprintf(s1,"Periodic boundary condition");
		sprintf(s2,"periodic boundary");
	}
	else{
		sprintf(s1,"Free surface boundary condition");
		sprintf(s2,"free surface boundary");
	}
	switch (ABTYPE_BA){
	case 0 :
		fprintf(fp," %s will be applied at the back boundary.\n\n",s1);
		break;
	case 1 :
		fprintf(fp," PML boundary condition \n");
		fprintf(fp," together with %s will be applied at the back boundary.\n",s2);
		fprintf(fp," Width of PML zone is %d grid points.\n\n",FWBA);
		break;
	case 2 :
		fprintf(fp," Exponential damping of the wavefield within the absorbing frame \n");
		fprintf(fp," together with %s will be applied at the back boundary. \n",s2);
		fprintf(fp," Width of absorbing frame is %d gridpoints.\n",FWBA);
		fprintf(fp," Percentage of amplitude decay: %f .\n\n",DAMPBA);
		break;
	default :
		warning(" Wrong integer value for ABTYPE_BA specified ! \n");
		warning(" No absorbing boundary condition will be applied at the back boundary. \n");
		ABTYPE_BA=0;
		break;
	}
	switch (ABTYPE_FR){
	case 0 :
		fprintf(fp," %s will be applied at the front boundary.\n\n",s1);
		break;
	case 1 :
		fprintf(fp," PML boundary condition \n");
		fprintf(fp," together with %s will be applied at the front boundary.\n",s2);
		fprintf(fp," Width of PML zone is %d grid points.\n\n",FWFR);
		break;
	case 2 :
		fprintf(fp," Exponential damping of the wavefield within the absorbing frame \n");
		fprintf(fp," together with %s will be applied at the front boundary. \n",s2);
		fprintf(fp," Width of absorbing frame is %d gridpoints.\n",FWFR);
		fprintf(fp," Percentage of amplitude decay: %f .\n\n",DAMPFR);
		break;
	default :
		warning(" Wrong integer value for ABTYPE_FR specified ! \n");
		warning(" No absorbing boundary condition will be applied at the front boundary. \n");
		ABTYPE_FR=0;
		break;
	}

	if (PERIODIC[2]){
		sprintf(s1,"Periodic boundary condition");
		sprintf(s2,"periodic boundary");
	}
	else{
		sprintf(s1,"Free surface boundary condition");
		sprintf(s2,"free surface boundary");
	}
	switch (ABTYPE_TO){
	case 0 :
		fprintf(fp," %s will be applied at the top boundary.\n\n",s1);
		break;
	case 1 :
		fprintf(fp," PML boundary condition \n");
		fprintf(fp," together with %s will be applied at the top boundary.\n",s2);
		fprintf(fp," Width of PML zone is %d grid points.\n\n",FWTO);
		break;
	case 2 :
		fprintf(fp," Exponential damping of the wavefield within the absorbing frame \n");
		fprintf(fp," together with %s will be applied at the top boundary. \n",s2);
		fprintf(fp," Width of absorbing frame is %d gridpoints.\n",FWTO);
		fprintf(fp," Percentage of amplitude decay: %f .\n\n",DAMPTO);
		break;
	default :
		warning(" Wrong integer value for ABTYPE_TO specified ! \n");
		warning(" No absorbing boundary condition will be applied at the top boundary. \n");
		ABTYPE_TO=0;
		break;
	}
	switch (ABTYPE_BO){
	case 0 :
		fprintf(fp," %s will be applied at the bottom boundary.\n\n",s1);
		break;
	case 1 :
		fprintf(fp," PML boundary condition \n");
		fprintf(fp," together with %s will be applied at the bottom boundary.\n",s2);
		fprintf(fp," Width of PML zone is %d grid points.\n\n",FWBO);
		break;
	case 2 :
		fprintf(fp," Exponential damping of the wavefield within the absorbing frame \n");
		fprintf(fp," together with %s will be applied at the bottom boundary. \n",s2);
		fprintf(fp," Width of absorbing frame is %d gridpoints.\n",FWBO);
		fprintf(fp," Percentage of amplitude decay: %f .\n\n",DAMPBO);
		break;
	default :
		warning(" Wrong integer value for ABTYPE_BO specified ! \n");
		warning(" No absorbing boundary condition will be applied at the bottom boundary. \n");
		ABTYPE_BO=0;
		break;
	}

	fprintf(fp,"\n");
	fprintf(fp," ------------------------- MODEL -------------------------------\n");
	if (READMOD){
		fprintf(fp," names of model-files: \n");
		fprintf(fp,"\t compressional wave velocities:\n\t %s.vp\n",MFILE);
		fprintf(fp,"\t density:\n\t %s.rho\n",MFILE);
		if (L){
			fprintf(fp,"\t Q for P-waves:\n\t %s.qp\n",MFILE);
		}
	}
	else{
		fprintf(fp," use model defined in source code (model_p.c) \n");
	}

	if (L){
		fprintf(fp,"\n");
		fprintf(fp," ------------------------- Q-APROXIMATION --------------------\n");
		fprintf(fp," Number of relaxation mechanisms (L): %i\n",L);
		fprintf(fp," The L relaxation frequencies are at:  \n");
		for (l=1;l<=L;l++)
			fprintf(fp,"\t%f",FL[l]);
		fprintf(fp," Hz\n");
		fprintf(fp," Global value for tau is : %f\n (is ignored if Qp is read from model file !)",TAU);
	}

	if (SNAP & 15){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
		fprintf(fp," Snapshots of\n");

		if (SNAP & 8)
			fprintf(fp,"\t x-, y- and z-component of particle acceleration\n");
		if (SNAP & 4)
			fprintf(fp,"\t divergence of the wavefield\n");
		if (SNAP & 2)
			fprintf(fp,"\t pressure field\n");
		if (SNAP & 1)
			fprintf(fp,"\t x-, y- and z-component of particle velocity\n");

		fprintf(fp,"\n");
		fprintf(fp," first (TSNAP1)            =%8.5f s\n",TSNAP1);
		fprintf(fp," last (TSNAP2)             =%8.5f s\n",TSNAP2);
		fprintf(fp," time increment (TSNAPINC) =%8.5f s\n\n",TSNAPINC);
		fprintf(fp," spatial increment in x-direction =%d grid points\n",IDX[0]);
		fprintf(fp," spatial increment in y-direction =%d grid points\n",IDX[1]);
		fprintf(fp," spatial increment in z-direction =%d grid points\n",IDX[2]);
		fprintf(fp,"\n basic name of output-file (SNAP_FILE):\n\t %s\n",SNAP_FILE);

		switch (SNAP_FORMAT){
		case 1 :
			error(" SU-Format not yet available !!");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary (native byte order) (4 bytes per sample)");
			break;
		default:
			error(" Format for the snapshot-data not known ! \n");
		}
	}

	if (SEISMO & 15){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SEISMOGRAMS  ----------------------\n");
		fprintf(fp," Seismograms of\n");

		if (SEISMO & 8)
			fprintf(fp,"\t x-, y- and z-component of particle acceleration\n");
		if (SEISMO & 4)
			fprintf(fp,"\t divergence of the wavefield\n");
		if (SEISMO & 2)
			fprintf(fp,"\t pressure field\n");
		if (SEISMO & 1)
			fprintf(fp,"\t x-, y- and z-component of particle velocity\n");

		fprintf(fp," \n basic name of output-file (SEIS_FILE):\n\t %s\n",SEIS_FILE);

		switch (SEIS_FORMAT){
		case 1 :
			fprintf(fp," The seismograms are saved in SU-format (native byte order) . \n");
			break;
		case 2 :
			fprintf(fp," The seismograms are written in ASCII format. \n");
			break;
		case 3 :
			fprintf(fp," The data is written in native binary format (4 byte per sample)");
			break;
		default:
			error(" Format for the seismic data not known ! \n");
		}
		fprintf(fp," samplingrate of seismic data: %f s\n",NDT*DT);
		fprintf(fp," Number of samples per trace: %i \n", iround(NT/NDT));
		fprintf(fp," ----------------------------------------------------------\n");
	}
	fprintf(fp,"\n");
	fprintf(fp," **************************************************************\n");
	fprintf(fp,"\n");

}
