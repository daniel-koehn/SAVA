/*  -------------------------------------------
 *   Acoustic forward modelling for one shot
 *
 *
 *   D. Koehn
 *   Kiel, 19.10.2017
 *
 *  -------------------------------------------*/

#include "fd.h"
#include "fd3d.h"	/* general include file for program FD3D */
#include "fd3dp.h"	/* include file for acoustic version */
#include "fd3diso.h"	/* include file for isotropic version */
#include "AC_struct.h"	/* data structures for AC forward modelling */

void forward_shot_AC(struct wave *wave, struct pmls *pmls, struct mat *mat, struct geom *geom, struct mpi *mpi, 
                     struct seis *seis, struct acq *acq, struct times *times, int ns){

        /* global variables */
	extern int NT, NX[3], SEISMO, SNAP, NDT, L, MYID, PML;
	extern int NTR, NTR_LOC, NSRC, NSRC_LOC;
	extern float TSNAP1, TSNAP2, TSNAPINC, DT;
	extern FILE *FP;

        /* local variables */
	int nt, infoout, lsamp, lsnap, nsnap=0;	

	(*times).time2 = MPI_Wtime();

	fprintf(FP,"\n\n\n *********** STARTING TIME STEPPING ***************\n");
	fprintf(FP," real time before starting time loop: %4.2f s.\n",(*times).time2-(*times).time1);

	lsamp = NDT;                    /* seismogram sampling rate */
	lsnap = iround(TSNAP1/DT);	/* first snapshot at this timestep */

	/*----------------------  loop over timesteps  ------------------*/

	for (nt=1;nt<=NT;nt++){

		infoout = !(nt%50);

		/* timing */
		(*times).time2 = MPI_Wtime();

		if (infoout){
			fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
		}
	
		/* Check if simulation is still stable */
		if (isnan((*wave).v[NX[0]/2][NX[1]/2][NX[2]/2].x)) error(" Simulation is unstable !");
	
		/* update of particle velocities */
		(*times).time_update[1] = update_v(1, NX[0], 1, NX[1], 1, NX[2], (*wave).v, (*wave).p, (*wave).a, (*mat).rhoijpkp, (*mat).rhoipjkp, (*mat).rhoipjpk, (*pmls).absorb_coeff,
						(*geom).dxp, (*geom).dyp, (*geom).dzp, infoout); 

		if (PML) (*times).time_update[2] = pml_update_v((*wave).v, (*wave).p, (*mat).rhoijpkp, (*mat).rhoipjkp, (*mat).rhoipjpk,
						(*pmls).pmlle, (*pmls).pmlri, (*pmls).pmlba, (*pmls).pmlfr, (*pmls).pmlto, (*pmls).pmlbo, 
						(*pmls).Pp_x_l, (*pmls).Pp_x_r, (*pmls).Pp_y_b, (*pmls).Pp_y_f, (*pmls).Pp_z_t, (*pmls).Pp_z_b,
						(*geom).dxp, (*geom).dyp, (*geom).dzp, infoout);

		/* apply body force source */
		fsource(nt,(*wave).v,(*mat).rhoijpkp,(*mat).rhoipjkp,(*mat).rhoipjpk,(*acq).srcpos_loc,(*acq).signals,NSRC_LOC,(*geom).dx,(*geom).dy,(*geom).dz,(*geom).dxp,(*geom).dyp,(*geom).dzp);

		/* exchange of particle velocities between PEs */
		(*times).time_update[3] = exchange_v((*wave).v, 
						  (*mpi).sbuf_lef_to_rig, (*mpi).sbuf_bac_to_fro, (*mpi).sbuf_top_to_bot, 
						  (*mpi).rbuf_lef_to_rig, (*mpi).rbuf_bac_to_fro, (*mpi).rbuf_top_to_bot, 
						  (*mpi).request_v,
						  infoout);

		/* update of stress tensor components */
		if (L)  /* visco-acoustic simulation */
			(*times).time_update[4] = update_s_va(1, NX[0], 1, NX[1], 1, NX[2], (*wave).v, (*wave).p, (*wave).rpp, (*wave).diverg,
						           (*mat).pi, (*mat).taup, (*mat).eta, 
							   (*pmls).absorb_coeff, 
							   (*geom).dx, (*geom).dy, (*geom).dz, infoout);
		else  /* acoustic simulation (no visco-elasticity) */
			(*times).time_update[4] = update_s_ac(1, NX[0], 1, NX[1], 1, NX[2], (*wave).v, (*wave).p, (*wave).diverg,
							   (*mat).pi, 
							   (*pmls).absorb_coeff,
							   (*geom).dx, (*geom).dy, (*geom).dz, infoout);

		if (PML) (*times).time_update[5] = pml_update_s((*wave).v, (*wave).p, (*mat).pi, 
						             (*pmls).pmlle, (*pmls).pmlri, (*pmls).pmlba, (*pmls).pmlfr, (*pmls).pmlto, (*pmls).pmlbo, 
							     (*pmls).Pvx_x_l, (*pmls).Pvx_x_r, (*pmls).Pvy_y_b, (*pmls).Pvy_y_f, (*pmls).Pvz_z_t, (*pmls).Pvz_z_b,
							     (*geom).dx, (*geom).dy, (*geom).dz, infoout);

		/* apply explosive source */
		psource(nt,(*wave).p,(*acq).srcpos_loc,(*acq).signals,NSRC_LOC,(*geom).dx,(*geom).dy,(*geom).dz);

		/* stress exchange between PEs */
		(*times).time_update[6] = exchange_s((*wave).p, 
						  (*mpi).sbuf_rig_to_lef, (*mpi).sbuf_fro_to_bac, (*mpi).sbuf_bot_to_top, 
						  (*mpi).rbuf_rig_to_lef, (*mpi).rbuf_fro_to_bac, (*mpi).rbuf_bot_to_top, 
						  (*mpi).request_s,
						  infoout);


		/* store amplitudes at receivers in section-arrays */
		if ((SEISMO) && (nt==lsamp) && (nt<NT) && (NTR_LOC>0)){
			seismo(lsamp, (*acq).recpos_loc, 
			       (*seis).sectionvx, (*seis).sectionvy, (*seis).sectionvz, (*seis).sectionp, 
			       (*seis).sectionax, (*seis).sectionay, (*seis).sectionaz, (*seis).sectiondiv, 
			       (*wave).v, (*wave).p, (*wave).a, (*wave).diverg);
			lsamp += NDT;
		}

		/* write snapshots to disk */
		if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
			snap(FP, nt, ++nsnap, (*wave).v, (*wave).p, (*wave).a, (*wave).diverg, (*geom).xp, (*geom).yp, (*geom).zp);
			lsnap += iround(TSNAPINC/DT);
		}

		/* timing */
		(*times).time3          = MPI_Wtime();
		(*times).time_update[7] = ((*times).time3-(*times).time2);
		if (infoout)
			fprintf(FP," Total real time for timestep %d : \t\t %4.2f s.\n",nt,(*times).time_update[7]);

		/* update timing statistics */
		update_timing(nt, (*times).time_update, (*times).time_sum, (*times).time_avg, (*times).time_std, 7);

	}
	/*-------------------- End of loop over timesteps ----------*/
	MPI_Barrier(MPI_COMM_WORLD);

	/* save seismograms */
	if (SEISMO){
		saveseis(FP,(*seis).sectionvx,(*seis).sectionvy,(*seis).sectionvz,(*seis).sectionp,(*seis).sectionax,(*seis).sectionay,
                 	 (*seis).sectionaz,(*seis).sectiondiv,(*seis).section_fulldata,(*acq).recpos,(*acq).recpos_loc,(*acq).srcpos,
                 	 NSRC,ns,(*geom).xg,(*geom).yg,(*geom).zg,(*geom).xpg,(*geom).ypg,(*geom).zpg,(*acq).recswitch);
	}

	/* merge snapshots */
	if ((!(MYID))&&(SNAP)){
		savesnap(FP, (*geom).xpg, (*geom).ypg, (*geom).zpg);
	}

	MPI_Barrier(MPI_COMM_WORLD);
       
}
