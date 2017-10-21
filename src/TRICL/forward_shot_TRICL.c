/*  ---------------------------------------------------------------------------
 *   Anisotropic elastic forward modelling (triclinic medium) for one shot
 *
 *
 *   D. Koehn
 *   Kiel, 20.10.2017
 *
 *  ---------------------------------------------------------------------------*/

#include "fd.h"		    /* general include file for all FD programs */
#include "fd3d.h"	    /* general include file for FD3D program */
#include "fd3daniso.h"	    /* include file for anisotropic version */
#include "fd3dtricl.h"	    /* include file for triclinic version */
#include "fd3del.h"	    /* include file for elastic version */
#include "TRICL_struct.h"   /* data structures for TRICL forward modelling */

void forward_shot_TRICL(struct wave *wave, struct pmls *pmls, struct mat *mat, struct geom *geom, struct mpi *mpi, 
                     struct seis *seis, struct acq *acq, struct times *times, int ns){

        /* global variables */
	extern int NT, NX[3], SEISMO, SNAP, NDT, L, MYID, PML;
	extern int NTR, NTR_LOC, NSRC, NSRC_LOC;
	extern float TSNAP1, TSNAP2, TSNAPINC, DT;
	extern FILE *FP;

        /* local variables */
	int nt, infoout, lsamp, lsnap, nsnap=0;	

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
		(*times).time_update[1] = update_v(1, NX[0], 1, NX[1], 1, NX[2], (*wave).v, (*wave).t, (*wave).a, 
					(*mat).rhoijpkp, (*mat).rhoipjkp, (*mat).rhoipjpk, (*pmls).absorb_coeff,
					(*geom).dx, (*geom).dxp, (*geom).dy, (*geom).dyp, (*geom).dz, (*geom).dzp, infoout); 

		if (PML) (*times).time_update[2] = pml_update_v((*wave).v, (*wave).t, (*mat).rhoijpkp, (*mat).rhoipjkp, (*mat).rhoipjpk,
					(*pmls).pmlle, (*pmls).pmlri, (*pmls).pmlba, (*pmls).pmlfr, (*pmls).pmlto, (*pmls).pmlbo, 
					(*pmls).Pt_x_l, (*pmls).Pt_x_r, (*pmls).Pt_y_b, (*pmls).Pt_y_f, (*pmls).Pt_z_t, (*pmls).Pt_z_b,
					(*geom).dx, (*geom).dxp, (*geom).dy, (*geom).dyp, (*geom).dz, (*geom).dzp, infoout);

		/* apply force source */
		fsource(nt,(*wave).v,(*mat).rhoijpkp,(*mat).rhoipjkp,(*mat).rhoipjpk,(*acq).srcpos_loc,(*acq).signals,NSRC_LOC,(*geom).dx,(*geom).dy,(*geom).dz,(*geom).dxp,(*geom).dyp,(*geom).dzp);

		/* exchange of particle velocities between PEs */
		(*times).time_update[3] = exchange_v((*wave).v, 
					(*mpi).sbuf_v_lef_to_rig, (*mpi).sbuf_v_rig_to_lef, (*mpi).sbuf_v_bac_to_fro, 
					(*mpi).sbuf_v_fro_to_bac, (*mpi).sbuf_v_top_to_bot, (*mpi).sbuf_v_bot_to_top, 
					(*mpi).rbuf_v_lef_to_rig, (*mpi).rbuf_v_rig_to_lef, (*mpi).rbuf_v_bac_to_fro, 
					(*mpi).rbuf_v_fro_to_bac, (*mpi).rbuf_v_top_to_bot, (*mpi).rbuf_v_bot_to_top, 
					(*mpi).request_v, 
					infoout);

		/* update of stress tensor components */
		if (L)  /* viscoelastic simulation */
			/*(*times).time_update[4] = update_s_ve(1, NX[0], 1, NX[1], 1, NX[2], (*wave).v, (*wave).e, (*wave).r, (*wave).w, 
					pi, u, uipjk, uijpk, uijkp,
					(*mat).taup, (*mat).taus, (*mat).tausipjk, (*mat).tausijpk, (*mat).tausijkp, (*mat).eta, 
					(*pmls).absorb_coeff,
					(*geom).dx, (*geom).dxp, (*geom).dy, (*geom).dyp, (*geom).dz, (*geom).dzp, 
					infoout);*/
			fprintf(FP,"not yet implemented\n\n");
		else  /* elastic simulation (no visco-elasticity) */
			(*times).time_update[4] = update_e_el(1, NX[0], 1, NX[1], 1, NX[2], (*wave).v, (*wave).e, (*wave).w,
							   (*geom).dx, (*geom).dxp, (*geom).dy, (*geom).dyp, (*geom).dz, (*geom).dzp, infoout);

		if (PML) (*times).time_update[5] = pml_update_e((*wave).v, (*wave).e,
						(*pmls).pmlle, (*pmls).pmlri, (*pmls).pmlba, (*pmls).pmlfr, (*pmls).pmlto, (*pmls).pmlbo, 
						(*pmls).Pv_x_l, (*pmls).Pv_x_r, (*pmls).Pv_y_b, (*pmls).Pv_y_f, (*pmls).Pv_z_t, (*pmls).Pv_z_b,
						(*geom).dx, (*geom).dxp, (*geom).dy, (*geom).dyp, (*geom).dz, (*geom).dzp, infoout);

		/* exchange of strain components between PEs */
		(*times).time_update[6] = exchange_e((*wave).e, 
				       (*mpi).sbuf_e_lef_to_rig, (*mpi).sbuf_e_rig_to_lef, (*mpi).sbuf_e_bac_to_fro, 
				       (*mpi).sbuf_e_fro_to_bac, (*mpi).sbuf_e_top_to_bot, (*mpi).sbuf_e_bot_to_top, 
				       (*mpi).rbuf_e_lef_to_rig, (*mpi).rbuf_e_rig_to_lef, (*mpi).rbuf_e_bac_to_fro, 
				       (*mpi).rbuf_e_fro_to_bac, (*mpi).rbuf_e_top_to_bot, (*mpi).rbuf_e_bot_to_top, 
				       (*mpi).request_e, 
				       infoout);

		/* update of stress components */
		(*times).time_update[7] = update_s_el(1, NX[0], 1, NX[1], 1, NX[2], (*wave).e, (*wave).t, 
				       (*mat).c1111,  (*mat).c1122,  (*mat).c1133,  (*mat).c1123,  (*mat).c1113,   (*mat).c1112, 
				       (*mat).c2222,  (*mat).c2233,  (*mat).c2223,  (*mat).c1322,   (*mat).c1222, 
				       (*mat).c3333,  (*mat).c2333,  (*mat).c1333,   (*mat).c1233, 
				       (*mat).c1123h, (*mat).c2223h, (*mat).c2333h, (*mat).c2323h, (*mat).c1323ha, (*mat).c1223ha, 
				       (*mat).c1113h, (*mat).c1322h, (*mat).c1333h, (*mat).c1323h, (*mat).c1313h,  (*mat).c1213ha, 
				       (*mat).c1112h, (*mat).c1222h, (*mat).c1233h, (*mat).c1223h, (*mat).c1213h,  (*mat).c1212h, 
				       (*pmls).absorb_coeff,
				       (*geom).wxp1, (*geom).wx1, (*geom).wxp2, (*geom).wx2, (*geom).wyp1, (*geom).wy1, (*geom).wyp2, (*geom).wy2, (*geom).wzp1, (*geom).wz1, (*geom).wzp2, (*geom).wz2,
				       infoout);

		/* apply explosive source */
		psource(nt,(*wave).t,(*acq).srcpos_loc,(*acq).signals,NSRC_LOC,(*geom).dx,(*geom).dy,(*geom).dz,(*geom).dxp,(*geom).dyp,(*geom).dzp);

		/* stress exchange between PEs */
		(*times).time_update[8] = exchange_s((*wave).t,
					(*mpi).sbuf_s_lef_to_rig, (*mpi).sbuf_s_rig_to_lef, (*mpi).sbuf_s_bac_to_fro, 
					(*mpi).sbuf_s_fro_to_bac, (*mpi).sbuf_s_top_to_bot, (*mpi).sbuf_s_bot_to_top, 
					(*mpi).rbuf_s_lef_to_rig, (*mpi).rbuf_s_rig_to_lef, (*mpi).rbuf_s_bac_to_fro, 
					(*mpi).rbuf_s_fro_to_bac, (*mpi).rbuf_s_top_to_bot, (*mpi).rbuf_s_bot_to_top, 
					(*mpi).request_s, 
					infoout);

		/* store amplitudes at receivers in section-arrays */
		if ((SEISMO) && (nt==lsamp) && (nt<NT) && (NTR_LOC>0)){
			seismo(lsamp, NTR_LOC, (*acq).recpos_loc, 
				(*seis).sectionvx, (*seis).sectionvy, (*seis).sectionvz, (*seis).sectionax, (*seis).sectionay, (*seis).sectionaz, 
				(*seis).sectiondiv, (*seis).sectioncurlx, (*seis).sectioncurly, (*seis).sectioncurlz, (*seis).sectionp,
				(*seis).sectiontxx, (*seis).sectiontxy, (*seis).sectiontxz, (*seis).sectiontyy, (*seis).sectiontyz, (*seis).sectiontzz, 
				(*wave).v, (*wave).t, (*wave).a, (*wave).w);
			lsamp += NDT;
		}

		/* write snapshots to disk */
		if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
			snap(FP, nt, ++nsnap, (*wave).v, (*wave).t, (*wave).a, (*wave).w, (*geom).xp, (*geom).yp, (*geom).zp);
			lsnap += iround(TSNAPINC/DT);
		}

		/* timing */
		(*times).time3          = MPI_Wtime();
		(*times).time_update[9] = ((*times).time3-(*times).time2);
		if (infoout)
			fprintf(FP," Total real time for timestep %d : \t\t %4.2f s.\n",nt,(*times).time_update[9]);

		/* update timing statistics */
		update_timing(nt, (*times).time_update, (*times).time_sum, (*times).time_avg, (*times).time_std, 9);

	}
	/*-------------------- End of loop over timesteps ----------*/
	MPI_Barrier(MPI_COMM_WORLD);

	/* save seismograms */
	if (SEISMO){
       		saveseis(FP,(*seis).sectionvx,(*seis).sectionvy,(*seis).sectionvz,(*seis).sectionax,(*seis).sectionay,(*seis).sectionaz,
			 (*seis).sectiondiv,(*seis).sectioncurlx,(*seis).sectioncurly,(*seis).sectioncurlz,(*seis).sectionp,
			 (*seis).sectiontxx,(*seis).sectiontxy,(*seis).sectiontxz,(*seis).sectiontyy,(*seis).sectiontyz,(*seis).sectiontzz,
                	 (*seis).section_fulldata,(*acq).recpos,(*acq).recpos_loc,(*acq).srcpos,NSRC,ns,(*geom).xg,(*geom).yg,(*geom).zg,
                	 (*geom).xpg,(*geom).ypg,(*geom).zpg,(*acq).recswitch);
	}

	/* merge snapshots */
	if ((!(MYID))&&(SNAP)){
		savesnap(FP, (*geom).xpg, (*geom).ypg, (*geom).zpg);
	}

	MPI_Barrier(MPI_COMM_WORLD);

}
