/*  -------------------------------------------
 *   Reset wavefields for elastic modelling
 *
 *
 *   D. Koehn
 *   Kiel, 22.10.2017
 *
 *  -------------------------------------------*/

#include "fd.h"
#include "TRICL_struct.h"	/* data structures for anisotropic elastic forward modelling (triclinic medium) */

void reset_wavefields_tricl(struct wave *wave, struct pmls *pmls, struct seis *seis, int ns){

        /* global variables */
	extern int NX[3], L;
	extern int PML_LE, FWLE;
	extern int PML_RI, FWRI;
	extern int PML_BA, FWBA;
	extern int PML_FR, FWFR;
	extern int PML_TO, FWTO;
	extern int PML_BO, FWBO;

        /* local variables */
	int i, j, k;	

	/* reset elastic wave fields */	
	for (i=0;i<=NX[0]+1;i++){
		for (j=0;j<=NX[1]+1;j++){
		        for (k=0;k<=NX[2]+1;k++){

				(*wave).v[i][j][k].x = 0.0;
				(*wave).v[i][j][k].y = 0.0;
				(*wave).v[i][j][k].z = 0.0;

				(*wave).t[i][j][k].xx = 0.0;
				(*wave).t[i][j][k].yy = 0.0;
				(*wave).t[i][j][k].zz = 0.0;
				(*wave).t[i][j][k].yz = 0.0;
				(*wave).t[i][j][k].xz = 0.0;
				(*wave).t[i][j][k].xy = 0.0;

			}
		}
	}

	/* reset visco-elastic wave fields */
	if(L){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
		        	for (k=1;k<=NX[2];k++){

					(*wave).r[i][j][k].xx = 0.0;
					(*wave).r[i][j][k].yy = 0.0;
					(*wave).r[i][j][k].zz = 0.0;
					(*wave).r[i][j][k].yz = 0.0;
					(*wave).r[i][j][k].xz = 0.0;
					(*wave).r[i][j][k].xy = 0.0;

				}
			}
		}
	
	}

	/* reset PML variables */
	if(PML_LE){
		for (i=1;i<=FWLE;i++){
			for (j=1;j<=NX[1];j++){
		        	for (k=1;k<=NX[2];k++){

					(*pmls).Pt_x_l[i][j][k].x = 0.0;
					(*pmls).Pt_x_l[i][j][k].y = 0.0;
					(*pmls).Pt_x_l[i][j][k].z = 0.0;

					(*pmls).Pv_x_l[i][j][k].x = 0.0;
					(*pmls).Pv_x_l[i][j][k].y = 0.0;
					(*pmls).Pv_x_l[i][j][k].z = 0.0;

				}
			}
		}			
	}

	if(PML_RI){
		for (i=NX[0]-FWRI+1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
		        	for (k=1;k<=NX[2];k++){

					(*pmls).Pt_x_r[i][j][k].x = 0.0;
					(*pmls).Pt_x_r[i][j][k].y = 0.0;
					(*pmls).Pt_x_r[i][j][k].z = 0.0;

					(*pmls).Pv_x_r[i][j][k].x = 0.0;
					(*pmls).Pv_x_r[i][j][k].y = 0.0;
					(*pmls).Pv_x_r[i][j][k].z = 0.0;

				}
			}
		}			
	}
	
	if(PML_BA){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=FWBA;j++){
		        	for (k=1;k<=NX[2];k++){

					(*pmls).Pt_y_b[i][j][k].x = 0.0;
					(*pmls).Pt_y_b[i][j][k].y = 0.0;
					(*pmls).Pt_y_b[i][j][k].z = 0.0;

					(*pmls).Pv_y_b[i][j][k].x = 0.0;
					(*pmls).Pv_y_b[i][j][k].y = 0.0;
					(*pmls).Pv_y_b[i][j][k].z = 0.0;

				}
			}
		}			
	}

	if(PML_FR){
		for (i=1;i<=NX[0];i++){
			for (j=NX[1]-FWFR+1;j<=NX[1];j++){
		        	for (k=1;k<=NX[2];k++){

					(*pmls).Pt_y_f[i][j][k].x = 0.0;
					(*pmls).Pt_y_f[i][j][k].y = 0.0;
					(*pmls).Pt_y_f[i][j][k].z = 0.0;

					(*pmls).Pv_y_f[i][j][k].x = 0.0;
					(*pmls).Pv_y_f[i][j][k].y = 0.0;
					(*pmls).Pv_y_f[i][j][k].z = 0.0;

				}
			}
		}			
	}

	if(PML_TO){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
		        	for (k=1;k<=FWTO;k++){

					(*pmls).Pt_z_t[i][j][k].x = 0.0;
					(*pmls).Pt_z_t[i][j][k].y = 0.0;
					(*pmls).Pt_z_t[i][j][k].z = 0.0;

					(*pmls).Pv_z_t[i][j][k].x = 0.0;
					(*pmls).Pv_z_t[i][j][k].y = 0.0;
					(*pmls).Pv_z_t[i][j][k].z = 0.0;

				}
			}
		}			
	}


	if(PML_BO){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
		        	for (k=NX[2]-FWBO+1;k<=NX[2];k++){

					(*pmls).Pt_z_b[i][j][k].x = 0.0;
					(*pmls).Pt_z_b[i][j][k].y = 0.0;
					(*pmls).Pt_z_b[i][j][k].z = 0.0;

					(*pmls).Pv_z_b[i][j][k].x = 0.0;
					(*pmls).Pv_z_b[i][j][k].y = 0.0;
					(*pmls).Pv_z_b[i][j][k].z = 0.0;

				}
			}
		}			
	}
       
}
