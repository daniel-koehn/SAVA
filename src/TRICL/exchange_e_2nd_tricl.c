/*------------------------------------------------------------------------
 *   write values of deformation tensor at the edges of the
 *   local grid into buffer arrays and exchange between
 *   processes.
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


double exchange_e(struct tensor3d ***e,
float *** sbuf_lef_to_rig, float *** sbuf_rig_to_lef, 
float *** sbuf_bac_to_fro, float *** sbuf_fro_to_bac, 
float *** sbuf_top_to_bot, float *** sbuf_bot_to_top, 
float *** rbuf_lef_to_rig, float *** rbuf_rig_to_lef, 
float *** rbuf_bac_to_fro, float *** rbuf_fro_to_bac, 
float *** rbuf_top_to_bot, float *** rbuf_bot_to_top, 
MPI_Request * request, 
int infoout){

	/* extern variables */
	extern int	NX[3];
	extern FILE	*FP;

	/* MPI variables */
	MPI_Status	status[12];

	/* local variables */
	int	i, j, k;
	double	time=0, time1=0, time2=0;


	/* timing */
	time1 = MPI_Wtime();

	/* exchange top-bottom and bottom-top */
	for (i=1;i<=NX[0];i++){
		for (j=1;j<=NX[1];j++){
			/* storage of top of local volume into buffer */
			sbuf_top_to_bot[i][j][1] = e[i][j][1].yz;
			sbuf_top_to_bot[i][j][2] = e[i][j][1].xz;

			/* storage of bottom of local volume into buffer */
			sbuf_bot_to_top[i][j][1] = e[i][j][NX[2]].xx;
			sbuf_bot_to_top[i][j][2] = e[i][j][NX[2]].yy;
			sbuf_bot_to_top[i][j][3] = e[i][j][NX[2]].zz;
			sbuf_bot_to_top[i][j][4] = e[i][j][NX[2]].xy;
		}
	}
	/* exchange back-front and front-back */
	for (i=1;i<=NX[0];i++){
		for (k=1;k<=NX[2];k++){
			/* storage of top of local volume into buffer */
			sbuf_bac_to_fro[i][k][1] = e[i][1][k].yz;
			sbuf_bac_to_fro[i][k][2] = e[i][1][k].xy;

			/* storage of bottom of local volume into buffer */
			sbuf_fro_to_bac[i][k][1] = e[i][NX[1]][k].xx;
			sbuf_fro_to_bac[i][k][2] = e[i][NX[1]][k].yy;
			sbuf_fro_to_bac[i][k][3] = e[i][NX[1]][k].zz;
			sbuf_fro_to_bac[i][k][4] = e[i][NX[1]][k].xz;
		}
	}
	/* exchange left-right and right-left */
	for (j=1;j<=NX[1];j++){
		for (k=1;k<=NX[2];k++){
			/* storage of left edge of local volume into buffer */
			sbuf_lef_to_rig[j][k][1] = e[1][j][k].xz;
			sbuf_lef_to_rig[j][k][2] = e[1][j][k].xy;

			/* storage of right edge of local volume into buffer */
			sbuf_rig_to_lef[j][k][1] = e[NX[0]][j][k].xx;
			sbuf_rig_to_lef[j][k][2] = e[NX[0]][j][k].yy;
			sbuf_rig_to_lef[j][k][3] = e[NX[0]][j][k].zz;
			sbuf_rig_to_lef[j][k][4] = e[NX[0]][j][k].yz;
		}
	}

	MPI_Startall(12,&request[0]);
	MPI_Waitall(12, &request[0],&status[0]);

	/* exchange top-bottom and bottom-top */
	for (i=1;i<=NX[0];i++){
		for (j=1;j<=NX[1];j++){
			e[i][j][0].xx = rbuf_bot_to_top[i][j][1];
			e[i][j][0].yy = rbuf_bot_to_top[i][j][2];
			e[i][j][0].zz = rbuf_bot_to_top[i][j][3];
			e[i][j][0].xy = rbuf_bot_to_top[i][j][4];

			e[i][j][NX[2]+1].yz = rbuf_top_to_bot[i][j][1];
			e[i][j][NX[2]+1].xz = rbuf_top_to_bot[i][j][2];
		}
	}
	/* exchange back-front and front-back */
	for (i=1;i<=NX[0];i++){
		for (k=1;k<=NX[2];k++){
			e[i][0][k].xx = rbuf_fro_to_bac[i][k][1];
			e[i][0][k].yy = rbuf_fro_to_bac[i][k][2];
			e[i][0][k].zz = rbuf_fro_to_bac[i][k][3];
			e[i][0][k].xz = rbuf_fro_to_bac[i][k][4];

			e[i][NX[1]+1][k].yz = rbuf_bac_to_fro[i][k][1];
			e[i][NX[1]+1][k].xy = rbuf_bac_to_fro[i][k][2];
		}
	}
	/* exchange left-right and right-left */
	for (j=1;j<=NX[1];j++){
		for (k=1;k<=NX[2];k++){
			e[0][j][k].xx = rbuf_rig_to_lef[j][k][1];
			e[0][j][k].yy = rbuf_rig_to_lef[j][k][2];
			e[0][j][k].zz = rbuf_rig_to_lef[j][k][3];
			e[0][j][k].yz = rbuf_rig_to_lef[j][k][4];

			e[NX[0]+1][j][k].xz = rbuf_lef_to_rig[j][k][1];
			e[NX[0]+1][j][k].xy = rbuf_lef_to_rig[j][k][2];
		}
	}



	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for strain exchange: \t\t %4.2f s.\n",time);

	return time;
}
