/*------------------------------------------------------------------------
 *   copy material parameters from bottom / right boundary to 
 *   top / left boundary before averaging.
 *
 *   T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void matcopy(float *** rho){

	/* extern variables */
	extern int		NX[3], INDEX[6];
	extern const int	TAG2, TAG4, TAG6;
	extern MPI_Comm		COMM_CART;

	/* MPI variables */
	MPI_Status	status[2];
	MPI_Request	request[2];

	/* local variables */
	int	i, j, k;
	float	** sbuf_rig_to_lef, ** sbuf_fro_to_bac, ** sbuf_bot_to_top;
	float	** rbuf_rig_to_lef, ** rbuf_fro_to_bac, ** rbuf_bot_to_top;


	/* memory allocation for buffer variables */
	sbuf_rig_to_lef = matrix(0,NX[1],0,NX[2]);
	rbuf_rig_to_lef = matrix(0,NX[1],0,NX[2]);
	sbuf_fro_to_bac = matrix(0,NX[0],0,NX[2]);
	rbuf_fro_to_bac = matrix(0,NX[0],0,NX[2]);
	sbuf_bot_to_top = matrix(0,NX[0],0,NX[1]);
	rbuf_bot_to_top = matrix(0,NX[0],0,NX[1]);

	/* exchange of material parameters at right / left edge of global grid */
	for (j=0;j<=NX[1];j++){
		for (k=0;k<=NX[2];k++){
			/* storage of right edge of local volume into buffer */
			sbuf_rig_to_lef[j][k] = rho[NX[0]][j][k];
		}
	}

	MPI_Irecv(&rbuf_rig_to_lef[0][0],(NX[1]+1)*(NX[2]+1),MPI_FLOAT,INDEX[0],TAG2,COMM_CART,&request[0]);
	MPI_Isend(&sbuf_rig_to_lef[0][0],(NX[1]+1)*(NX[2]+1),MPI_FLOAT,INDEX[1],TAG2,COMM_CART,&request[1]);
	
	MPI_Waitall(2,&request[0], &status[0]);
	
	/* exchange of material parameters at right / left edge of global grid */
	for (j=0;j<=NX[1];j++){
		for (k=0;k<=NX[2];k++){
			rho[0][j][k] = rbuf_rig_to_lef[j][k];
		}
	}



	/* exchange of material parameters at front / back edge of global grid */
	for (i=0;i<=NX[0];i++){
		for (k=0;k<=NX[2];k++){
			/* storage of right edge of local volume into buffer */
			sbuf_fro_to_bac[i][k] = rho[i][NX[1]][k];
		}
	}

	MPI_Irecv(&rbuf_fro_to_bac[0][0],(NX[0]+1)*(NX[2]+1),MPI_FLOAT,INDEX[2],TAG4,COMM_CART,&request[0]);
	MPI_Isend(&sbuf_fro_to_bac[0][0],(NX[0]+1)*(NX[2]+1),MPI_FLOAT,INDEX[3],TAG4,COMM_CART,&request[1]);

	MPI_Waitall(2,&request[0], &status[0]);

	/* exchange of material parameters at front / back edge of global grid */
	for (i=0;i<=NX[0];i++){
		for (k=0;k<=NX[2];k++){
			rho[i][0][k] = rbuf_fro_to_bac[i][k];
		}
	}



	/* exchange of material parameters at bottom / top of global grid */
	for (i=0;i<=NX[0];i++){
		for (j=0;j<=NX[1];j++){
			/* storage of bottom of local volume into buffer */
			sbuf_bot_to_top[i][j] = rho[i][j][NX[2]];
		}
	}

	MPI_Irecv(&rbuf_bot_to_top[0][0],(NX[0]+1)*(NX[1]+1),MPI_FLOAT,INDEX[4],TAG6,COMM_CART,&request[0]);
	MPI_Isend(&sbuf_bot_to_top[0][0],(NX[0]+1)*(NX[1]+1),MPI_FLOAT,INDEX[5],TAG6,COMM_CART,&request[1]);

	MPI_Waitall(2,&request[0], &status[0]);

	/* exchange of material parameters at bottom / top of global grid */
	for (i=0;i<=NX[0];i++){
		for (j=0;j<=NX[1];j++){
			rho[i][j][0] = rbuf_bot_to_top[i][j];
		}
	}

	/* free buffer variables */
	free_matrix(sbuf_rig_to_lef,0,NX[1],0,NX[2]);
	free_matrix(rbuf_rig_to_lef,0,NX[1],0,NX[2]);
	free_matrix(sbuf_fro_to_bac,0,NX[0],0,NX[2]);
	free_matrix(rbuf_fro_to_bac,0,NX[0],0,NX[2]);
	free_matrix(sbuf_bot_to_top,0,NX[0],0,NX[1]);
	free_matrix(rbuf_bot_to_top,0,NX[0],0,NX[1]);

}
