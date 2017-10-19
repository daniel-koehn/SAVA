/*------------------------------------------------------------------------
 *   Perform arithmetic or harmonic averaging of wavefields
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void av_wavefield(float *** u, float *** ua, int cc, 
float * wx1, float * wx2, float * wy1, float * wy2, float * wz1, float * wz2){

	/* extern variables */
	extern int	NX[3];

	/* local variables */
	int	i, j, k;
	float	wa, wb, wc, wd;
      
        
	/* directions of averaging are specified by cc:
		c=1: right, i.e. (i)     and (i+1)   -> i+1/2 
		c=2: left, i.e.  (i-1/2) and (i+1/2) -> i
		c=3: front, i.e. (j)     and (j+1)   -> j+1/2 
		c=4: back, i.e.  (j-1/2) and (j+1/2) -> j
		c=5: down, i.e.  (k)     and (k+1)   -> k+1/2 
		c=6: up, i.e.    (k-1/2) and (k+1/2) -> k
		c=7: back, up, i.e.     (j-1/2,k-1/2), (j-1/2,k+1/2), (j+1/2,k+1/2) and (j+1/2,k-1/2) -> (j,k)
		c=8: left, up, i.e.     (i-1/2,k-1/2), (i-1/2,k+1/2), (i+1/2,k+1/2) and (i+1/2,k-1/2) -> (i,k)
		c=9: left, back, i.e.   (i-1/2,j-1/2), (i-1/2,j+1/2), (i+1/2,j+1/2) and (i+1/2,j-1/2) -> (i,j)
		c=10: front, down, i.e. (j,k), (j,k+1), (j+1,k+1) and (j+1,k) -> (j+1/2,k+1/2) 
		c=11: right, down, i.e. (i,k), (i,k+1), (i+1,k+1) and (i+1,k) -> (i+1/2,k+1/2) 	
		c=12: right, front, i.e.(i,j), (i,j+1), (i+1,j+1) and (i+1,j) -> (i+1/2,j+1/2)
		c=13: back, down, i.e.  (j-1/2,k), (j-1/2,k+1), (j+1/2,k+1) and (j+1/2,k) -> (j,k+1/2)
		c=14: left, down, i.e.  (i-1/2,k), (i-1/2,k+1), (i+1/2,k+1) and (i+1/2,k) -> (i,k+1/2)
		c=15: left, front, i.e. (i-1/2,j), (i-1/2,j+1), (i+1/2,j+1) and (i+1/2,j) -> (i,j+1/2)
		c=16: front, up, i.e.   (j,k-1/2), (j,k+1/2), (j+1,k+1/2) and (j+1,k-1/2) -> (j+1/2,k) 
		c=17: right, up, i.e.   (i,k-1/2), (i,k+1/2), (i+1,k+1/2) and (i+1,k-1/2) -> (i+1/2,k)
		c=18: right, back, i.e. (i,j-1/2), (i,j+1/2), (i+1,j+1/2) and (i+1,j-1/2) -> (i+1/2,j)		
	*/

	/* weighted arithmetic averaging using weights computed with fd_coeff.c */
	/* weights computed outside of loops and saved in vectors (faster) */
	switch (cc){
		case 1 : 	/* right -> wx */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = (wa*u[i][j][k] + wb*u[i+1][j][k]);
					}
				}
			} 
			break;
		case 2 :	/* left -> wxp */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
			        		ua[i][j][k] = (wa*u[i-1][j][k] + wb*u[i][j][k]);
					}
				}
			}
			break;
		case 3 :	/* front -> wy */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					wa = wy1[j];
					wb = wy2[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = (wa*u[i][j][k] + wb*u[i][j+1][k]);
					}
				}
			}
			break;
		case 4 :	/* back -> wyp */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					wa = wy1[j];
					wb = wy2[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = (wa*u[i][j-1][k] + wb*u[i][j][k]);
					}
				}
			}
			break;	
		case 5 :	/* down -> wz */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = (wz1[k]*u[i][j][k] + wz2[k]*u[i][j][k+1]);
					}
				}
			}
			break;
		case 6 :	/* up -> wzp */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = (wz1[k]*u[i][j][k-1] + wz2[k]*u[i][j][k]);
					}
				}
			}
			break;
		case 7 :	/* back and up -> wyp, wzp */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					wa = wy1[j];
					wb = wy2[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i][j][k];
						ua[i][j][k] += wz1[k]*wb*u[i][j+1][k];
						ua[i][j][k] += wz2[k]*wb*u[i][j+1][k+1];
						ua[i][j][k] += wz2[k]*wa*u[i][j][k+1];
					}
				}
			}
			break;
		case 8 :	/* left and up -> wxp, wzp */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i][j][k];
						ua[i][j][k] += wz1[k]*wb*u[i+1][j][k];
						ua[i][j][k] += wz2[k]*wb*u[i+1][j][k+1];
						ua[i][j][k] += wz2[k]*wa*u[i][j][k+1];
					}
				}
			}
			break;
		case 9 :	/* left and back -> wxp, wyp */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					wc  = wy2[j]*wb;
					wd  = wy2[j]*wa;
					wa *= wy1[j];
					wb *= wy1[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wa*u[i][j][k];
						ua[i][j][k] += wb*u[i+1][j][k];
						ua[i][j][k] += wc*u[i+1][j+1][k];
						ua[i][j][k] += wd*u[i][j+1][k];
					}
				}
			}
			break;
		case 10 :	/* front and down -> wy, wz */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					wa = wy1[j];
					wb = wy2[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i][j-1][k-1];
						ua[i][j][k] += wz1[k]*wb*u[i][j][k-1];
						ua[i][j][k] += wz2[k]*wb*u[i][j][k];
						ua[i][j][k] += wz2[k]*wa*u[i][j-1][k];
					}
				}
			}
			break;
		case 11 :	/* right and down -> wx, wz */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i-1][j][k-1];
						ua[i][j][k] += wz1[k]*wb*u[i][j][k-1];
						ua[i][j][k] += wz2[k]*wb*u[i][j][k];
						ua[i][j][k] += wz2[k]*wa*u[i-1][j][k];
					}
				}
			}
			break;
		case 12 :	/* right and front -> wx, wy */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					wc  = wy2[j]*wb;
					wd  = wy2[j]*wa;
					wa *= wy1[j];
					wb *= wy1[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wa*u[i-1][j-1][k];
						ua[i][j][k] += wb*u[i][j-1][k];
						ua[i][j][k] += wc*u[i][j][k];
						ua[i][j][k] += wd*u[i-1][j][k];
					}
				}
			}
			break;
		case 13 :	/* back and down -> wyp, wz */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					wa = wy1[j];
					wb = wy2[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i][j][k-1];
						ua[i][j][k] += wz1[k]*wb*u[i][j+1][k-1];
						ua[i][j][k] += wz2[k]*wb*u[i][j+1][k];
						ua[i][j][k] += wz2[k]*wa*u[i][j][k];
					}
				}
			}
			break;
		case 14 :	/* left and down -> wxp, wz */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i][j][k-1];
						ua[i][j][k] += wz1[k]*wb*u[i+1][j][k-1];
						ua[i][j][k] += wz2[k]*wb*u[i+1][j][k];
						ua[i][j][k] += wz2[k]*wa*u[i][j][k];
					}
				}
			}
			break;
		case 15 :	/* left and front -> wxp, wy */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					wc  = wy2[j]*wb;
					wd  = wy2[j]*wa;
					wa *= wy1[j];
					wb *= wy1[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wa*u[i][j-1][k];
						ua[i][j][k] += wb*u[i+1][j-1][k];
						ua[i][j][k] += wc*u[i+1][j][k];
						ua[i][j][k] += wd*u[i][j][k];
					}
				}
			}
			break;
		case 16 :	/* front and up -> wy, wzp*/
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					wa = wy1[j];
					wb = wy2[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i][j-1][k];
						ua[i][j][k] += wz1[k]*wb*u[i][j][k];
						ua[i][j][k] += wz2[k]*wb*u[i][j][k+1];
						ua[i][j][k] += wz2[k]*wa*u[i][j-1][k+1];
					}
				}
			}
			break;
		case 17 :	/* right and up -> wx, wzp */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wz1[k]*wa*u[i-1][j][k];
						ua[i][j][k] += wz1[k]*wb*u[i][j][k];
						ua[i][j][k] += wz2[k]*wb*u[i][j][k+1];
						ua[i][j][k] += wz2[k]*wa*u[i-1][j][k+1];
					}
				}
			}
			break;
		case 18 :	/* right and back -> wx, wyp */
			for (i=1;i<=NX[0];i++){
				wa = wx1[i];
				wb = wx2[i];
				for (j=1;j<=NX[1];j++){
					wc  = wy2[j]*wb;
					wd  = wy2[j]*wa;
					wa *= wy1[j];
					wb *= wy1[j];
					for (k=1;k<=NX[2];k++){
						ua[i][j][k]  = wa*u[i-1][j][k];
						ua[i][j][k] += wb*u[i][j][k];
						ua[i][j][k] += wc*u[i][j+1][k];
						ua[i][j][k] += wd*u[i-1][j+1][k];
					}
				}
			}
			break;
	}
	

}
