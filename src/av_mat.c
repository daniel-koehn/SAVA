/*------------------------------------------------------------------------
 *   Perform arithmetic or harmonic averaging of material properties
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void av_mat(float *** u, float *** ua, int cc, char t, float * x, float * xp, float * y, float * yp, float * z, float * zp){

	/* extern variables */
	extern int	NX[3];
	extern FILE	*FP;

	/* local variables */
	int	i, j, k;
	float	a, b, c, d;
	float	dx1, dx2, dy1, dy2, dz1, dz2;

	/* directions of averaging are specified by c:
		c=1: right, i.e. i+1/2
		c=2: left, i.e. i-1/2
		c=3: front, i.e. j+1/2
		c=4: back, i.e. j-1/2
		c=5: down, i.e. k+1/2
		c=6: up, i.e. k-1/2
		c=7: center left, up, i.e. i-1/2,k-1/2 
		c=8: center back, left, i.e. j-1/2,i-1/2 
		c=9: center up, back, i.e. k-1/2,j-1/2
	*/


	/* type of averaging is determined by t 
		t=h: harmonic averaging
		t=a: arithmetic averaging
	*/


	fprintf(FP," averaging of material parameters: c=%d \t t=%c \n", cc,t);

	if (t=='a'){ /* arithmetic averaging */
		switch (cc){
		case 1 : 	/* right */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i+1]-xp[i])/(xp[i+1]-xp[i]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = dx1*u[i][j][k] + dx2*u[i+1][j][k];
					}
				}
			} 
			break;
		case 2 :	/* left */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i]-xp[i-1])/(xp[i]-xp[i-1]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = dx1*u[i-1][j][k] + dx2*u[i][j][k];
					}
				}
			}
			break;
		case 3 :	/* front */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j+1]-yp[j])/(yp[j+1]-yp[j]);
					dy2 = 1.0 - dy1;
					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = dy1*u[i][j][k] + dy2*u[i][j+1][k];
					}
				}
			}
			break;
		case 4 :	/* back */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j]-yp[j-1])/(yp[j]-yp[j-1]);
					dy2 = 1.0 - dy1;
					for (k=1;k<=NX[2];k++){
		        			ua[i][j][k] = dy1*u[i][j-1][k] + dy2*u[i][j][k];
					}
				}
			}
			break;
		case 5 :	/* down */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						dz1 = (z[k+1]-zp[k])/(zp[k+1]-zp[k]);
						dz2 = 1.0 - dz1;
						ua[i][j][k] = dz1*u[i][j][k] + dz2*u[i][j][k+1];
					}
				}
			}
			break;
		case 6 :	/* up */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						dz1 = (z[k]-zp[k-1])/(zp[k]-zp[k-1]);
						dz2 = 1.0 - dz1;
						ua[i][j][k] = dz1*u[i][j][k-1] + dz2*u[i][j][k];
					}
				}
			}
			break;
		case 7 :	/* up and left */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i]-xp[i-1])/(xp[i]-xp[i-1]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						dz1 = (z[k]-zp[k-1])/(zp[k]-zp[k-1]);
						dz2 = 1.0 - dz1;

						a = dx1*dz1;
						b = dx2*dz1;
						c = dx1*dz2;
						d = dx2*dz2;

						ua[i][j][k] = a*u[i-1][j][k-1]+b*u[i][j][k-1]+c*u[i-1][j][k]+d*u[i][j][k];
					}
				}
			}
			break;
		case 8 :	/* left and back */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i]-xp[i-1])/(xp[i]-xp[i-1]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j]-yp[j-1])/(yp[j]-yp[j-1]);
					dy2 = 1.0 - dy1;

					a = dx1*dy1;
					b = dx2*dy1;
					c = dx1*dy2;
					d = dx2*dy2;

					for (k=1;k<=NX[2];k++){
						ua[i][j][k] = a*u[i-1][j-1][k]+b*u[i][j-1][k]+c*u[i-1][j][k]+d*u[i][j][k];
					}
				}
			}
			break;
		case 9 :	/* back and up */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j]-yp[j-1])/(yp[j]-yp[j-1]);
					dy2 = 1.0 - dy1;
					for (k=1;k<=NX[2];k++){
						dz1 = (z[k]-zp[k-1])/(zp[k]-zp[k-1]);
						dz2 = 1.0 - dz1;

						a = dy1*dz1;
						b = dy2*dz1;
						c = dy1*dz2;
						d = dy2*dz2;

						ua[i][j][k] = a*u[i][j-1][k-1]+b*u[i][j][k-1]+c*u[i][j-1][k]+d*u[i][j][k];
					}
				}
			}
			break;
		}
	}

	if (t=='h'){ /* harmonic averaging */
		switch (cc){
		case 1 : 	/* right */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i+1]-xp[i])/(xp[i+1]-xp[i]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						if ((u[i][j][k]==0.0) || (u[i+1][j][k]==0.0)) 
							ua[i][j][k] = 0.0;
						else
							ua[i][j][k] = u[i][j][k]*u[i+1][j][k]/(dx2*u[i][j][k] + dx1*u[i+1][j][k]);
					}
				}
			}
			break;
		case 2 :	/* left */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i]-xp[i-1])/(xp[i]-xp[i-1]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						if ((u[i][j][k]==0.0) || (u[i-1][j][k]==0.0)) 
							ua[i][j][k] = 0.0;
						else
							ua[i][j][k] = u[i][j][k]*u[i-1][j][k]/(dx1*u[i][j][k] + dx2*u[i-1][j][k]);
					}
				}
			}
			break;
		case 3 : 	/* front */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j+1]-yp[j])/(yp[j+1]-yp[j]);
					dy2 = 1.0 - dy1;
					for (k=1;k<=NX[2];k++){
						if ((u[i][j][k]==0.0) || (u[i][j+1][k]==0.0)) 
							ua[i][j][k] = 0.0;
						else
							ua[i][j][k] = u[i][j][k]*u[i][j+1][k]/(dy2*u[i][j][k] + dy1*u[i][j+1][k]);
					}
				}
			}
			break;
		case 4 :	/* back */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j]-yp[j-1])/(yp[j]-yp[j-1]);
					dy2 = 1.0 - dy1;
					for (k=1;k<=NX[2];k++){
						if ((u[i][j][k]==0.0) || (u[i][j-1][k]==0.0)) 
							ua[i][j][k] = 0.0;
						else
							ua[i][j][k] = u[i][j][k]*u[i][j-1][k]/(dy1*u[i][j][k] + dy2*u[i][j-1][k]);
					}
				}
			} 
			break;
		case 5 :	/* down */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						if ((u[i][j][k]==0.0) || (u[i][j][k+1]==0.0)) 
							ua[i][j][k] = 0.0;
						else{
							dz1 = (z[k+1]-zp[k])/(zp[k+1]-zp[k]);
							dz2 = 1.0 - dz1;

							ua[i][j][k] = u[i][j][k]*u[i][j][k+1]/(dz2*u[i][j][k] + dz1*u[i][j][k+1]);
						}
					}
				}
			}
			break;
		case 6 :	/* up */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						if ((u[i][j][k]==0.0) || (u[i][j][k-1]==0.0))
							ua[i][j][k] = 0.0;
						else{
							dz1 = (z[k]-zp[k-1])/(zp[k]-zp[k-1]);
							dz2 = 1.0 - dz1;

							ua[i][j][k] = u[i][j][k]*u[i][j][k-1]/(dz1*u[i][j][k] + dz2*u[i][j][k-1]);
						}
					}
				}
			}
			break;
		case 7 :	/* up and left */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i]-xp[i-1])/(xp[i]-xp[i-1]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					for (k=1;k<=NX[2];k++){
						if ((u[i-1][j][k-1]==0.0) || (u[i][j][k-1]==0.0) || (u[i-1][j][k]==0.0) || (u[i][j][k]==0.0)){
							ua[i][j][k] = 0.0;
						}
						else{
							dz1 = (z[k]-zp[k-1])/(zp[k]-zp[k-1]);
							dz2 = 1.0 - dz1;

							a = dx1*dz1;
							b = dx2*dz1;
							c = dx1*dz2;
							d = dx2*dz2;

							ua[i][j][k] = 1.0/(a/u[i-1][j][k-1] + b/u[i][j][k-1] + c/u[i-1][j][k] + d/u[i][j][k]);
						}
					}
				}
			}
			break;
		case 8 :	/* left and back */
			for (i=1;i<=NX[0];i++){
				dx1 = (x[i]-xp[i-1])/(xp[i]-xp[i-1]);
				dx2 = 1.0 - dx1;
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j]-yp[j-1])/(yp[j]-yp[j-1]);
					dy2 = 1.0 - dy1;

					a = dx1*dy1;
					b = dx2*dy1;
					c = dx1*dy2;
					d = dx2*dy2;

					for (k=1;k<=NX[2];k++){
						if ((u[i-1][j-1][k]==0.0) || (u[i][j-1][k]==0.0) || (u[i-1][j][k]==0.0) || (u[i][j][k]==0.0)){
							ua[i][j][k] = 0.0;
						}
						else{
							ua[i][j][k] = 1.0/(a/u[i-1][j-1][k] + b/u[i][j-1][k] + c/u[i-1][j][k] + d/u[i][j][k]);
						}
					}
				}
			}
			break;
		case 9 :	/* back and up */
			for (i=1;i<=NX[0];i++){
				for (j=1;j<=NX[1];j++){
					dy1 = (y[j]-yp[j-1])/(yp[j]-yp[j-1]);
					dy2 = 1.0 - dy1;
					for (k=1;k<=NX[2];k++){
						if ((u[i][j-1][k-1]==0.0) || (u[i][j][k-1]==0.0) || (u[i][j-1][k]==0.0) || (u[i][j][k]==0.0)){
							ua[i][j][k] = 0.0;
						}
						else{
							dz1 = (z[k]-zp[k-1])/(zp[k]-zp[k-1]);
							dz2 = 1.0 - dz1;

							a = dy1*dz1;
							b = dy2*dz1;
							c = dy1*dz2;
							d = dy2*dz2;

							ua[i][j][k] = 1.0/(a/u[i][j-1][k-1] + b/u[i][j][k-1] + c/u[i][j-1][k] + d/u[i][j][k]);
						}
					}
				}
			}
			break;
		}
	}


}
