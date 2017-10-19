/*-------------------------------------------------------------
 *   Check FD-Grid for stability and grid dispersion.
 *   If the stability criterion is not fullfilled the program will
 *   terminate.
 *
 *   T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------*/


#include "fd.h"

void checkfd(FILE *fp, float *** rho, float *** pi, float *** taup, float *eta, float *x, float *y, float *z){

	/* extern variables */
	extern float	DT, TS;
        extern int	NX[3], L, MYID;

	/* local variables */
	float		vpmax, vpmin;
	float		dx, dy, dz;
	float		vmax, vmin;
	float		dhmax;
	float		vpmax_loc=0.0, vpmin_loc=1e9;
	float		dxmax_loc=0.0, dxmin_loc=1e9;
	float		dymax_loc=0.0, dymin_loc=1e9;
	float		dzmax_loc=0.0, dzmin_loc=1e9;
	float		vmax_loc,  vmin_loc;
	float		dhmax_loc, dhmin_loc;
	float		vmax_glob, vmin_glob;
	float		dhmax_glob,dhmin_glob;
	float		dtstab, dhstab;
	float		dtstab_loc=1e9, dhstab_loc=1e9;
	float		dtstab_glob, dhstab_glob;
	float		sum, ws;
	const float	w    = 2.0*PI*DT/TS;	/* center frequency of source */
	const float	fmax = 3.0/TS;		/* max. frequency */
	int		i, j, k, l;


	fprintf(fp,"\n **Message from checkfd (printed by PE %d):\n",MYID);

	sum = 0.0;
	for (l=1;l<=L;l++){
		ws   = eta[l]/w;
		sum += (1.0/(1.0+ws*ws));
	}

	for (i=1;i<=NX[0];i++){
		for (j=1;j<=NX[1];j++){
			for (k=1;k<=NX[2];k++){
				/* maximum phase velocity of P-waves at infinite frequency */
				vpmax = sqrt((pi[i][j][k]/rho[i][j][k])*(1.0+L*taup[i][j][k]));

				/* minimum phase velocity of P-waves at center frequency of the source */
				vpmin = sqrt((pi[i][j][k]/rho[i][j][k])*(1.0+sum*taup[i][j][k]));

				/* min. and max. grid spacing */
				dx = x[i+1]-x[i];
				dy = y[j+1]-y[j];
				dz = z[k+1]-z[k];

				vmax  = vpmax;
				vmin  = vpmin;
				dhmax = max(dx,dy);
				dhmax = max(dhmax,dz);

				vpmax_loc = max(vpmax_loc,vpmax);
				vpmin_loc = min(vpmin_loc,vpmin);
				dxmax_loc = max(dxmax_loc,dx);
				dxmin_loc = min(dxmin_loc,dx);
				dymax_loc = max(dymax_loc,dy);
				dymin_loc = min(dymin_loc,dy);
				dzmax_loc = max(dzmax_loc,dz);
				dzmin_loc = min(dzmin_loc,dz);

				dtstab     = 1.0/(vmax*sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));
				dtstab_loc = min(dtstab_loc,dtstab);
				if (vmin != 0.0){
					dhstab     = vmin/(10.0*fmax*dhmax);
					dhstab_loc = min(dhstab_loc,dhstab);
				}
			}
		}
	}
	vmax_loc  = vpmax_loc;
	vmin_loc  = vpmin_loc;
	dhmax_loc = max(dxmax_loc,dymax_loc);
	dhmax_loc = max(dhmax_loc,dzmax_loc);
	dhmin_loc = min(dxmin_loc,dymin_loc);
	dhmin_loc = min(dhmin_loc,dzmin_loc);

	fprintf(fp," Minimum and maximum P-wave velocity within subvolumes (MYID=%d): \n",MYID);
	fprintf(fp," vP_max(f=inf)=%e m/s\t vP_min(f=fc)=%e m/s\n\n", vpmax_loc,vpmin_loc);

	fprintf(fp," Minimum and maximum grid spacing within subvolumes (MYID=%d): \n",MYID);
	fprintf(fp," DX_max=%e m\t DX_min=%e m\n",  dxmax_loc,dxmin_loc);
	fprintf(fp," DY_max=%e m\t DY_min=%e m\n",  dymax_loc,dymin_loc);
	fprintf(fp," DZ_max=%e m\t DZ_min=%e m\n\n",dzmax_loc,dzmin_loc);
	MPI_Barrier(MPI_COMM_WORLD);

	/* find global velocity maximum and minimum */
	MPI_Allreduce(&vmax_loc,&vmax_glob,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&vmin_loc,&vmin_glob,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);

	/* find global maximum and global minimum grid spacing */
	MPI_Allreduce(&dhmax_loc,&dhmax_glob,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&dhmin_loc,&dhmin_glob,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);

	/* find out whether grid dispersion occurs and whether simulation is stable */
	MPI_Allreduce(&dhstab_loc,&dhstab_glob,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	MPI_Allreduce(&dtstab_loc,&dtstab_glob,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);

	if (!(MYID)){
		fprintf(fp," Global values for entire model: \n");
		fprintf(fp," V_max = %e m/s \t V_min =%e m/s \n",  vmax_glob, vmin_glob);
		fprintf(fp," DH_max= %e m   \t DH_min=%e m   \n\n",dhmax_glob,dhmin_glob);
		fprintf(fp," ------------------ CHECK FOR GRID DISPERSION --------------------\n");
		fprintf(fp," To satisfactorily limit grid dispersion, the number of gridpoints \n");
		fprintf(fp," per minimum wavelength should be 10 (better more).\n");
		fprintf(fp," Here the minimum wavelength is assumed to be minimum model phase velocity \n");
		fprintf(fp," at maximum frequency of the source\n");
		fprintf(fp," devided by maximum frequency of the source.\n");
		fprintf(fp," The maximum frequency of the source signal (approx. 3 times the \n");
		fprintf(fp," center frequency of the signal defined in SOURCE_FILE) is %8.2f Hz.\n",fmax);
		fprintf(fp," The minimum wavelength in the following simulation will\n");
		fprintf(fp," be %e meter.\n", vmin_glob/fmax);
		fprintf(fp," Thus, the recommended maximum value for the spatial increments in the \n");
		fprintf(fp," x-, y- and z-direction is %e meter.\n", vmin_glob/fmax/10.0);
		fprintf(fp," You have specified a maximum spatial increment of DH= %e meter \n\n", dhmax_glob);

		if (1 > dhstab_glob)
			warning(" Grid dispersion will influence wave propagation, choose smaller grid spacing.");

		/* gamma for second order spatial FD operators */
/*		switch(M){
		case 0: gamma=1.15; break;
		case 1: gamma=1.14; break;
		case 2: gamma=1.34; break;
		case 3: gamma=1.15; break;
		case 4: gamma=2.37; break;
		case 5: gamma=1.15; break;
		case 6: gamma=3.35; break;
		default: gamma=1.15;
			 warning(" Assuming gamma=1.15 for stability limit, DT may be choosen smaller.");
		         break;
		}
*/

		/*gamma=sqrt(1.0+(M*M/8.0)); gamma=2.37;*/
		fprintf(fp," ----------------------- CHECK FOR STABILITY ---------------------\n");
		fprintf(fp," The following simulation is stable provided that\n\n");
		fprintf(fp," \t DT < 1/(v_max*sqrt(1/(DX*DX)+1/(DY*DY)+1/(DZ*DZ)))\n\n");
		fprintf(fp," where DT is the time step, DH is the grid spacing and \n");
		fprintf(fp," v_max is the maximum phase velocity at infinite frequency.\n\n");
		/*fprintf(fp," and gamma = %4.2f .\n", gamma);*/

		fprintf(fp," You have specified DT= %e s.\n", DT);
		fprintf(fp," In this simulation the stability limit for timestep DT is %e s.\n",dtstab_glob);

		/*dtstab=dh/(gamma*sqrt(2.0)*cmax); dtstab*=1.0/sqrt(1.0+M*M*dh/8.0);*/
		if (DT > dtstab_glob)
			warning(" The simulation MAY get unstable, choose smaller DT! \n");
		else 
			fprintf(fp," The simulation will be stable.\n");

	}
}

