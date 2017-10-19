/*------------------------------------------------------------------------
 *   Check FD-Grid for stability and grid dispersion.
 *   If the stability criterion is not fullfilled the program will
 *   terminate.
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/


#include "fd.h"
#include "fd3daniso.h"

void checkfd(FILE *fp, float *** rho,
	float *** c1111, float *** c1122, float *** c1133, float *** c1123, float *** c1113, float *** c1112, 
			 float *** c2222, float *** c2233, float *** c2223, float *** c1322, float *** c1222, 
					  float *** c3333, float *** c2333, float *** c1333, float *** c1233, 
							   float *** c2323, float *** c1323, float *** c1223, 
									    float *** c1313, float *** c1213, 
											     float *** c1212,
	float *** taus, float *** taup, float * eta, float * x, float * y, float * z){


	/* extern variables */
	extern float	DT, TS;
        extern int	NX[3], L, MYID;

	/* local variables */
	float	dx, dy, dz;
	float	vmax, vmin;
	float	vmin1, vmax1;
	float	dhmax, dhmin;
	float	vmax_loc =0.0, vmin_loc =1e9;
	float	dxmax_loc=0.0, dxmin_loc=1e9;
	float	dymax_loc=0.0, dymin_loc=1e9;
	float	dzmax_loc=0.0, dzmin_loc=1e9;
	float	dhmax_loc, dhmin_loc;
	float	vmax_glob, vmin_glob;
	float	dhmax_glob,dhmin_glob;
	float	dtstab, dhstab;
	float	dtstab_loc=1e9, dhstab_loc=1e9;
	float	dtstab_glob,    dhstab_glob;
	float	sum, ws;
//	float	** g, * e;
	float	* g, * e;
	float	l11, l22, l33, l23, l13, l12;
	float	g11, g22, g33, g23, g13, g12;
	float	c11, c22, c33, c23, c13, c12;
	float	phi, theta;
	float	dphi, dtheta;
	float	*sphi, *cphi, *stheta, *ctheta;

	const int	nphi=18, ntheta=18;
	const float	w    = 2.0*PI*DT/TS;	/* center frequency of source */
	const float	fmax = 3.0/TS;		/* max. frequency */

	int	i, j, k, l, m;


	fprintf(fp,"\n **Message from checkfd (printed by PE %d):\n",MYID);

	/* memory allocation */
//	g = matrix(1,3,1,3);
	g = vector(1,6);
	e = vector(1,3);

	sum = 0.0;
	for (l=1;l<=L;l++){
		ws  = eta[l]/w;
		sum = sum+(1.0/(1.0+ws*ws));
	}

	/* vectors for angles */
	sphi   = vector(1,nphi);
	cphi   = vector(1,nphi);
	stheta = vector(1,ntheta);
	ctheta = vector(1,ntheta);

	/* 0<phi<PI, 0<theta<=PI */
	dphi   = PI/((float)(nphi));
	dtheta = PI/((float)(ntheta));

	/* compute angles */
	for (l=1;l<nphi;l++){
		/* inclination */
		phi = (float)(l)*dphi;

		cphi[l] = cos(phi);
		cphi[l] = cphi[l]*cphi[l];
		sphi[l] = 1.0-cphi[l];
	}
	for (m=1;m<=ntheta;m++){
		/* azimuth */
		theta = (float)(m)*dtheta;

		ctheta[m] = cos(theta);
		ctheta[m] = ctheta[m]*ctheta[m];
		stheta[m] = 1.0-ctheta[m];
	}


	/* loop over grid points */
	for (i=1;i<=NX[0];i++){
		fprintf(fp," ... check grid points %d of %d (printed by PE %d)\n",i,NX[0],MYID);
		for (j=1;j<=NX[1];j++){
			for (k=1;k<=NX[2];k++){
				/* find minimum and maximum velocities */

				/* scan direction of wave propagation phi=0 */
				/* (l1==l2=0, l3=1, l23=l13=l12=0) */
				/* Kelvin-Christoffel matrix */
//				g[1][1] = c1313[i][j][k];
//				g[2][2] = c2323[i][j][k];
//				g[3][3] = c3333[i][j][k];
				g[1] = c1313[i][j][k];
				g[2] = c2323[i][j][k];
				g[3] = c3333[i][j][k];
  
//				g[2][3] = c2333[i][j][k];
//				g[1][3] = c1333[i][j][k];
//				g[1][2] = c1323[i][j][k];
				g[4] = c2333[i][j][k];
				g[5] = c1333[i][j][k];
				g[6] = c1323[i][j][k];
  
//				g[3][2] = g[2][3];
//				g[3][1] = g[1][3];
//				g[2][1] = g[1][2];

				/* min./max. eigenvalues */
//				eigenvalue(g,e,3);
				eigenvalue(g,e);

				/* maximum model phase velocity at infinite frequency */
				vmax = sqrt((e[1]/rho[i][j][k])*(1.0+L*taup[i][j][k]));

				/* minimum model phase velocity at center frequency of the source */
				vmin = sqrt((e[3]/rho[i][j][k])*(1.0+sum*taus[i][j][k]));

				/* combined material parameters for Kelvin-Christoffel matrix */
				c11 = c1123[i][j][k] + c1213[i][j][k];
				c22 = c1223[i][j][k] + c1322[i][j][k];
				c33 = c1323[i][j][k] + c1233[i][j][k];
				c23 = c2323[i][j][k] + c2233[i][j][k];
				c13 = c1133[i][j][k] + c1313[i][j][k];
				c12 = c1122[i][j][k] + c1212[i][j][k];

				/* scan different directions of wave propagation (phi>0) */
				for (l=1;l<nphi;l++){
					/* inclination */
					l33  = cphi[l];	/* cos(phi)^2 */

					/* Kelvin-Christoffel matrix (l3-terms) */
					g11 = c1313[i][j][k]*l33;
					g22 = c2323[i][j][k]*l33;
					g33 = c3333[i][j][k]*l33;
					g23 = c2333[i][j][k]*l33;
					g13 = c1333[i][j][k]*l33;
					g12 = c1323[i][j][k]*l33;

					for (m=1;m<=ntheta;m++){
						/* azimuth */
						l11 = sphi[l]*ctheta[m];	/* sin(phi)^2*cos(theta)^2 */
						l22 = sphi[l]*stheta[m];	/* sin(phi)^2*sin(theta)^2 */
						l23 = sqrt(l22*l33);		/* sin(phi)*cos(phi)*sin(theta) */
						l13 = sqrt(l11*l33);		/* sin(phi)*cos(phi)*cos(theta) */
						l12 = sqrt(l11*l22);		/* sin(phi)^2*sin(theta)*cos(theta) */

						/* Kelvin-Christoffel matrix */
//						g[1][1] = c1111[i][j][k]*l11 + c1212[i][j][k]*l22 + g11 + 2.0*(c1213[i][j][k]*l23 + c1113[i][j][k]*l13 + c1112[i][j][k]*l12);
//						g[2][2] = c1212[i][j][k]*l11 + c2222[i][j][k]*l22 + g22 + 2.0*(c2223[i][j][k]*l23 + c1223[i][j][k]*l13 + c1222[i][j][k]*l12);
//						g[3][3] = c1313[i][j][k]*l11 + c2323[i][j][k]*l22 + g33 + 2.0*(c2333[i][j][k]*l23 + c1333[i][j][k]*l13 + c1323[i][j][k]*l12);

						g[1] = c1111[i][j][k]*l11 + c1212[i][j][k]*l22 + g11 + 2.0*(c1213[i][j][k]*l23 + c1113[i][j][k]*l13 + c1112[i][j][k]*l12);
						g[2] = c1212[i][j][k]*l11 + c2222[i][j][k]*l22 + g22 + 2.0*(c2223[i][j][k]*l23 + c1223[i][j][k]*l13 + c1222[i][j][k]*l12);
						g[3] = c1313[i][j][k]*l11 + c2323[i][j][k]*l22 + g33 + 2.0*(c2333[i][j][k]*l23 + c1333[i][j][k]*l13 + c1323[i][j][k]*l12);

//						g[2][3] = c1213[i][j][k]*l11 + c2223[i][j][k]*l22 + g23 + c23*l23 + c33*l13 + c22*l12;
//						g[1][3] = c1113[i][j][k]*l11 + c1223[i][j][k]*l22 + g13 + c33*l23 + c13*l13 + c11*l12;
//						g[1][2] = c1112[i][j][k]*l11 + c1222[i][j][k]*l22 + g12 + c22*l23 + c11*l13 + c12*l12;

						g[4] = c1213[i][j][k]*l11 + c2223[i][j][k]*l22 + g23 + c23*l23 + c33*l13 + c22*l12;
						g[5] = c1113[i][j][k]*l11 + c1223[i][j][k]*l22 + g13 + c33*l23 + c13*l13 + c11*l12;
						g[6] = c1112[i][j][k]*l11 + c1222[i][j][k]*l22 + g12 + c22*l23 + c11*l13 + c12*l12;
  
//						g[3][2] = g[2][3];
//						g[3][1] = g[1][3];
//						g[2][1] = g[1][2];

						/* min./max. eigenvalues */
//						eigenvalue(g,e,3);
						eigenvalue(g,e);

						/* maximum model phase velocity at infinite frequency at this grid point */
						vmax1 = sqrt((e[1]/rho[i][j][k])*(1.0+L*taup[i][j][k]));

						/* minimum model phase velocity at center frequency of the source at this grid point */
						vmin1 = sqrt((e[3]/rho[i][j][k])*(1.0+sum*taus[i][j][k]));

						if (vmax1 > vmax) vmax = vmax1;
						if (vmin1 < vmin) vmin = vmin1;
					}
				}
				/* min. and max. grid spacing */
				dx = x[i+1]-x[i];
				dy = y[j+1]-y[j];
				dz = z[k+1]-z[k];

				dhmax = max(dx,dy);
				dhmax = max(dhmax,dz);
				dhmin = min(dx,dy);
				dhmin = min(dhmin,dz);

				vmax_loc  = max(vmax_loc,vmax);
				vmin_loc  = min(vmin_loc,vmin);
				dxmax_loc = max(dxmax_loc,dx);
				dxmin_loc = min(dxmin_loc,dx);
				dymax_loc = max(dymax_loc,dy);
				dymin_loc = min(dymin_loc,dy);
				dzmax_loc = max(dzmax_loc,dz);
				dzmin_loc = min(dzmin_loc,dz);

				dtstab     = dhmin/(vmax*sqrt(3.0));
				dtstab_loc = min(dtstab_loc,dtstab);
				if (vmin != 0.0){
					dhstab     = vmin/(10.0*fmax*dhmax);
					dhstab_loc = min(dhstab_loc,dhstab);
				}
			}
		}
	}
	
	/* deallocate memory for angle vectors */
	free_vector(sphi,1,nphi);
	free_vector(cphi,1,nphi);
	free_vector(stheta,1,ntheta);
	free_vector(ctheta,1,ntheta);


	dhmax_loc = max(dxmax_loc,dymax_loc);
	dhmax_loc = max(dhmax_loc,dzmax_loc);
	dhmin_loc = min(dxmin_loc,dymin_loc);
	dhmin_loc = min(dhmin_loc,dzmin_loc);
	
	/* deallocate memory */
//	free_matrix(g,1,3,1,3);
	free_vector(g,1,6);
	free_vector(e,1,3);


	fprintf(fp," Minimum and maximum velocities within subvolume (MYID=%d): \n",MYID);
	fprintf(fp," v_max(f=inf)=%e m/s\t v_min(f=fc)=%e m/s\n\n", vmax_loc,vmin_loc);

	fprintf(fp," Minimum and maximum grid spacing within subgrid (MYID=%d): \n",MYID);
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

	/* find out whether grid dispersion occurs and whether simulation is stable*/
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
		fprintf(fp," \t DT < DH_min/(v_max*sqrt(3))\n\n");
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