/*------------------------------------------------------------------------
 *  Parallel 3-D Viscoelastic Finite Difference Modelling 
 *  in Cartesian coordinates,
 *  merging snapshot files
 * 
 * In case of questions or if you require further information
 * you are welcome to contact the authors:                       
 *
 *  ----------------------------------------------
 *  Olaf Hellwg
 *  TU Bergakademie Freiberg
 *  Institut fuer Geophysik
 *  Gustav-Zeuner-Str. 12
 *  D-09596 Freiberg, Germany
 *  Phone: +49 (0)3731 392233
 *  Fax: +49 (0)3731 392636
 *  e-mail: olaf.hellwig@geophysik.tu-freiberg.de
 *  http://www.geophysik.tu-freiberg.de
 *
 *  Prof. Dr. Thomas Bohlen
 *  Karlsruher Institut fuer Technologie (KIT)
 *  Geophysikalisches Institut (GPI)
 *  Hertzstrasse 16, Geb. 6.42
 *  D-76187 Karlsruhe, Germany
 *  Phone: +49 (0)721 608 4416
 *  Fax: +49 (0)721 71173
 *  e-mail:thomas.bohlen@kit.edu
 *  http://www.gpi.kit.edu
 *  ----------------------------------------------
 */

#include "fd.h"		/* general include file for all FD programs */
#include "fd3d.h"	/* general include file for FD3D program */

#include "globvar.h"	/* definition of global variables  */

int main(int argc, char **argv){

/* local variables */
float	*xpg, *ypg, *zpg;
int	nsnap;


/* read parameters from parameter-file (stdin) */
read_par(fopen(argv[1],"r")); 

/* size of local grid */
NX[0] = NXG[0]/NPROCX[0];
NX[1] = NXG[1]/NPROCX[1];
NX[2] = NXG[2]/NPROCX[1];

FP=stdout;

/*allocate memory for x-, y- and z-grid point positions*/
xpg = vector(0,NXG[0]+1);
ypg = vector(0,NXG[1]+1);
zpg = vector(0,NXG[2]+1);

/*read files with grid-spacing*/
read_grid_glob(FP,DX_FILE,0,xpg);
read_grid_glob(FP,DY_FILE,1,ypg);
read_grid_glob(FP,DZ_FILE,2,zpg);

nsnap=1+iround((TSNAP2-TSNAP1)/TSNAPINC);


if (SNAP & 16){
	/* stress components */
	merge(nsnap,1,xpg,ypg,zpg);
	merge(nsnap,2,xpg,ypg,zpg);
	merge(nsnap,3,xpg,ypg,zpg);
	merge(nsnap,4,xpg,ypg,zpg);
	merge(nsnap,5,xpg,ypg,zpg);
	merge(nsnap,6,xpg,ypg,zpg);
}
if (SNAP & 8){
	/* particle acceleration */
	merge(nsnap,7,xpg,ypg,zpg);
	merge(nsnap,8,xpg,ypg,zpg);
        merge(nsnap,9,xpg,ypg,zpg);
}
if (SNAP & 4){
	/* div and curl */
	merge(nsnap,10,xpg,ypg,zpg);
	merge(nsnap,11,xpg,ypg,zpg);
	merge(nsnap,12,xpg,ypg,zpg);
	merge(nsnap,13,xpg,ypg,zpg);
}
if (SNAP & 2){
	/* pressure */
	merge(nsnap,14,xpg,ypg,zpg);
}
if (SNAP & 1){
	/* particle velocity */
	merge(nsnap,15,xpg,ypg,zpg);
	merge(nsnap,16,xpg,ypg,zpg);
	merge(nsnap,17,xpg,ypg,zpg);
}

/*deallocate memory for x-, y- and z-grid point positions*/
free_vector(xpg,0,NXG[0]+1);
free_vector(ypg,0,NXG[1]+1);
free_vector(zpg,0,NXG[2]+1);

return 0;

}
