/*------------------------------------------------------------------------
 *   Write program name etc to stdout
 *
 *   D. Köhn, O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," ******************************************************************************\n");
	fprintf(fp," This is program FD-3D.  \n");
	fprintf(fp," Parallel 3-D Anisotropic (triclinic) Elastic Finite Difference Modelling in \n");
	fprintf(fp," Cartesian coordinates for seismic wave simulation \n");
	fprintf(fp," (version with variable grid and PMLs)             \n\n");
	fprintf(fp," Institute of Geophysics and Geoinformatics, TU Bergakademie Freiberg, Germany, and \n");
	fprintf(fp," Institute of Geophysics, Christian-Albrechts-University Kiel, Germany \n\n");
	fprintf(fp," ******************************************************************************\n");
	fprintf(fp,"\n");

}
