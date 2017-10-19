/*------------------------------------------------------------------------
 *   Write program name etc to stdout
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," ******************************************************************************\n");
	fprintf(fp," This is program FD-3D.  \n");
	fprintf(fp," Parallel 3-D Viscoacoustic Finite Difference Modelling in      \n");
	fprintf(fp," Cartesian coordinates for seismic wave simulation    \n");
	fprintf(fp," (version with variable grid and PMLs)                           \n\n");
	fprintf(fp," Institute of Geophysics and Geoinformatics, TU Bergakademie Freiberg, Germany \n\n");
	fprintf(fp," ******************************************************************************\n");
	fprintf(fp,"\n");

}
