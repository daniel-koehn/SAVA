/*
 *  SAVA_TRICL: 3D parallel triclinic elastic time domain FD modelling, 
 *  FWI and RTM code  
 *
 *  Daniel Koehn
 *  Kiel, 14/08/2017
 *
 *  The original forward modelling code fd3dtricl was developed by Olaf Hellwig, 
 *  Denise De Nil and Daniel Köhn:
 *  ------------------------------------------------------------------------
 *  Dr. Olaf Hellwig
 *  TU Bergakademie Freiberg
 *  Institut fuer Geophysik
 *  Gustav-Zeuner-Str. 12
 *  D-09596 Freiberg, Germany
 *  Phone: +49 (0)3731 392233
 *  Fax: +49 (0)3731 392636
 *  e-mail: olaf.hellwig[at]geophysik.tu-freiberg.de
 *  http://www.geophysik.tu-freiberg.de
 *  ------------------------------------------------------------------------
 *
 *  The FWI and RTM code was developed by Denise De Nil and Daniel Köhn
 *  ------------------------------------------------------------------------
 *  Dr. Daniel Koehn 
 *  Christian-Albrechts Universität Kiel
 *  Institute of Geoscience, Department of Geophysics
 *  Otto-Hahn-Platz 1
 *  D-24098 Kiel, Germany 
 *  Phone: +49 431 880 4566
 *  e-mail:dkoehn[at]geophysik.uni-kiel.de,
 *  http://www.geophysik.uni-kiel.de/~dkoehn
 *  ------------------------------------------------------------------------
 *
 *  SAVA is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, version 2.0 of the License only. 
 *  
 *  SAVA is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 *  
 *  You should have received a copy of the GNU General Public License 
 *  along with SAVA (see file LICENSE.md) 
 *
 *  If you show modelling/inversion results in a paper or presentation please 
 *  give a reference to the following papers:
 *
 *  Zehner, B., Hellwig, O., Linke, M., Görz, I., Buske, S. (2016) 
 *  "Rasterizing geological models for parallel finite difference simulation 
 *  using seismic simulation as an example." 
 *  Computers & Geosciences 86, 83-91. 
 *
 *  Köhn, Daniel, Olaf Hellwig, Denise De Nil, and Wolfgang Rabbel (2015) 
 *  "Waveform inversion in triclinic anisotropic media -- A resolution study." 
 *  Geophysical Journal International, 201(3), 1642-1656. 
 * 
 *  
 *  Thank you for your co-operation, 
 *  Daniel Koehn
 */

#include "fd.h"		/* general include file for all FD programs */
#include "fd3d.h"	/* general include file for FD3D program */
#include "fd3daniso.h"	/* include file for anisotropic version */
#include "fd3dtricl.h"	/* include file for triclinic version */
#include "fd3del.h"	/* include file for elastic version */

#include "globvar.h"	/* definition of global variables  */

int main(int argc, char **argv){

	/* global variables */
	// extern int MODE;

	/* Initialize MPI environment */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&NP);
	MPI_Comm_rank(MPI_COMM_WORLD,&MYID);
        
	setvbuf(stdout, NULL, _IONBF, 0);

	if (!(MYID)){

	    /* print program name, version etc to stdout */
	    info(stdout);

	    /* PE 0 is reading the parameters from the input file */
	    read_par(fopen(argv[1],"r"));

	    /* PE 0 is checking the parameters from input file */
	    check_par(); 
	}

        /* PE 0 will broadcast the parameters to all others PEs */
        exchange_par(); 

	/* 3D triclinic elastic forward problem */
	//if(MODE==0){
	   FD_TRICL();
	//}

	/* 3D triclinic elastic Full Waveform Inversion */
	/*if(MODE==1){
	   FWI_TRICL();
	}*/

        /* 3D triclinic elastic Reverse Time Migration */
	/*if(MODE==2){
	   RTM_TRICL();
	}*/

	MPI_Finalize();
	return 0;

}
