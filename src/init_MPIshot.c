/*
 * Initiate MPI parallelization by distributing shots over MPI processes
 *
 * Daniel Koehn
 * Kiel, 11.12.2015
 */

#include "fd.h"

void init_MPIshot(){

        /* global variables */
	extern int NPROCSHOT, COLOR, NSHOT1, NSHOT2, NSRC;
	/* extern int MYID, NP, MYID_SHOT; */

	/* local variables */
	int i, j;

	/* count number of shot positions */
	countsrc();	

	/* distribute shots over MPI processes */
	NSHOT1 = (NSRC / NPROCSHOT) * COLOR;

	if (NSRC % NPROCSHOT > COLOR){
	  NSHOT1 += COLOR;
	  NSHOT2 = NSHOT1 + (NSRC / NPROCSHOT) + 1;
	}else{
	  NSHOT1 += NSRC % NPROCSHOT;
	  NSHOT2 = NSHOT1 + (NSRC / NPROCSHOT);
	}

	NSHOT1++;
	NSHOT2++;

	/*printf("WORLD RANK/SIZE: %d/%d \t shot_comm COLOR/RANK/SIZE/NSHOT1/NSHOT2: %d/%d/%d/%d/%d\n", MYID, NP, COLOR, MYID_SHOT, NPROCSHOT, NSHOT1, NSHOT2);*/

}



