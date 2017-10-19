/*-------------------------------------------------------------
 * Cat seismograms before output
 *  
 * Andre Kurzmann
 * Karlsruhe, 10.11.2009
 *-------------------------------------------------------------
*/

#include "fd.h"

void	catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns) {

	extern int      MYID;
	int		i, j, k;
	float		**fulldata2;

	/* temporary global data array for MPI-exchange */
	fulldata2 = matrix(1,ntr_glob,1,ns);

	k = 0;	/* trace counter for local data array */

	/* loop over global traces: copy traces of local array	*/
	/* to appropriate locations in the global array		*/
	for(i=1;i<=ntr_glob;i++)
	{
		
		if (recswitch[i]) {
			k++;
			for(j=1;j<=ns;j++)	fulldata2[i][j] = data[k][j];
		}
	}

	MPI_Allreduce(&fulldata2[1][1], &fulldata[1][1], ntr_glob*ns, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	free_matrix(fulldata2, 1,ntr_glob,1,ns);
}
