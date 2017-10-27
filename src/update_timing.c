/*------------------------------------------------------------------------
 *   recursive computation of statistics (average, standard deviation) 
 *   for timing information (wavefield update and exchange)
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_timing(int nt, double * time_update, double * time_sum, double * time_avg, double * time_std, int ntime){

	/* extern variables */
	extern int	OUT_TIMING, MYID;
	extern FILE	*FP; 
	extern char	LOG_FILE[STRING_SIZE];
  
	/* local variables */
	char	T_FILE[STRING_SIZE];
	FILE	*fp;
	int	ii;
	double	fak1, fak2;


	/* compute updates for average update and exchange times and their standard deviations */
	fak1 = ((double)(nt-1))/((double)nt);

	for (ii=1;ii<=ntime;ii++){
		fak2          = time_update[ii] - time_avg[ii];

		time_std[ii] += fak1*fak2*fak2;
		time_sum[ii] += time_update[ii];
		time_avg[ii]  = time_sum[ii]/((double)(nt));
	}

	/* output of timings to ASCII file (for later performance analysis) */
	if ((OUT_TIMING) && (!(MYID))){
		sprintf(T_FILE,"%s.timings",LOG_FILE);  
		/* fprintf(FP," PE 0 is writing timing information in ASCII format to %s \n", T_FILE); */
		if (nt > 1)
			fp=fopen(T_FILE,"a");
		else{
			fp=fopen(T_FILE,"w");
			if (ntime==7)
				fprintf(fp," %% timestep \t velocity update [s] \t PML velocity update [s] \t velocity exchange [s] \t stress update [s] \t PML stress update [s] \t stress exchange [s] \t total [s]\n");
			if (ntime==8)
				fprintf(fp," %% timestep \t velocity update [s] \t PML velocity update [s] \t velocity exchange [s] \t strain update [s] \t PML strain update [s] \t strain exchange [s] \t stress update [s] \t total [s]\n");
			if (ntime==9)
				fprintf(fp," %% timestep \t velocity update [s] \t PML velocity update [s] \t velocity exchange [s] \t strain update [s] \t PML strain update [s] \t strain exchange [s] \t stress update [s] \t stress exchange [s] \t total [s]\n");
		}
		fprintf(fp," %d \t\t ",nt);
		for (ii=1;ii<ntime;ii++)
			fprintf(fp," %e \t\t ",time_update[ii]); 
		fprintf(fp,"%e \n",time_update[ntime]);
		fclose(fp);
	}

}
