/*------------------------------------------------------------------------
 *   output timing information, i.e. real times for wavefield updates and exchange
 *   and some statitics (average, standard deviation)
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void timing(double * time_avg, double * time_std, int ishot){

	/* extern variables */
	extern int	NT;
	extern int	PML;
	extern int	NPROCX[3], MYID;
	extern FILE	*FP; 

	/* local variables */
	int	ii;


	/* compute 3 times standard deviation */
	for (ii=1;ii<=8;ii++){
		time_std[ii] = 3.0*sqrt(time_std[ii]/(NT-1));
	}

	/* log information */
	fprintf(FP,"\n **Message from function timing (printed by PE %d):\n",MYID);

	fprintf(FP," Average times and standard deviations (in percent) for \n");
	fprintf(FP," velocity update: \t %6.3f s +/- %6.3f s (%6.1f per cent of total timestep)\n",time_avg[1], time_std[1], 100.0*time_avg[1]/time_avg[8]);
	fprintf(FP," strain update:   \t %6.3f s +/- %6.3f s (%6.1f per cent of total timestep)\n",time_avg[4], time_std[4], 100.0*time_avg[4]/time_avg[8]);
	fprintf(FP," stress update:   \t %6.3f s +/- %6.3f s (%6.1f per cent of total timestep)\n",time_avg[6], time_std[6], 100.0*time_avg[6]/time_avg[8]);
	if (PML){
		fprintf(FP," PML velocity update: \t %6.3f s +/- %6.3f s (%6.1f per cent of total timestep)\n",time_avg[2], time_std[2], 100.0*time_avg[2]/time_avg[8]);
		fprintf(FP," PML stress update:   \t %6.3f s +/- %6.3f s (%6.1f per cent of total timestep)\n",time_avg[5], time_std[5], 100.0*time_avg[5]/time_avg[8]);
	}
	if (NPROCX[0]*NPROCX[1] > 1){
		fprintf(FP," velocity exchange: \t %6.3f s +/- %6.3f s (%6.1f per cent of total timestep)\n",time_avg[3], time_std[3], 100.0*time_avg[3]/time_avg[8]);
		fprintf(FP," stress exchange:   \t %6.3f s +/- %6.3f s (%6.1f per cent of total timestep)\n",time_avg[7], time_std[7], 100.0*time_avg[7]/time_avg[8]);
	}
	fprintf(FP," timestep:  \t\t %6.3f s +/- %6.3f s  \n",time_avg[8], time_std[8]);


	/* output of timings to ASCII file (for later performance analysis) */
//	if (OUT_TIMING){
//		sprintf(T_FILE,"%s.timings",LOG_FILE);  
//		fprintf(FP," PE 0 is writing timing information in ASCII format to %s \n", T_FILE);
//		fp = fopen(T_FILE,"w");
//		fprintf(fp," %% timestep \t  velocity update [s] \t strain update [s] \t stress update [s] \t PML velocity update [s] \t PML strain update [s] \t velocity exchange [s] \t stress exchange [s]\t total [s]\n");
//		for (nt=1;nt<=NT;nt++)
//			fprintf(fp," %d \t\t %e \t\t %e \t\t %e \t\t %e \t\t %e \t\t %e \t\t %e \t\t %e \n", 
//				nt, time_v_update[nt],time_e_update[nt],time_s_update[nt],time_v_pml_update[nt],time_e_pml_update[nt],time_v_exchange[nt],time_s_exchange[nt],time_timestep[nt]);
//		fclose(fp);
//	}
}
