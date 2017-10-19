/*------------------------------------------------------------------------
*   Calculating source signal at different source positions with different
*   time-shift, centre frequency and amplitude (as specified in SOURCE_FILE).
*   Source signals are written to array signals 
*
*   O. Hellwig, S. Wenk
*  ----------------------------------------------------------------------*/

#include "fd.h"
#include "fd3d.h"

float ** wavelet(float ** srcpos_loc, int nsrc_loc){

	/* extern variables */
	extern int	SRCSIGNAL, NT, MYID;
	extern float	DT;
	extern FILE	*FP;

	/*local variables */
	int	k, nt, nsamp=0, ns;
	float	tshift, amp=0.0, a, fc, tau, t, ts;
	float	**source=NULL, **signals;


	if (SRCSIGNAL>=4) 
		source = read_wavelet(srcpos_loc, nsrc_loc,&nsamp);

	signals = matrix(1,nsrc_loc,1,NT);

	for (k=1;k<=nsrc_loc;k++) {
		tshift = srcpos_loc[k][4];
		fc     = srcpos_loc[k][5];
		a      = srcpos_loc[k][6];
		ts     = 2.0/fc;
		ns     = 1;

		for (nt=1;nt<=NT;nt++){
			t = (float)nt*DT;

			switch (SRCSIGNAL){
				case 1 : /* Ricker signal shifted by 2/fc */
					tau = PI*((t-tshift)*fc-2.0)/sqrt(2.0);
					amp = ((1.0-4.0*tau*tau)*exp(-2.0*tau*tau));
					break;
				case 2 : 
					if ((t<tshift) || (t>(tshift+ts))) 
						amp=0.0;
					else 
						amp = ((sin(2.0*PI*(t-tshift)*fc) - 0.5*sin(4.0*PI*(t-tshift)*fc)));
					break;
				case 3 : /* sinus raised to the power of three */
					if ((t<tshift) || (t>(tshift+ts))) 
						amp = 0.0;
					else 
						amp = pow(sin(PI*(t-tshift)*fc),3.0);
					break; 
				case 4 : /* source wavelet from ASCII-file */
					if ((t<tshift) || (t>srcpos_loc[k][8]*DT+tshift)) 
						amp=0.0;
					else{
						amp = source[k][ns];
						ns++;
					}
					break;  
				case 5 :  /* source wavelet from BIN-file */
					if ((t<tshift) || (t>srcpos_loc[k][8]*DT+tshift)) 
						amp=0.0;
					else{
						amp = source[k][ns];
						ns++;
					}
					break; 
				default : 
					error("Which source-wavelet ? ");
			}
			signals[k][nt] = amp*a;
		}
	}
	fprintf(FP," Message from function wavelet written by PE %d \n",MYID);
	fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc_loc,MYID);
	fprintf(FP," have been assigned with a source signal. \n");
	
	if (SRCSIGNAL>=4)
		free_matrix(source,1,nsrc_loc,1,nsamp);


	return signals;

}
