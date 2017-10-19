/*------------------------------------------------------------------------
 *   Read FD-Parameters from File or stdin
 *
 *   T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void read_par(FILE *fp_in){

	/* extern variables */
	extern int	NXG[3], SRCSIGNAL, MODEL, MODEL_FORMAT, SNAP, SNAP_FORMAT, L;
	extern float	TIME, DT, *FL, TAU, DAMPLE, DAMPRI, DAMPFR, DAMPBA, DAMPTO, DAMPBO;
	extern float	XREC1[3], XREC2[3];
	extern float	N_SIGMA_LE, SIGMA_LE, N_KAPPA_LE, KAPPA_LE, N_ALPHA_LE, ALPHA_LE;
	extern float	N_SIGMA_RI, SIGMA_RI, N_KAPPA_RI, KAPPA_RI, N_ALPHA_RI, ALPHA_RI;
	extern float	N_SIGMA_BA, SIGMA_BA, N_KAPPA_BA, KAPPA_BA, N_ALPHA_BA, ALPHA_BA;
	extern float	N_SIGMA_FR, SIGMA_FR, N_KAPPA_FR, KAPPA_FR, N_ALPHA_FR, ALPHA_FR;
	extern float	N_SIGMA_TO, SIGMA_TO, N_KAPPA_TO, KAPPA_TO, N_ALPHA_TO, ALPHA_TO;
	extern float	N_SIGMA_BO, SIGMA_BO, N_KAPPA_BO, KAPPA_BO, N_ALPHA_BO, ALPHA_BO;
	extern int	SEISMO, NDT, DXREC[3], SEIS_FORMAT, READMOD, READREC;
	extern int	LOG;
	extern int	ABTYPE_LE, ABTYPE_RI, ABTYPE_BA, ABTYPE_FR, ABTYPE_TO, ABTYPE_BO, FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern int	PERIODIC[3];
	extern int	OUT_TIMING;
	extern int	FFID;
	extern float	TSNAP1, TSNAP2, TSNAPINC, REFMOD[3], REFSRC[3], REFREC[3];
	extern float	XSNAPMIN, XSNAPMAX, YSNAPMIN, YSNAPMAX, ZSNAPMIN, ZSNAPMAX;
	extern char	DX_FILE[STRING_SIZE], DY_FILE[STRING_SIZE], DZ_FILE[STRING_SIZE];
	extern char	MFILE[STRING_SIZE], MODEL_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
	extern char	SNAP_FILE[STRING_SIZE], SRC_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char	SEIS_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern int	NPROCX[3], MODEL_IDX[3], IDX[3], CHECKPTREAD, CHECKPTWRITE; 

	/* local variables */
	char	s[74];
	int	c=0, lineno=0, l;


	while ((c=getc(fp_in)) != EOF){
		if ((c=='\n') && (getc(fp_in)!='#')){     
			lineno++;
			/* printf(" reading line %d \n",lineno);*/
			switch (lineno){
			case 1 :
				fscanf(fp_in,"%s =%i",s,&NPROCX[0]);
				break;
			case 2 :
				fscanf(fp_in,"%s =%i",s,&NPROCX[1]);
				break;
			case 3 :
				fscanf(fp_in,"%s =%i",s,&NPROCX[2]);
				break;
			case 4 :
				fscanf(fp_in,"%s =%i",s,&NXG[0]);
				break;
			case 5 :
				fscanf(fp_in,"%s =%i",s,&NXG[1]);
				break;
			case 6 :
				fscanf(fp_in,"%s =%i",s,&NXG[2]);
				break;
			case 7 :
				fscanf(fp_in,"%s =%s",s,DX_FILE);  
				break;
			case 8 :
				fscanf(fp_in,"%s =%s",s,DY_FILE);  
				break;
			case 9 :
				fscanf(fp_in,"%s =%s",s,DZ_FILE);
				break;
			case 10 :  
				fscanf(fp_in,"%s =%f ,%f ,%f",s,&REFMOD[0],&REFMOD[1],&REFMOD[2]);
				break;
			case 11 :
				fscanf(fp_in,"%s =%f",s,&TIME);
				break;
			case 12 :
				fscanf(fp_in,"%s =%f",s,&DT);
				break;
			case 13 :
				fscanf(fp_in,"%s =%f ,%f ,%f",s,&REFSRC[0],&REFSRC[1],&REFSRC[2]);
				break;
			case 14 :
				fscanf(fp_in,"%s =%i",s,&SRCSIGNAL);
				break;
			case 15 :
				fscanf(fp_in,"%s =%s",s,SIGNAL_FILE);  
				break;
			case 16 :
				fscanf(fp_in,"%s =%s",s,SRC_FILE);
				break;
			case 17 :
				fscanf(fp_in,"%s =%i",s,&READMOD);
				break;
			case 18 :
				fscanf(fp_in,"%s =%s",s,MFILE);
				break;
			case 19 :
				fscanf(fp_in,"%s =%i",s,&MODEL);
				break;
			case 20 :
				fscanf(fp_in,"%s =%i",s,&MODEL_IDX[0]);
				break;
			case 21 :
				fscanf(fp_in,"%s =%i",s,&MODEL_IDX[1]);
				break;
			case 22 :
				fscanf(fp_in,"%s =%i",s,&MODEL_IDX[2]);
				break;
			case 23 :
				fscanf(fp_in,"%s =%i",s,&MODEL_FORMAT);
				break;
			case 24 :
				fscanf(fp_in,"%s =%s",s,MODEL_FILE);
				break;
			case 25 :
				fscanf(fp_in,"%s =%i",s,&L);
				break;
			case 26 :
				FL=vector(1,L);
				fscanf(fp_in,"%s =%f",s,&FL[1]);
				for (l=2;l<=L;l++)
					fscanf(fp_in,"%f",&FL[l]);
				break;
			case 27 :
				fscanf(fp_in,"%s =%f",s,&TAU);
				TAU=2.0/TAU; /* Q-factor is input in parameter file */
				break;
			case 28 :
				fscanf(fp_in,"%s =%i",s,&PERIODIC[0]);
				if (PERIODIC[0]<0)
					PERIODIC[0]=0;
				break;
			case 29 :
				fscanf(fp_in,"%s =%i",s,&ABTYPE_LE);
				break;
			case 30 :
				fscanf(fp_in,"%s =%i",s,&FWLE);
				if (FWLE<0) 
					FWLE=0;
				break;
			case 31 :
				fscanf(fp_in,"%s =%f",s,&DAMPLE);
				break;
			case 32 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_SIGMA_LE,&SIGMA_LE);
				break;
			case 33 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_KAPPA_LE,&KAPPA_LE);
				break;
			case 34 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_ALPHA_LE,&ALPHA_LE);
				break;
			case 35 :
				fscanf(fp_in,"%s =%i",s,&ABTYPE_RI);
				break;
			case 36 :
				fscanf(fp_in,"%s =%i",s,&FWRI);
				if (FWRI<0) 
					FWRI=0;
				break;
			case 37 :
				fscanf(fp_in,"%s =%f",s,&DAMPRI);
				break;
			case 38 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_SIGMA_RI,&SIGMA_RI);
				break;
			case 39 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_KAPPA_RI,&KAPPA_RI);
				break;
			case 40 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_ALPHA_RI,&ALPHA_RI);
				break;
			case 41 :
				fscanf(fp_in,"%s =%i",s,&PERIODIC[1]);
				if (PERIODIC[1]<0)
					PERIODIC[1]=0;
				break;
			case 42 :
				fscanf(fp_in,"%s =%i",s,&ABTYPE_BA);
				break;
			case 43 :
				fscanf(fp_in,"%s =%i",s,&FWBA);
				if (FWBA<0) 
					FWBA=0;
				break;
			case 44 :
				fscanf(fp_in,"%s =%f",s,&DAMPBA);
				break;
			case 45 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_SIGMA_BA,&SIGMA_BA);
				break;
			case 46 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_KAPPA_BA,&KAPPA_BA);
				break;
			case 47 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_ALPHA_BA,&ALPHA_BA);
				break;
			case 48 :
				fscanf(fp_in,"%s =%i",s,&ABTYPE_FR);
				break;
			case 49 :
				fscanf(fp_in,"%s =%i",s,&FWFR);
				if (FWFR<0) 
					FWFR=0;
				break;
			case 50 :
				fscanf(fp_in,"%s =%f",s,&DAMPFR);
				break;
			case 51 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_SIGMA_FR,&SIGMA_FR);
				break;
			case 52 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_KAPPA_FR,&KAPPA_FR);
				break;
			case 53 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_ALPHA_FR,&ALPHA_FR);
				break;
			case 54 :
				fscanf(fp_in,"%s =%i",s,&PERIODIC[2]);
				if (PERIODIC[2]<0)
					PERIODIC[2]=0;
				break;
			case 55 :
				fscanf(fp_in,"%s =%i",s,&ABTYPE_TO);
				break;
			case 56 :
				fscanf(fp_in,"%s =%i",s,&FWTO);
				if (FWTO<0) 
					FWTO=0;
				break;
			case 57 :
				fscanf(fp_in,"%s =%f",s,&DAMPTO);
				break;
			case 58 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_SIGMA_TO,&SIGMA_TO);
				break;
			case 59 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_KAPPA_TO,&KAPPA_TO);
				break;
			case 60 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_ALPHA_TO,&ALPHA_TO);
				break;
			case 61 :
				fscanf(fp_in,"%s =%i",s,&ABTYPE_BO);
				break;
			case 62 :
				fscanf(fp_in,"%s =%i",s,&FWBO);
				if (FWBO<0) 
					FWBO=0;
				break;
			case 63 :
				fscanf(fp_in,"%s =%f",s,&DAMPBO);
				break;
			case 64 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_SIGMA_BO,&SIGMA_BO);
				break;
			case 65 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_KAPPA_BO,&KAPPA_BO);
				break;
			case 66 :
				fscanf(fp_in,"%s =%f ,%f",s,&N_ALPHA_BO,&ALPHA_BO);
				break;
			case 67 :
				fscanf(fp_in,"%s =%i",s,&SNAP);
				break;
			case 68 :
				fscanf(fp_in,"%s =%f",s,&TSNAP1);
				break;
			case 69 :
				fscanf(fp_in,"%s =%f",s,&TSNAP2);
				break;
			case 70 :
				fscanf(fp_in,"%s =%f",s,&TSNAPINC);
				break;
			case 71 :
				fscanf(fp_in,"%s =%f",s,&XSNAPMIN);
				break;
			case 72 :
				fscanf(fp_in,"%s =%f",s,&XSNAPMAX);
				break;
			case 73 :
				fscanf(fp_in,"%s =%i",s,&IDX[0]);
				break;
			case 74 :
				fscanf(fp_in,"%s =%f",s,&YSNAPMIN);
				break;
			case 75 :
				fscanf(fp_in,"%s =%f",s,&YSNAPMAX);
				break;
			case 76 :
				fscanf(fp_in,"%s =%i",s,&IDX[1]);
				break;
			case 77 :
				fscanf(fp_in,"%s =%f",s,&ZSNAPMIN);
				break;
			case 78 :
				fscanf(fp_in,"%s =%f",s,&ZSNAPMAX);
				break;   
			case 79 :
				fscanf(fp_in,"%s =%i",s,&IDX[2]);
				break;
			case 80 :
				fscanf(fp_in,"%s =%i",s,&SNAP_FORMAT);
				break;
			case 81 :
				fscanf(fp_in,"%s =%s",s,SNAP_FILE);
				break;
			case 82 :
				fscanf(fp_in,"%s =%i",s,&READREC);
				break;
			case 83 :
				fscanf(fp_in,"%s =%s",s,REC_FILE);
				break;
			case 84 :
				fscanf(fp_in,"%s =%f ,%f ,%f",s,&REFREC[0],&REFREC[1],&REFREC[2]);
				break;
			case 85 :
				fscanf(fp_in,"%s =%f ,%f ,%f",s,&XREC1[0],&XREC1[1],&XREC1[2]);
				break;
			case 86 :
				fscanf(fp_in,"%s =%f ,%f ,%f",s,&XREC2[0],&XREC2[1],&XREC2[2]);
				break;
			case 87 :
				fscanf(fp_in,"%s =%i ,%i ,%i",s,&DXREC[0],&DXREC[1],&DXREC[2]);
				break;
			case 88 :
				fscanf(fp_in,"%s =%i",s,&SEISMO);
				break;
			case 89 :
				fscanf(fp_in,"%s =%i",s,&NDT);
				break;
			case 90 :
				fscanf(fp_in,"%s =%i",s,&SEIS_FORMAT);
				break;
			case 91 :
				fscanf(fp_in,"%s =%s",s,SEIS_FILE);
				break;
			case 92 :
				fscanf(fp_in,"%s =%i",s,&FFID);
				break;
			case 93 :
				fscanf(fp_in,"%s =%s",s,LOG_FILE);
				break;
			case 94 :
				fscanf(fp_in,"%s =%i",s,&LOG);
				break;
			case 95 :
				fscanf(fp_in,"%s =%i",s,&OUT_TIMING);
				break;
			case 96 :
				fscanf(fp_in,"%s =%i",s,&CHECKPTREAD);
				break;
			case 97 :
				fscanf(fp_in,"%s =%i",s,&CHECKPTWRITE);
				break;
			case 98 :
				fscanf(fp_in,"%s =%s",s,CHECKPTFILE);
				break;
			default:
				break;
			}
		}

	}
	fclose(fp_in);

}
