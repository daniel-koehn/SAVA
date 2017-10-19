/*------------------------------------------------------------------------
 *   merge model files written by the different processes to 
 *   a single file                                 
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void mergemod(char modfile[STRING_SIZE], int format){

	/* extern variables */
	extern int	MYID, NPROCX[3];
	extern int	NXG[3], NX[3], NPROC, MODEL_IDX[3];
	extern FILE	*FP;

	/* local variables */
	int	i, j, k, ip, jp, kp;
	float	a;
	char	file[STRING_SIZE];
	char    fileoutxy[STRING_SIZE], fileoutyz[STRING_SIZE], fileoutxz[STRING_SIZE];
	FILE	*fp[NPROCX[0]][NPROCX[1]][NPROCX[2]], *fpout;
	FILE	*fpoutyz[NPROCX[0]], *fpoutxz[NPROCX[1]], *fpoutxy[NPROCX[2]];



	fprintf(FP,"\n **Message from mergemod (printed by PE %d):\n",MYID);
	fprintf(FP," Starting to merge %d model files \n",NPROC);
	fprintf(FP," Writing merged model file to \n\t %s \n",modfile);

	fpout   = fopen(modfile,"w");
	for (ip=0;ip<NPROCX[0];ip++){
		sprintf(fileoutyz,"%s_yz%i",modfile,ip);
		fpoutyz[ip] = fopen(fileoutyz,"w");
	}
	for (jp=0;jp<NPROCX[1];jp++){
		sprintf(fileoutxz,"%s_xz%i",modfile,jp);
		fpoutxz[jp] = fopen(fileoutxz,"w");
	}
	for (kp=0;kp<NPROCX[2];kp++){
		sprintf(fileoutxy,"%s_xy%i",modfile,kp);
		fpoutxy[kp] = fopen(fileoutxy,"w");
	}

	fprintf(FP," Opening model files: %s.??? ",modfile);
	for (ip=0;ip<NPROCX[0]; ip++){
		for (jp=0;jp<NPROCX[1]; jp++){
			for (kp=0;kp<NPROCX[2]; kp++){
				sprintf(file,"%s.%i.%i.%i",modfile,ip,jp,kp);
				fp[ip][jp][kp] = fopen(file,"r");
				if (fp[ip][jp][kp]==NULL) 
					error("merge: can't read model file !");
			}
		}
	}

	fprintf(FP," ... finished. \n");



	fprintf(FP," Copying ...");

	/* loop over processes */
	for (ip=0;ip<NPROCX[0];ip++){
		/* loop over gridpoints */
		for (i=1;i<=NX[0];i+=MODEL_IDX[0]){
			/* loop over processes */
			for (jp=0;jp<NPROCX[1];jp++){
				/* loop over gridpoints */
				for (j=1;j<=NX[1];j+=MODEL_IDX[1]){
					/* loop over processes */
					for (kp=0;kp<NPROCX[2];kp++){
						/* loop over gridpoints */
						for (k=1;k<=NX[2];k+=MODEL_IDX[2]){
							a = readdsk(fp[ip][jp][kp],format);
							if (i == (NX[0]/(2*MODEL_IDX[0]))*MODEL_IDX[0]+1)
								writedsk(fpoutyz[ip],a,format);
							if (j == (NX[1]/(2*MODEL_IDX[1]))*MODEL_IDX[1]+1)
								writedsk(fpoutxz[jp],a,format);
							if (k == (NX[2]/(2*MODEL_IDX[2]))*MODEL_IDX[2]+1)
								writedsk(fpoutxy[kp],a,format);
							writedsk(fpout,a,format);
						}
					}
				}
			}
		}
	}
	fprintf(FP," finished. \n");

	for (ip=0;ip<NPROCX[0]; ip++){
		for (jp=0;jp<NPROCX[1]; jp++){
			for (kp=0;kp<NPROCX[2]; kp++){
				fclose(fp[ip][jp][kp]);
			}
		}
	}
	fclose(fpout);
	for (ip=0;ip<NPROCX[0];ip++)
		fclose(fpoutyz[ip]);
	for (jp=0;jp<NPROCX[1];jp++)
		fclose(fpoutxz[jp]);
	for (kp=0;kp<NPROCX[2];kp++)
		fclose(fpoutxy[kp]);

	for (ip=0;ip<NPROCX[0];ip++){
		fprintf(FP," You may use the SU program \n\n");
		fprintf(FP," ximage n1=%d < %s_yz%i  label1=Z label2=Y title=%s_yz%i \n\n",((NXG[2]-1)/MODEL_IDX[2])+1,modfile,ip,modfile,ip);
		fprintf(FP," to visualize/check the model in the y-z-plane at ix=%i .\n",NX[0]*ip+(NX[0]/(2*MODEL_IDX[0]))*MODEL_IDX[0]+1);
	}
	for (jp=0;jp<NPROCX[1];jp++){
		fprintf(FP," You may use the SU program \n\n");
		fprintf(FP," ximage n1=%d < %s_xz%i  label1=Z label2=X title=%s_xz%i \n\n",((NXG[2]-1)/MODEL_IDX[2])+1,modfile,jp,modfile,jp);
		fprintf(FP," to visualize/check the model in the x-z-plane at iy=%i .\n",NX[1]*jp+(NX[1]/(2*MODEL_IDX[1]))*MODEL_IDX[1]+1);
	}
	for (kp=0;kp<NPROCX[2];kp++){
		fprintf(FP," You may use the SU program \n\n");
		fprintf(FP," ximage n1=%d < %s_xy%i  label1=Y label2=X title=%s_xy%i \n\n",((NXG[1]-1)/MODEL_IDX[1])+1,modfile,kp,modfile,kp);
		fprintf(FP," to visualize/check the model in the x-y-plane at iz=%i .\n",NX[2]*kp+(NX[2]/(2*MODEL_IDX[2]))*MODEL_IDX[2]+1);
	}

}


