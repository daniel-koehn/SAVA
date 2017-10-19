/*------------------------------------------------------------------------
 *   Read extern source wavelet                               
 *
 *   O. Hellwig, S. Wenk
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float **read_wavelet(float **srcpos_loc, int nsrc_loc, int *nsamp_max){

	/* extern variables */
	extern int	SRCSIGNAL;
	extern int	NT;
	extern char	SIGNAL_FILE[STRING_SIZE];
	
	/* local variables */
	float	**source=NULL;
	int	i, k, c, nsamp;
	FILE	*fp_source;


	nsamp = 0;

	switch (SRCSIGNAL){
		case 4:
			fp_source = fopen(SIGNAL_FILE,"r");
			if (fp_source==NULL) 
				error(" Source file (ASCII) could not be opened !");

			if (srcpos_loc[1][8]==0) {
				/* number of samples */
				while ((c = fgetc(fp_source)) != EOF)
					if (c=='\n')
						nsamp++;
				rewind(fp_source);
				printf(" Number of samples (nsamp) in source file: %i\n",nsamp);

				if (nsamp > NT) 
					nsamp=NT;
				source = matrix(1,nsrc_loc,1,nsamp);

				/* read one source signal from ascii file */
				for (i=1;i<=nsrc_loc;i++){
					for (k=1;k<=nsamp;k++) 
						fscanf(fp_source,"%e\n",&source[i][k]);
					rewind(fp_source);
					srcpos_loc[i][8] = (float)(nsamp);
				}
			}
			else{
				/* find maximum number of samples */
				for (i=1;i<=nsrc_loc;i++){
					if(nsamp  >= (int)(srcpos_loc[i][8])) 
						continue;
					else 
						nsamp = (int)(srcpos_loc[i][8]);
				}
				printf(" Max. number of samples (nsamp) per trace in source file: %i\n",nsamp);

				source=matrix(1,nsrc_loc,1,nsamp);

				/* read each source signal from ascii file */
				for (i=1;i<=nsrc_loc;i++){
					for (k=1;k<=(int)(srcpos_loc[i][9]);k++)
						fscanf(fp_source,"%e\n",&source[i][1]);
					for (k=1;k<=(int)(srcpos_loc[i][8]);k++)
						fscanf(fp_source,"%e\n",&source[i][k]);
					for (k=(int)(srcpos_loc[i][8])+1;k<=nsamp;k++)
						source[i][k] = 0.0;
					rewind(fp_source);
				}
			}
			fclose(fp_source);
			break;
		case 5:
			fp_source = fopen(SIGNAL_FILE,"rb");
			if (fp_source==NULL) 
				error(" Source file (BIN) could not be opened !");

			if (srcpos_loc[1][8]==0){
				/* number of samples */
				fseek(fp_source,0,SEEK_END);
				nsamp = (int)(ftell(fp_source)/sizeof(float));
				rewind(fp_source);
				printf(" Number of samples (nsamp) in source file: %i\n",nsamp);

				if (nsamp > NT) 
					nsamp=NT;
				source = matrix(1,nsrc_loc,1,nsamp);

				/* read one source signal from bin file */
				for (i=1;i<=nsrc_loc;i++){
					for (k=1;k<=nsamp;k++) 
						fread(&source[i][k],sizeof(float),1,fp_source);
					rewind(fp_source);
					srcpos_loc[i][8] = (float)(nsamp);
				}
			}
			else{
				/* find maximum number of samples */
				for (i=1;i<=nsrc_loc;i++){
					if(nsamp  >= (int)(srcpos_loc[i][8])) 
						continue;
					else 
						nsamp = (int)(srcpos_loc[i][8]);
				}
				printf(" Max. number of samples (nsamp) per trace in source file: %i\n",nsamp);

				source=matrix(1,nsrc_loc,1,nsamp);

				/* read each source signal from bin file */
				for (i=1;i<=nsrc_loc;i++){
					fseek(fp_source,(int)(srcpos_loc[i][9])*sizeof(float),SEEK_SET);
					for (k=1;k<=(int)(srcpos_loc[i][8]);k++) 
						fread(&source[i][k],sizeof(float),1,fp_source);
					for (k=(int)(srcpos_loc[i][8])+1;k<=nsamp;k++) 
						source[i][k] = 0.0;
				}

			}
			fclose(fp_source);
			break;
		default : 
			error("Which source-wavelet ? ");
	}
	*nsamp_max = nsamp;
	return source;
}
