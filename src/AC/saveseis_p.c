/*------------------------------------------------------------------------
 *   write seismograms to files
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "fd3d.h"

void saveseis(FILE *fp, float **sectionvx, float **sectionvy, float **sectionvz, float **sectionp, 
float **sectionax, float **sectionay, float **sectionaz, float **sectiondiv, float **section_fulldata, 
int  **recpos, int  **recpos_loc, float ** srcpos, int nsrc, int ns, float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg, int *recswitch){

	/* extern variables */
	extern int	SEISMO, SEIS_FORMAT, MYID;
	extern int	FFID, NTR;
	extern char	SEIS_FILE[STRING_SIZE];

	/* local variables */
	char	seisfile[STRING_SIZE];
	char	ext[8];


	/* log information */
	fprintf(fp,"\n **Message from function saveseis (printed by PE %d):\n",MYID);
	
	switch(SEIS_FORMAT){
	case 1:
		sprintf(ext,".su");
		break;
	case 2:
		sprintf(ext,".asc");
		break;
	case 3:
		sprintf(ext,".bin");
		break;
	}

	fprintf(fp," Writing %d seismogram traces to\n",NTR);

	if (SEISMO & 8){
		/* particle acceleration */                
		catseis(sectionax, section_fulldata, recswitch, ns);
		if(!MYID){		   
		   sprintf(seisfile,"%s%.4i_ax%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xg,ypg,zpg);
		}

                catseis(sectionay, section_fulldata, recswitch, ns);
		if(!MYID){
		   sprintf(seisfile,"%s%.4i_ay%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,yg,zpg);
                }

                catseis(sectionaz, section_fulldata, recswitch, ns);
		if(!MYID){
		   sprintf(seisfile,"%s%.4i_az%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n\n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zg);
		}
	}

	if (SEISMO & 4){
		/* div */
	        catseis(sectiondiv, section_fulldata, recswitch, ns);	        
		if(!MYID){
		   sprintf(seisfile,"%s%.4i_div%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n\n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);
		}
	}
	if (SEISMO & 2){
		/* pressure */
                catseis(sectionp, section_fulldata, recswitch, ns);
		if(!MYID){
		   sprintf(seisfile,"%s%.4i_p%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n\n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);
		}
	}
	if (SEISMO & 1){
		/* particle velocity */ 
                catseis(sectionvx, section_fulldata, recswitch, ns);
		if(!MYID){
 		   sprintf(seisfile,"%s%.4i_vx%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xg,ypg,zpg);
		}

                catseis(sectionvy, section_fulldata, recswitch, ns);
		if(!MYID){
		   sprintf(seisfile,"%s%.4i_vy%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,yg,zpg);
		}

                catseis(sectionvz, section_fulldata, recswitch, ns);
		if(!MYID){
		   sprintf(seisfile,"%s%.4i_vz%s",SEIS_FILE,FFID,ext);
		   fprintf(fp,"\t %s \n\n",seisfile);
		   outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zg);
		}
	}
}
