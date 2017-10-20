/*------------------------------------------------------------------------
 *   write seismograms to files
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "fd3d.h"

void saveseis(FILE *fp, float **sectionvx, float **sectionvy, float **sectionvz, float **sectionax, float **sectionay, float **sectionaz, 
float **sectiondiv, float **sectioncurlx, float **sectioncurly, float **sectioncurlz, float **sectionp,
float **sectiontxx, float **sectiontxy, float **sectiontxz, float **sectiontyy, float **sectiontyz, float **sectiontzz, float **section_fulldata,
int  **recpos, int  **recpos_loc, float ** srcpos, int nsrc, int ns, float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg, 
int *recswitch){

	/* extern variables */
	extern int	SEISMO, SEIS_FORMAT, MYID, NTR;
	extern int	FFID;
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

	if (SEISMO & 16){
		/* stress */
                catseis(sectiontxx, section_fulldata, recswitch, ns);
		sprintf(seisfile,"%s%.4i_txx%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectiontxy, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_txy%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectiontxz, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_txz%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectiontyy, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_tyy%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectiontyz, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_tyz%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectiontzz, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_tzz%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n\n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);
	}
	if (SEISMO & 8){
		/* particle acceleration */
		catseis(sectionax, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_ax%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectionay, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_ay%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectionaz, section_fulldata, recswitch, ns);
		sprintf(seisfile,"%s%.4i_az%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n\n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);
	}
	if (SEISMO & 4){
		/* div and curl */
		catseis(sectiondiv, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_div%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectioncurlx, section_fulldata, recswitch, ns);
		sprintf(seisfile,"%s%.4i_curlx%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectioncurly, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_curly%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectioncurlz, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_curlz%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n\n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);
	}
	if (SEISMO & 2){
		/* pressure */
		catseis(sectionp, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_p%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n\n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);
	}
	if (SEISMO & 1){
		/* particle velocity */ 
		catseis(sectionvx, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_vx%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectionvy, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_vy%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);

		catseis(sectionvz, section_fulldata, recswitch,ns);
		sprintf(seisfile,"%s%.4i_vz%s",SEIS_FILE,FFID,ext);
		fprintf(fp,"\t %s \n\n",seisfile);
		outseis(fp,fopen(seisfile,"w"),section_fulldata,recpos,recpos_loc,srcpos,nsrc,ns,SEIS_FORMAT,xpg,ypg,zpg);
	}
}
