/*------------------------------------------------------------------------
 *   Read elastic model properties (Cijkl, density) from files 
 *
 *   D. Köhn 
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "fd3d.h"

void readmod(float *** rho, 
	     float *** c1111, float *** c1112, float *** c1113, float *** c1122, float *** c1133, float *** c1123, float *** c1212, 
	     float *** c1213, float *** c1222, float *** c1223, float *** c1233, float *** c1313, float *** c1322, float *** c1323, 
	     float *** c1333, float *** c2222, float *** c2223, float *** c2233, float *** c2323, float *** c2333, float *** c3333,
	     float *** taus, float *** taup, float * eta){

	/* extern variables */
	extern float	DT, *FL, TAU;
	extern int	NX[3], NXG[3], POS[3], L, MYID;
	extern int	READMOD;
	extern int	MODEL, MODEL_FORMAT;
	extern char	MFILE[STRING_SIZE], MODEL_FILE[STRING_SIZE];
	extern FILE	*FP;

	/* local variables */
	float	rhov;
	float	cv1111, cv1112, cv1113, cv1122, cv1133, cv1123, cv1212;
	float	cv1213, cv1222, cv1223, cv1233, cv1313, cv1322, cv1323;
	float	cv1333, cv2222, cv2223, cv2233, cv2323, cv2333, cv3333;
	float	qp, qs;
	float	ts, tp, sumu, sumpi, ws;
	int	i, j, k, l, ii, jj, kk;
	char	filename[STRING_SIZE], ext[8];
	FILE	*fp_rho;
	FILE	*fp_c1111, *fp_c1112, *fp_c1113, *fp_c1122, *fp_c1133, *fp_c1123, *fp_c1212;
	FILE	*fp_c1213, *fp_c1222, *fp_c1223, *fp_c1233, *fp_c1313, *fp_c1322, *fp_c1323;
	FILE	*fp_c1333, *fp_c2222, *fp_c2223, *fp_c2233, *fp_c2323, *fp_c2333, *fp_c3333;
	FILE	*fp_qp=NULL, *fp_qs=NULL;

	/*-----------------------------------------------------------------------*/

	fprintf(FP,"\n **Message from function readmod (printed by PE %d):\n",MYID);

	/* vector for maxwellbodies */
	tp = TAU;
	ts = TAU;

	sumpi = 0.0;
	sumu  = 0.0;
	for (l=1;l<=L;l++){
		eta[l] = 2.0*PI*DT*FL[l];
		ws     = FL[l]/FL[1];
		sumpi += 1.0/(1.0+ws*ws);
		sumu  += 1.0/(1.0+ws*ws);
	}
	sumpi *= tp;
	sumu  *= ts;


	fprintf(FP," Reading model information from model-files...\n");

	sprintf(filename,"%s.c1111",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1111=fopen(filename,"r");
	if (fp_c1111==NULL) error(" Could not open model file for component c1111 ! ");

	sprintf(filename,"%s.c1112",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1112=fopen(filename,"r");
	if (fp_c1112==NULL) error(" Could not open model file for component c1112 ! ");

	sprintf(filename,"%s.c1113",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1113=fopen(filename,"r");
	if (fp_c1113==NULL) error(" Could not open model file for component c1113 ! ");

	sprintf(filename,"%s.c1122",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1122=fopen(filename,"r");
	if (fp_c1122==NULL) error(" Could not open model file for component c1122 ! ");

	sprintf(filename,"%s.c1133",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1133=fopen(filename,"r");
	if (fp_c1133==NULL) error(" Could not open model file for component c1133 ! ");

	sprintf(filename,"%s.c1123",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1123=fopen(filename,"r");
	if (fp_c1123==NULL) error(" Could not open model file for component c1123 ! ");

	sprintf(filename,"%s.c1212",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1212=fopen(filename,"r");
	if (fp_c1212==NULL) error(" Could not open model file for component c1212 ! ");

	sprintf(filename,"%s.c1213",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1213=fopen(filename,"r");
	if (fp_c1213==NULL) error(" Could not open model file for component c1213 ! ");

	sprintf(filename,"%s.c1222",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1222=fopen(filename,"r");
	if (fp_c1222==NULL) error(" Could not open model file for component c1222 ! ");

	sprintf(filename,"%s.c1223",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1223=fopen(filename,"r");
	if (fp_c1223==NULL) error(" Could not open model file for component c1223 ! ");

	sprintf(filename,"%s.c1233",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1233=fopen(filename,"r");
	if (fp_c1233==NULL) error(" Could not open model file for component c1233 ! ");

	sprintf(filename,"%s.c1313",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1313=fopen(filename,"r");
	if (fp_c1313==NULL) error(" Could not open model file for component c1313 ! ");

	sprintf(filename,"%s.c1322",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1322=fopen(filename,"r");
	if (fp_c1322==NULL) error(" Could not open model file for component c1322 ! ");

	sprintf(filename,"%s.c1323",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1323=fopen(filename,"r");
	if (fp_c1323==NULL) error(" Could not open model file for component c1323 ! ");

	sprintf(filename,"%s.c1333",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c1333=fopen(filename,"r");
	if (fp_c1333==NULL) error(" Could not open model file for component c1333 ! ");

	sprintf(filename,"%s.c2222",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c2222=fopen(filename,"r");
	if (fp_c2222==NULL) error(" Could not open model file for component c2222 ! ");

	sprintf(filename,"%s.c2223",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c2223=fopen(filename,"r");
	if (fp_c2223==NULL) error(" Could not open model file for component c2223 ! ");

	sprintf(filename,"%s.c2233",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c2233=fopen(filename,"r");
	if (fp_c2233==NULL) error(" Could not open model file for component c2233 ! ");

	sprintf(filename,"%s.c2323",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c2323=fopen(filename,"r");
	if (fp_c2323==NULL) error(" Could not open model file for component c2323 ! ");

	sprintf(filename,"%s.c2333",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c2333=fopen(filename,"r");
	if (fp_c2333==NULL) error(" Could not open model file for component c2333 ! ");

	sprintf(filename,"%s.c3333",MFILE);
	fprintf(FP,"\t Elasticity tensor:\n\t %s \n\n",filename);
	fp_c3333=fopen(filename,"r");
	if (fp_c3333==NULL) error(" Could not open model file for component c3333 ! ");

	sprintf(filename,"%s.rho",MFILE);
	fprintf(FP,"\t Density:\n\t %s \n\n",filename);
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) error(" Could not open model file for densities ! ");

	if (L){
		sprintf(filename,"%s.qp",MFILE);
		fprintf(FP,"\t Qp:\n\t %s \n\n",filename);
		fp_qp=fopen(filename,"r");
		if (fp_qp==NULL) error(" Could not open model file for Qp-values ! ");

		sprintf(filename,"%s.qs",MFILE);
		fprintf(FP,"\t Qs:\n\t %s \n\n",filename);
		fp_qs=fopen(filename,"r");
		if (fp_qs==NULL) error(" Could not open model file for Qs-values ! ");
	}

	/* loop over global grid */
	for (i=1;i<=NXG[0];i++){
		if (!(i%100))
			fprintf(FP,"\n Reading model: step %d of %d \n",i,NXG[0]);
		for (j=1;j<=NXG[1];j++){
			if (READMOD==2){
				/* jump over 240 byte SU header */
				fseek(fp_rho,240,SEEK_CUR);
				fseek(fp_c1111,240,SEEK_CUR);
				fseek(fp_c1112,240,SEEK_CUR);
				fseek(fp_c1113,240,SEEK_CUR);
				fseek(fp_c1122,240,SEEK_CUR);
				fseek(fp_c1133,240,SEEK_CUR);
				fseek(fp_c1123,240,SEEK_CUR);
				fseek(fp_c1212,240,SEEK_CUR);
				fseek(fp_c1213,240,SEEK_CUR);
				fseek(fp_c1222,240,SEEK_CUR);
				fseek(fp_c1223,240,SEEK_CUR);
				fseek(fp_c1233,240,SEEK_CUR);
				fseek(fp_c1313,240,SEEK_CUR);
				fseek(fp_c1322,240,SEEK_CUR);
				fseek(fp_c1323,240,SEEK_CUR);
				fseek(fp_c1333,240,SEEK_CUR);
				fseek(fp_c2222,240,SEEK_CUR);
				fseek(fp_c2223,240,SEEK_CUR);
				fseek(fp_c2233,240,SEEK_CUR);
				fseek(fp_c2323,240,SEEK_CUR);
				fseek(fp_c2333,240,SEEK_CUR);
				fseek(fp_c3333,240,SEEK_CUR);
				if (L){
					fseek(fp_qp,240,SEEK_CUR);
					fseek(fp_qs,240,SEEK_CUR);
				}
			}
			for (k=1;k<=NXG[2];k++){
				fread(&cv1111, sizeof(float), 1, fp_c1111);
				fread(&cv1112, sizeof(float), 1, fp_c1112);
				fread(&cv1113, sizeof(float), 1, fp_c1113);
				fread(&cv1122, sizeof(float), 1, fp_c1122);
				fread(&cv1133, sizeof(float), 1, fp_c1133);
				fread(&cv1123, sizeof(float), 1, fp_c1123);
				fread(&cv1212, sizeof(float), 1, fp_c1212);
				fread(&cv1213, sizeof(float), 1, fp_c1213);
				fread(&cv1222, sizeof(float), 1, fp_c1222);
				fread(&cv1223, sizeof(float), 1, fp_c1223);
				fread(&cv1233, sizeof(float), 1, fp_c1233);
				fread(&cv1313, sizeof(float), 1, fp_c1313);
				fread(&cv1322, sizeof(float), 1, fp_c1322);
				fread(&cv1323, sizeof(float), 1, fp_c1323);
				fread(&cv1333, sizeof(float), 1, fp_c1333);
				fread(&cv2222, sizeof(float), 1, fp_c2222);
				fread(&cv2223, sizeof(float), 1, fp_c2223);
				fread(&cv2233, sizeof(float), 1, fp_c2233);
				fread(&cv2323, sizeof(float), 1, fp_c2323);
				fread(&cv2333, sizeof(float), 1, fp_c2333);
				fread(&cv3333, sizeof(float), 1, fp_c3333);
				fread(&rhov, sizeof(float), 1, fp_rho);
				if (L){
					fread(&qp, sizeof(float), 1, fp_qp);
					fread(&qs, sizeof(float), 1, fp_qs);
				}

				cv1111=cv1111/(1.0+sumpi);
				cv1112=cv1112/(1.0+sumu); /* ??? */
				cv1113=cv1113/(1.0+sumu); /* ??? */
				cv1122=cv1122/(1.0+sumpi);
				cv1133=cv1133/(1.0+sumpi);
				cv1123=cv1123/(1.0+sumu); /* ??? */
				cv1212=cv1212/(1.0+sumu);
				cv1213=cv1213/(1.0+sumu);
				cv1222=cv1222/(1.0+sumu); /* ??? */
				cv1223=cv1223/(1.0+sumu);
				cv1233=cv1233/(1.0+sumu); /* ??? */
				cv1313=cv1313/(1.0+sumu);
				cv1322=cv1322/(1.0+sumu); /* ??? */
				cv1323=cv1323/(1.0+sumu);
				cv1333=cv1333/(1.0+sumu); /* ??? */
				cv2222=cv2222/(1.0+sumpi);
				cv2223=cv2223/(1.0+sumu); /* ??? */
				cv2233=cv2233/(1.0+sumpi);
				cv2323=cv2323/(1.0+sumu);
				cv2333=cv2333/(1.0+sumu); /* ??? */
				cv3333=cv3333/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint 
				is saving model parameters in its local arrays */
				if ((POS[0]==((i-1)/NX[0])) && (POS[1]==((j-1)/NX[1])) && (POS[2]==((k-1)/NX[2]))){
					ii = i-POS[0]*NX[0];
					jj = j-POS[1]*NX[1];
					kk = k-POS[2]*NX[2];

					if (L){
						taus[ii][jj][kk]=2.0/qs;
						taup[ii][jj][kk]=2.0/qp;
					}
					c1111[ii][jj][kk] = cv1111;
					c1112[ii][jj][kk] = cv1112;
					c1113[ii][jj][kk] = cv1113;
					c1122[ii][jj][kk] = cv1122;
					c1133[ii][jj][kk] = cv1133;
					c1123[ii][jj][kk] = cv1123;
					c1212[ii][jj][kk] = cv1212;
					c1213[ii][jj][kk] = cv1213;
					c1222[ii][jj][kk] = cv1222;
					c1223[ii][jj][kk] = cv1223;
					c1233[ii][jj][kk] = cv1233;
					c1313[ii][jj][kk] = cv1313;
					c1322[ii][jj][kk] = cv1322;
					c1323[ii][jj][kk] = cv1323;
					c1333[ii][jj][kk] = cv1333;
					c2222[ii][jj][kk] = cv2222;
					c2223[ii][jj][kk] = cv2223;
					c2233[ii][jj][kk] = cv2233;
					c2323[ii][jj][kk] = cv2323;
					c2333[ii][jj][kk] = cv2333;
					c3333[ii][jj][kk] = cv3333;
					rho[ii][jj][kk]   = rhov;
				}
			}
		}
	}

	fclose(fp_c1111);
	fclose(fp_c1112);
	fclose(fp_c1113);
	fclose(fp_c1122);
	fclose(fp_c1133);
	fclose(fp_c1123);
	fclose(fp_c1212);
	fclose(fp_c1213);
	fclose(fp_c1222);
	fclose(fp_c1223);
	fclose(fp_c1233);
	fclose(fp_c1313);
	fclose(fp_c1322);
	fclose(fp_c1323);
	fclose(fp_c1333);
	fclose(fp_c2222);
	fclose(fp_c2223);
	fclose(fp_c2233);
	fclose(fp_c2323);
	fclose(fp_c2333);
	fclose(fp_c3333);

	fclose(fp_rho);
	if (L){
		fclose(fp_qp);
		fclose(fp_qs);
	}

	/* write model to disk */

	/* different data formats of output:
	   MODEL_FORMAT=1:  SU (native byte order)
	   MODEL_FORMAT=2:  ASCII
	   MODEL_FORMAT=3:  BINARY (native byte order) */

	switch(MODEL_FORMAT){
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

	if (MODEL & 4){
		sprintf(filename,"%s_c1212%s",MODEL_FILE,ext);

		writemod(filename,c1212,MODEL_FORMAT);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!(MYID))
			mergemod(filename,MODEL_FORMAT);
	}
	if (MODEL & 2){
		sprintf(filename,"%s_c1111%s",MODEL_FILE,ext);

		writemod(filename,c1111,MODEL_FORMAT);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!(MYID))
			mergemod(filename,MODEL_FORMAT);
	}
	if (MODEL & 1){
		sprintf(filename,"%s_rho%s",MODEL_FILE,ext);

		writemod(filename,rho,MODEL_FORMAT);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!(MYID))
			mergemod(filename,MODEL_FORMAT);
	}

}




