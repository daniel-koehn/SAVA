/*------------------------------------------------------------------------ *   
 *   Non-splitting CFS-PML boundary condition at the numerical edge of the grid,
 *   initilize absorbing profile in the PML
 *
 *   F.H. Drossaert, A. Giannopoulos (2007), A nonsplit complex frequency-shifted
 *   PML based on recursive integration for FDTD modeling of elastic waves. 
 *   Geophysics 72, T9-T17.
 *
 *   F.H. Drossaert, A. Giannopoulos (2007), Complex frequency shifted
 *   convolution PML for FDTD modeling of elastic waves. Wave Motion 44, 593-604.
 *
 *   O. Hellwig
 *------------------------------------------------------------------------ */

#include "fd.h"

void pml_profile(struct pml *pmlle, struct pml *pmlri, struct pml *pmlba, struct pml *pmlfr, struct pml *pmlto, struct pml *pmlbo, 
		float * xg, float * yg, float * zg, float * xpg, float * ypg, float *zpg){

	/* extern variables */
	extern float	N_SIGMA_LE, SIGMA_LE, N_KAPPA_LE, KAPPA_LE, N_ALPHA_LE, ALPHA_LE;
	extern float	N_SIGMA_RI, SIGMA_RI, N_KAPPA_RI, KAPPA_RI, N_ALPHA_RI, ALPHA_RI;
	extern float	N_SIGMA_BA, SIGMA_BA, N_KAPPA_BA, KAPPA_BA, N_ALPHA_BA, ALPHA_BA;
	extern float	N_SIGMA_FR, SIGMA_FR, N_KAPPA_FR, KAPPA_FR, N_ALPHA_FR, ALPHA_FR;
	extern float	N_SIGMA_TO, SIGMA_TO, N_KAPPA_TO, KAPPA_TO, N_ALPHA_TO, ALPHA_TO;
	extern float	N_SIGMA_BO, SIGMA_BO, N_KAPPA_BO, KAPPA_BO, N_ALPHA_BO, ALPHA_BO;
	extern float	DT;
	extern int	PML_LE, PML_RI, PML_BA, PML_FR, PML_TO, PML_BO;
	extern int	NXG[3];
	extern int	FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern int	MYID;
	extern FILE	*FP;

	/* local variables */
	int	i;
	float	fwx_l, fwx_r, fwy_b, fwy_f, fwz_t, fwz_b;
	float	sigma, kappa, alpha;
	float	k, l;

	/* parameters for PML */
	/* coordinate streching: 1/(kappa+sigma/(alpha+i*omega)) */
	/* SIGMA=0, KAPPA=1, ALPHA=0 ... no PML; 
	   SIGMA>0, KAPPA=1, ALPHA=0 ... ordinary PML;
	   SIGMA>0, KAPPA>1, ALPHA>0 ... CFS-PML */
	/* SIGMA=(N+1)/2*v_PML*ln(1/R); R ... level of reflectivity
	   KAPPA>=1                       ... improves wave absorption
	   ALPHA>=0                       ... improves wave absorption */

	/* top boundary */
	/* exponent N_SIGMA_TO and maximum sigma value SIGMA_TO [m/s] */
	/* exponent N_KAPPA_TO and maximum kappa value KAPPA_TO */
	/* exponent N_ALPHA_TO and maximum alpha value ALPHA_TO [m/s] */

	/* bottom boundary */
	/* exponent N_SIGMA_BO and maximum sigma value SIGMA_BO [m/s] */
	/* exponent N_KAPPA_BO and maximum kappa value KAPPA_BO */
	/* exponent N_ALPHA_BO and maximum alpha value ALPHA_BO [m/s] */

	/* back boundary */
	/* exponent N_SIGMA_BA and maximum sigma value SIGMA_BA [m/s] */
	/* exponent N_KAPPA_BA and maximum kappa value KAPPA_BA */
	/* exponent N_ALPHA_BA and maximum alpha value ALPHA_BA [m/s] */

	/* front boundary */
	/* exponent N_SIGMA_FR and maximum sigma value SIGMA_FR [m/s] */
	/* exponent N_KAPPA_FR and maximum kappa value KAPPA_FR */
	/* exponent N_ALPHA_FR and maximum alpha value ALPHA_FR [m/s] */

	/* left boundary */
	/* exponent N_SIGMA_LE and maximum sigma value SIGMA_LE [m/s] */
	/* exponent N_KAPPA_LE and maximum kappa value KAPPA_LE */
	/* exponent N_ALPHA_LE and maximum alpha value ALPHA_LE [m/s] */

	/* right boundary */
	/* exponent N_SIGMA_RI and maximum sigma value SIGMA_RI [m/s] */
	/* exponent N_KAPPA_RI and maximum kappa value KAPPA_RI */
	/* exponent N_ALPHA_RI and maximum alpha value ALPHA_RI [m/s] */


	/* width of PML */
	fwx_l  = xg[FWLE+1]-xg[1];
	fwx_r  = xpg[NXG[0]]-xpg[NXG[0]-FWRI];
	fwy_b  = yg[FWBA+1]-yg[1];
	fwy_f  = ypg[NXG[1]]-ypg[NXG[1]-FWFR];
	fwz_t  = zg[FWTO+1]-zg[1];
	fwz_b  = zpg[NXG[2]]-zpg[NXG[2]-FWBO];

	/* maximum value of PML-coefficients */
	SIGMA_TO = SIGMA_TO/fwz_t;
	ALPHA_TO = ALPHA_TO/fwz_t;
	SIGMA_BO = SIGMA_BO/fwz_b;
	ALPHA_BO = ALPHA_BO/fwz_b;
	SIGMA_BA = SIGMA_BA/fwy_b;
	ALPHA_BA = ALPHA_BA/fwy_b;
	SIGMA_FR = SIGMA_FR/fwy_f;
	ALPHA_FR = ALPHA_FR/fwy_f;
	SIGMA_LE = SIGMA_LE/fwx_l;
	ALPHA_LE = ALPHA_LE/fwx_l;
	SIGMA_RI = SIGMA_RI/fwx_r;
	ALPHA_RI = ALPHA_RI/fwx_r;

	fprintf(FP,"\n ----------------------- PML -------------------------------------\n");
	fprintf(FP," **Message from pml_profile (printed by PE %d):\n",MYID);
	fprintf(FP," Coefficients for PML frame are now calculated.\n");
	if (PML_LE){
		fprintf(FP," PML parameters applied at left boundary\n");
		fprintf(FP," Width of PML boundary (grid points)= %d\n\n",FWLE);
		fprintf(FP," Width of PML boundary L=%e m \n",fwx_l);
		fprintf(FP," n=%e, sigma_max/L=%e Hz \n",N_SIGMA_LE,SIGMA_LE);
		fprintf(FP," n=%e, kappa_max  =%e    \n",N_KAPPA_LE,KAPPA_LE);
		fprintf(FP," n=%e, alpha_max/L=%e Hz \n\n",N_ALPHA_LE,ALPHA_LE);

		/* CFS-PML at left boundary */
		for (i=1;i<=FWLE;i++){
			k = ((float)i-0.5)/(float)FWLE;
			l = 1.0-k;

			sigma = SIGMA_LE*pow(k,N_SIGMA_LE);
			kappa = (KAPPA_LE-1.0)*pow(k,N_KAPPA_LE)+1.0;
			alpha = ALPHA_LE*pow(l,N_ALPHA_LE);
			//sigma = SIGMA_LE*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_LE-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_LE*0.5*(1.0+cos(PI*k));

			pmlle[i].xp2  = 1.0/kappa - 1.0;
			pmlle[i].expp = exp(-DT*(sigma/kappa+alpha));
			pmlle[i].xp1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlle[i].expp-1.0);	/* PML version 1 */
			//pmlle[i].xp1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */


			k = (float)i/(float)FWLE;
			l = 1.0-k;

			sigma = SIGMA_LE*pow(k,N_SIGMA_LE);
			kappa = (KAPPA_LE-1.0)*pow(k,N_KAPPA_LE)+1.0;
			alpha = ALPHA_LE*pow(l,N_ALPHA_LE);
			//igma = SIGMA_LE*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_LE-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_LE*0.5*(1.0+cos(PI*k));

			pmlle[i].x2  = 1.0/kappa - 1.0;
			pmlle[i].exp = exp(-DT*(sigma/kappa+alpha));
			pmlle[i].x1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlle[i].exp-1.0);	/* PML version 1 */
			//pmlle[i].x1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */
		}
		fprintf(FP,"\n PML damping profile in x-direction (left boundary)\n");
		fprintf(FP," # \t DT*sigma/kappa^2 \t (1-kappa)/kappa \t exp(-DT*(sigma/kappa+alpha))\n");
		for (i=1;i<=FWLE;i++)
			fprintf(FP," %d \t %e \t %e \t %e\n", i, 2.0*pmlle[i].x1, pmlle[i].x2, pmlle[i].exp);
	}
	if (PML_RI){
		fprintf(FP," PML parameters applied at right boundary\n");
		fprintf(FP," Width of PML boundary (grid points)= %d\n\n",FWRI);
		fprintf(FP," Width of PML boundary L=%e m \n",fwx_r);
		fprintf(FP," n=%e, sigma_max/L=%e Hz \n",N_SIGMA_RI,SIGMA_RI);
		fprintf(FP," n=%e, kappa_max  =%e    \n",N_KAPPA_RI,KAPPA_RI);
		fprintf(FP," n=%e, alpha_max/L=%e Hz \n\n",N_ALPHA_RI,ALPHA_RI);

		/* CFS-PML at right boundary */
		for (i=1;i<=FWRI;i++){
			k = ((float)i-0.5)/(float)FWRI;
			l = 1.0-k;

			sigma = SIGMA_RI*pow(k,N_SIGMA_RI);
			kappa = (KAPPA_RI-1.0)*pow(k,N_KAPPA_RI)+1.0;
			alpha = ALPHA_RI*pow(l,N_ALPHA_RI);
			//sigma = SIGMA_RI*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_RI-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_RI*0.5*(1.0+cos(PI*k));

			pmlri[i].x2  = 1.0/kappa - 1.0;
			pmlri[i].exp = exp(-DT*(sigma/kappa+alpha));
			pmlri[i].x1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlri[i].exp-1.0);	/* PML version 1 */
			//pmlri[i].x1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */


			k = (float)i/(float)FWRI;
			l = 1.0-k;

			sigma = SIGMA_RI*pow(k,N_SIGMA_RI);
			kappa = (KAPPA_RI-1.0)*pow(k,N_KAPPA_RI)+1.0;
			alpha = ALPHA_RI*pow(l,N_ALPHA_RI);
			//sigma = SIGMA_RI*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_RI-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_RI*0.5*(1.0+cos(PI*k));//

			pmlri[i].xp2  = 1.0/kappa - 1.0;
			pmlri[i].expp = exp(-DT*(sigma/kappa+alpha));
			pmlri[i].xp1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlri[i].expp-1.0);	/* PML version 1 */
			//pmlri[i].xp1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */
		}
		fprintf(FP,"\n PML damping profile in x-direction (right boundary)\n");
		fprintf(FP," # \t DT*sigma/kappa^2 \t (1-kappa)/kappa \t exp(-DT*(sigma/kappa+alpha))\n");
		for (i=1;i<=FWRI;i++)
			fprintf(FP," %d \t %e \t %e \t %e\n", i, 2.0*pmlri[i].x1, pmlri[i].x2, pmlri[i].exp);
	}
	if (PML_BA){
		fprintf(FP," PML parameters applied at back boundary\n");
		fprintf(FP," Width of PML boundary (grid points)= %d\n\n",FWBA);
		fprintf(FP," Width of PML boundary L=%e m \n",fwy_b);
		fprintf(FP," n=%e, sigma_max/L=%e Hz \n",N_SIGMA_BA,SIGMA_BA);
		fprintf(FP," n=%e, kappa_max  =%e    \n",N_KAPPA_BA,KAPPA_BA);
		fprintf(FP," n=%e, alpha_max/L=%e Hz \n\n",N_ALPHA_BA,ALPHA_BA);

		/* CFS-PML at back boundary */
		for (i=1;i<=FWBA;i++){
			k = ((float)i-0.5)/(float)FWBA;
			l = 1.0-k;

			sigma = SIGMA_BA*pow(k,N_SIGMA_BA);
			kappa = (KAPPA_BA-1.0)*pow(k,N_KAPPA_BA)+1.0;
			alpha = ALPHA_BA*pow(l,N_ALPHA_BA);
			//sigma = SIGMA_BA*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_BA-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_BA*0.5*(1.0+cos(PI*k));*/

			pmlba[i].xp2  = 1.0/kappa - 1.0;
			pmlba[i].expp = exp(-DT*(sigma/kappa+alpha));
			pmlba[i].xp1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlba[i].expp-1.0);	/* PML version 1 */
			//pmlba[i].xp1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */


			k = (float)i/(float)FWBA;
			l = 1.0-k;

			sigma = SIGMA_BA*pow(k,N_SIGMA_BA);
			kappa = (KAPPA_BA-1.0)*pow(k,N_KAPPA_BA)+1.0;
			alpha = ALPHA_BA*pow(l,N_ALPHA_BA);
			//sigma = SIGMA_BA*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_BA-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_BA*0.5*(1.0+cos(PI*k));

			pmlba[i].x2  = 1.0/kappa - 1.0;
			pmlba[i].exp = exp(-DT*(sigma/kappa+alpha));
			pmlba[i].x1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlba[i].exp-1.0);	/* PML version 1 */
			//pmlba[i].x1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */
		}
		fprintf(FP,"\n PML damping profile in y-direction (back boundary)\n");
		fprintf(FP," # \t DT*sigma/kappa^2 \t (1-kappa)/kappa \t exp(-DT*(sigma/kappa+alpha))\n");
		for (i=1;i<=FWBA;i++)
			fprintf(FP," %d \t %e \t %e \t %e\n", i, 2.0*pmlba[i].x1, pmlba[i].x2, pmlba[i].exp);
	}
	if (PML_FR){
		fprintf(FP," PML parameters applied at front boundary\n");
		fprintf(FP," Width of PML boundary (grid points)= %d\n\n",FWFR);
		fprintf(FP," Width of PML boundary L=%e m \n",fwy_f);
		fprintf(FP," n=%e, sigma_max/L=%e Hz \n",N_SIGMA_FR,SIGMA_FR);
		fprintf(FP," n=%e, kappa_max  =%e    \n",N_KAPPA_FR,KAPPA_FR);
		fprintf(FP," n=%e, alpha_max/L=%e Hz \n\n",N_ALPHA_FR,ALPHA_FR);

		/* CFS-PML at front boundary */
		for (i=1;i<=FWFR;i++){
			k = ((float)i-0.5)/(float)FWFR;
			l = 1.0-k;

			sigma = SIGMA_FR*pow(k,N_SIGMA_FR);
			kappa = (KAPPA_FR-1.0)*pow(k,N_KAPPA_FR)+1.0;
			alpha = ALPHA_FR*pow(l,N_ALPHA_FR);
			//sigma = SIGMA_FR*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_FR-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_FR*0.5*(1.0+cos(PI*k));

			pmlfr[i].x2  = 1.0/kappa - 1.0;
			pmlfr[i].exp = exp(-DT*(sigma/kappa+alpha));
			pmlfr[i].x1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlfr[i].exp-1.0);	/* PML version 1 */
			//pmlfr[i].x1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */


			k = (float)i/(float)FWFR;
			l = 1.0-k;

			sigma = SIGMA_FR*pow(k,N_SIGMA_FR);
			kappa = (KAPPA_FR-1.0)*pow(k,N_KAPPA_FR)+1.0;
			alpha = ALPHA_FR*pow(l,N_ALPHA_FR);
			//sigma = SIGMA_FR*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_FR-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_FR*0.5*(1.0+cos(PI*k));

			pmlfr[i].xp2  = 1.0/kappa - 1.0;
			pmlfr[i].expp = exp(-DT*(sigma/kappa+alpha));
			pmlfr[i].xp1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlfr[i].expp-1.0);	/* PML version 1 */
			//pmlfr[i].xp1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */
		}
		fprintf(FP,"\n PML damping profile in y-direction (front boundary)\n");
		fprintf(FP," # \t DT*sigma/kappa^2 \t (1-kappa)/kappa \t exp(-DT*(sigma/kappa+alpha))\n");
		for (i=1;i<=FWFR;i++)
			fprintf(FP," %d \t %e \t %e \t %e\n", i, 2.0*pmlfr[i].x1, pmlfr[i].x2, pmlfr[i].exp);
	}
	if (PML_TO){
		fprintf(FP," PML parameters applied at top boundary\n");
		fprintf(FP," Width of PML boundary (grid points)= %d\n\n",FWTO);
		fprintf(FP," Width of PML boundary L=%e m \n",fwz_t);
		fprintf(FP," n=%e, sigma_max/L=%e Hz \n",N_SIGMA_TO,SIGMA_TO);
		fprintf(FP," n=%e, kappa_max  =%e    \n",N_KAPPA_TO,KAPPA_TO);
		fprintf(FP," n=%e, alpha_max/L=%e Hz \n\n",N_ALPHA_TO,ALPHA_TO);

		/* CFS-PML at top boundary */
		for (i=1;i<=FWTO;i++){
			k = ((float)i-0.5)/(float)FWTO;
			l = 1.0-k;

			sigma = SIGMA_TO*pow(k,N_SIGMA_TO);
			kappa = (KAPPA_TO-1.0)*pow(k,N_KAPPA_TO)+1.0;
			alpha = ALPHA_TO*pow(l,N_ALPHA_TO);
			//sigma = SIGMA_TO*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_TO-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_TO*0.5*(1.0+cos(PI*k));

			pmlto[i].xp2  = 1.0/kappa - 1.0;
			pmlto[i].expp = exp(-DT*(sigma/kappa+alpha));
			pmlto[i].xp1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlto[i].expp-1.0);	/* PML version 1 */
			//pmlto[i].xp1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */


			k = (float)i/(float)FWTO;
			l = 1.0-k;

			sigma = SIGMA_TO*pow(k,N_SIGMA_TO);
			kappa = (KAPPA_TO-1.0)*pow(k,N_KAPPA_TO)+1.0;
			alpha = ALPHA_TO*pow(l,N_ALPHA_TO);
			//sigma = SIGMA_TO*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_TO-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_TO*0.5*(1.0+cos(PI*k));

			pmlto[i].x2  = 1.0/kappa - 1.0;
			pmlto[i].exp = exp(-DT*(sigma/kappa+alpha));
			pmlto[i].x1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlto[i].exp-1.0);	/* PML version 1 */ 
			//pmlto[i].x1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */
		}
		fprintf(FP," PML damping profile in z-direction (top boundary)\n");
		fprintf(FP," # \t DT*sigma/kappa^2 \t (1-kappa)/kappa \t exp(-DT*(sigma/kappa+alpha))\n");
		for (i=1;i<=FWTO;i++)
			fprintf(FP," %d \t %e \t %e \t %e\n", i, 2.0*pmlto[i].x1, pmlto[i].x2, pmlto[i].exp);
	}
	if (PML_BO){
		fprintf(FP," PML parameters applied at bottom boundary\n");
		fprintf(FP," Width of PML boundary (grid points)= %d\n\n",FWBO);
		fprintf(FP," Width of PML boundary L=%e m \n",fwz_b);
		fprintf(FP," n=%e, sigma_max/L=%e Hz \n",N_SIGMA_BO,SIGMA_BO);
		fprintf(FP," n=%e, kappa_max  =%e    \n",N_KAPPA_BO,KAPPA_BO);
		fprintf(FP," n=%e, alpha_max/L=%e Hz \n\n",N_ALPHA_BO,ALPHA_BO);

		/* CFS-PML at bottom boundary */
		for (i=1;i<=FWBO;i++){
			k = ((float)i-0.5)/(float)FWBO;
			l = 1.0-k;

			sigma = SIGMA_BO*pow(k,N_SIGMA_BO);
			kappa = (KAPPA_BO-1.0)*pow(k,N_KAPPA_BO)+1.0;
			alpha = ALPHA_BO*pow(l,N_ALPHA_BO);
			//sigma = SIGMA_BO*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_BO-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_BO*0.5*(1.0+cos(PI*k));

			pmlbo[i].x2  = 1.0/kappa - 1.0;
			pmlbo[i].exp = exp(-DT*(sigma/kappa+alpha));
			pmlbo[i].x1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlbo[i].exp-1.0);	/* PML version 1 */
			//pmlbo[i].x1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */


			k = (float)i/(float)FWBO;
			l = 1.0-k;

			sigma = SIGMA_BO*pow(k,N_SIGMA_BO);
			kappa = (KAPPA_BO-1.0)*pow(k,N_KAPPA_BO)+1.0;
			alpha = ALPHA_BO*pow(l,N_ALPHA_BO);
			//sigma = SIGMA_BO*0.5*(1.0-cos(PI*k));
			//kappa = (KAPPA_BO-1.0)*0.5*(1.0-cos(PI*k))+1.0;
			//alpha = ALPHA_BO*0.5*(1.0+cos(PI*k));

			pmlbo[i].xp2  = 1.0/kappa - 1.0;
			pmlbo[i].expp = exp(-DT*(sigma/kappa+alpha));
			pmlbo[i].xp1  = sigma/(kappa*(sigma+kappa*alpha))*(pmlbo[i].expp-1.0);	/* PML version 1 */
			//pmlbo[i].xp1  = -0.5*DT*sigma/(kappa*kappa);	/* PML version 2 */
		}
		fprintf(FP,"\n PML damping profile in z-direction (bottom boundary)\n");
		fprintf(FP," # \t DT*sigma/kappa^2 \t (1-kappa)/kappa \t exp(-DT*(sigma/kappa+alpha))\n");
		for (i=1;i<=FWBO;i++)
			fprintf(FP," %d \t %e \t %e \t %e\n", i, 2.0*pmlbo[i].x1, pmlbo[i].x2, pmlbo[i].exp);
	}

}
