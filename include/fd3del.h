/*------------------------------------------------------------------------
 *  fdel.h - include file for program FD3D (visco-elastic)         
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* declaration of functions */
double exchange_s(struct tensor3d ***t,
float *** sbuf_lef_to_rig, float ** sbuf_rig_to_lef, 
float *** sbuf_bac_to_fro, float ** sbuf_fro_to_bac, 
float *** sbuf_top_to_bot, float ** sbuf_bot_to_top, 
float *** rbuf_lef_to_rig, float ** rbuf_rig_to_lef, 
float *** rbuf_bac_to_fro, float ** rbuf_fro_to_bac, 
float *** rbuf_top_to_bot, float ** rbuf_bot_to_top, 
MPI_Request * request,
int infoout);

double exchange_v(struct vector3d ***v,
float ** sbuf_lef_to_rig, float *** sbuf_rig_to_lef, 
float ** sbuf_bac_to_fro, float *** sbuf_fro_to_bac, 
float ** sbuf_top_to_bot, float *** sbuf_bot_to_top, 
float ** rbuf_lef_to_rig, float *** rbuf_rig_to_lef, 
float ** rbuf_bac_to_fro, float *** rbuf_fro_to_bac, 
float ** rbuf_top_to_bot, float *** rbuf_bot_to_top, 
MPI_Request * request,
int infoout);

void FD_ISO();

double pml_update_v(struct vector3d ***v, struct tensor3d ***t, float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk, 
struct pml *pmlle, struct pml *pmlri, 
struct pml *pmlba, struct pml *pmlfr, 
struct pml *pmlto, struct pml *pmlbo, 
struct vector3d *** Pt_x_l, struct vector3d *** Pt_x_r,
struct vector3d *** Pt_y_b, struct vector3d *** Pt_y_f,
struct vector3d *** Pt_z_t, struct vector3d *** Pt_z_b,
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, int infoout);

void psource(int nt, struct tensor3d ***t,
float ** srcpos_loc, float ** signals, int nsrc_loc, 
float * dx, float * dy, float * dz, float * dxp, float * dyp, float * dzp);

void saveseis(FILE *fp, float **sectionvx, float **sectionvy, float **sectionvz, float **sectionax, float **sectionay, float **sectionaz, 
float **sectiondiv, float **sectioncurlx, float **sectioncurly, float **sectioncurlz, float **sectionp,
float **sectiontxx, float **sectiontxy, float **sectiontxz, float **sectiontyy, float **sectiontyz, float **sectiontzz, float **section_fulldata,
int  **recpos, int  **recpos_loc, float ** srcpos, int nsrc, int ns, float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg, 
int *recswitch);

void savesnap(FILE *fp, float * xpg, float * ypg, float * zpg);

void seismo(int lsamp, int ntr_loc, int **recpos_loc, 
float **sectionvx, float **sectionvy, float **sectionvz, float **sectionax, float **sectionay, float ** sectionaz, 
float **sectiondiv, float **sectioncurlx, float **sectioncurly, float **sectioncurlz, float **sectionp,
float **sectiontxx, float **sectiontxy, float ** sectiontxz, float ** sectiontyy, float ** sectiontyz, float ** sectiontzz, 
struct vector3d ***v, struct tensor3d ***t, struct vector3d ***a, struct divcurl3d ***w);

void snap(FILE *fp,int nt, int nsnap, struct vector3d ***v, struct tensor3d ***t, struct vector3d ***a, struct divcurl3d ***w, 
float *x, float *y, float *z);

void snapmerge(int nsnap);

float **sources(float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg);

double update_v(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, struct tensor3d ***t, struct vector3d ***a, float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk,
float *** absorb_coeff,
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, int infoout);
