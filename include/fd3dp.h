/*------------------------------------------------------------------------
 *  fdp.h - include file for program FD3D (visco-acoustic)         
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* declaration of functions */
void checkfd(FILE *fp, float *** rho, float *** pi, float *** taup, float *eta, float *x, float *y, float *z);

double exchange_s(float *** p, 
float ** sbuf_rig_to_lef, float ** sbuf_fro_to_bac, 
float ** sbuf_bot_to_top, float ** rbuf_rig_to_lef, 
float ** rbuf_fro_to_bac, float ** rbuf_bot_to_top, 
MPI_Request * request,
int infoout);

double exchange_v(struct vector3d ***v, 
float ** sbuf_lef_to_rig, float ** sbuf_bac_to_fro, 
float ** sbuf_top_to_bot, float ** rbuf_lef_to_rig, 
float ** rbuf_bac_to_fro, float ** rbuf_top_to_bot, 
MPI_Request * request,
int infoout);

void FD_AC();

void fd_coeff(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, 
float * x, float * y, float * z, float * xp, float * yp, float * zp, 
float * dx, float *dxp, float * dy, float *dyp, float * dz, float * dzp,
float *** pi, float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk);

/*void forward_shot_AC(struct wave *wave, struct pmls *pmls, struct mat *mat, struct geom *geom, struct mpi *mpi, 
                     struct seis *seis, struct acq *acq, struct times *times, int nsrc, int nsrc_loc, int ns, int ntr, 
		     int ntr_loc);*/

void info(FILE *fp);

void model(float *** rho, float *** pi, float *** taup, float * eta, float * x, float * y, float * z);

double pml_update_s(struct vector3d ***v, float *** p, float *** pi, 
struct pml *pmlle, struct pml *pmlri, 
struct pml *pmlba, struct pml *pmlfr, 
struct pml *pmlto, struct pml *pmlbo, 
float *** Pvx_x_l, float *** Pvx_x_r, 
float *** Pvy_y_b, float *** Pvy_y_f, 
float *** Pvz_z_t, float *** Pvz_z_b,
float * dx, float * dy, float * dz, int infoout);

double pml_update_v(struct vector3d ***v, float *** p, float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk, 
struct pml *pmlle, struct pml *pmlri, 
struct pml *pmlba, struct pml *pmlfr, 
struct pml *pmlto, struct pml *pmlbo, 
float *** Pp_x_l, float *** Pp_x_r, 
float *** Pp_y_b, float *** Pp_y_f, 
float *** Pp_z_t, float *** Pp_z_b,
float * dxp, float * dyp, float * dzp, int infoout);

void psource(int nt, float *** p, 
float **  srcpos_loc, float ** signals, int nsrc_loc, 
float * dx, float * dy, float * dz);

void readmod(float *** rho, float *** pi, float *** taup, float * eta);

void saveseis(FILE *fp, float **sectionvx, float **sectionvy, float **sectionvz, float **sectionp, 
float **sectionax, float **sectionay, float **sectionaz, float **sectiondiv, float **section_fulldata, 
int  **recpos, int  **recpos_loc, float ** srcpos, int nsrc, int ns, float *xg, float *yg, 
float *zg, float *xpg, float *ypg, float *zpg, int *recswitch);

void savesnap(FILE *fp, float * xpg, float * ypg, float * zpg);

void seismo(int lsamp, int **recpos_loc, 
float **sectionvx, float **sectionvy, float **sectionvz, float **sectionp, 
float **sectionax, float **sectionay, float **sectionaz, float **sectiondiv, 
struct vector3d ***v, float ***p, struct vector3d ***a, float ***diverg);

void snap(FILE *fp,int nt, int nsnap, struct vector3d ***v, float ***p, struct vector3d ***a, float ***diverg, float *x, float *y, float *z);

void snapmerge(int nsnap);

float **sources(float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg);

double update_s_ac(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, float *** p, float *** diverg,
float *** pi, 
float *** absorb_coeff,
float * dx, float * dy, float * dz, int infoout);

double update_s_va(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, float *** p, float *** rpp, float *** diverg,
float *** pi, float *** taup, float * eta, 
float *** absorb_coeff, 
float * dx, float * dy, float * dz, int infoout);

double update_v(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, float *** p, struct vector3d ***a, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk,
float *** absorb_coeff,
float * dxp, float * dyp, float * dzp, int infoout);

void write_par(FILE *fp);
