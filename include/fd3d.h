/*------------------------------------------------------------------------
 *  fd3d.h - include file for program FD3D
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* declaration of functions */
void absorb(float *** absorb_coeff);

void av_mat(float *** u, float *** ua, int cc, char t, float * x, float * xp, float * y, float * yp, float * z, float * zp);

void check_par(void);

void countsrc(void);

void exchange_par(void);

void fsource(int nt, struct vector3d ***v, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk,
float ** srcpos_loc, float ** signals, int nsrc, 
float * dx, float * dy, float * dz, float * dxp, float * dyp, float * dzp);

void init_MPIshot();

void initproc(void);

void matcopy(float *** rho);

void merge(int nsnap, int type, float *x, float *y, float *z);

void mergemod(char modfile[STRING_SIZE], int format);

void  output_source_signal(FILE *fp, float **signals, float **srcpos_loc, int ns, int seis_form, 
float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg);

void  outseis(FILE *fp, FILE *fpdata, float **section,
int **recpos, int **recpos_loc, float ** srcpos_loc,
int nsrc, int ns, int seis_form, float *x, float *y, float *z);

void pml_profile(struct pml *pmlle, struct pml *pmlri, struct pml *pmlba, struct pml *pmlfr, struct pml *pmlto, struct pml *pmlbo,
float * xg, float * yg, float * zg, float * xpg, float * ypg, float * zpg);

void read_grid_glob(FILE *fp, char *dx_file, int dim, float *xpg);

void read_par(FILE *fp_in);

float **read_wavelet(float **srcpos_loc, int nsrc_loc, int *nsamp_max);

int **receiver(FILE *fp, float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg);

void shotno_glob2loc(struct acq *acq);

int **splitrec(int **recpos, int *recswitch);

float **splitsrc(float **srcpos);

float ** wavelet(float ** srcpos_loc);

void writemod(char modfile[STRING_SIZE], float *** array, int format);
