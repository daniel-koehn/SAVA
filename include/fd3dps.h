/*------------------------------------------------------------------------
 *  fdps.h - include file for program FD3D (visco-elastic)         
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* declaration of functions */
void checkfd(FILE *fp, float *** rho, float *** pi, float *** u, float *** taus, float *** taup, float *eta, 
float *x, float *y, float *z);

void fd_coeff(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, 
float * x, float * y, float * z, float * xp, float * yp, float * zp, 
float * dx, float *dxp, float * dy, float *dyp, float * dz, float * dzp,
float *** pi, float *** u, float *** uipjk, float *** uijpk, float *** uijkp, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk);

void info(FILE *fp);

void model(float *** rho, float *** pi, float *** u, float *** taus, float *** taup, float * eta, 
float * x, float * y, float * z);

double pml_update_s(struct vector3d ***v, struct tensor3d ***t,
float *** pi, float *** u, float *** uipjk, float *** uijpk, float *** uijkp, 
struct pml *pmlle, struct pml *pmlri, 
struct pml *pmlba, struct pml *pmlfr, 
struct pml *pmlto, struct pml *pmlbo, 
struct vector3d *** Pv_x_l, struct vector3d *** Pv_x_r,
struct vector3d *** Pv_y_b, struct vector3d *** Pv_y_f,
struct vector3d *** Pv_z_t, struct vector3d *** Pv_z_b,
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, int infoout);

void readmod(float *** rho, float *** pi, float *** u, float *** taus, float *** taup, float * eta);

double update_s_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, struct tensor3d ***t, struct divcurl3d ***w,
float *** pi, float *** u, float *** uipjk, float *** uijpk, float *** uijkp,
float *** absorb_coeff,
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, int infoout);

double update_s_ve(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, struct tensor3d ***t, struct tensor3d ***r, struct divcurl3d ***w,
float *** pi, float *** u, float *** uipjk, float *** uijpk, float *** uijkp,
float *** taup, float *** taus, float *** tausipjk, float *** tausijpk, float *** tausijkp, float * eta, 
float *** absorb_coeff, 
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, int infoout);

void write_par(FILE *fp);
