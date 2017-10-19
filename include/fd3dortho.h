/*------------------------------------------------------------------------
 *  fdortho.h - include file for program FD3D (elastic, orthotropic)
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* declaration of functions */

void checkfd(FILE *fp, float *** rho,
float *** c1111, float *** c1122, float *** c1133, float *** c2222, float *** c2233, float *** c3333,
float *** c2323, float *** c1313, float *** c1212,
float *** taus, float *** taup, float * eta, float * x, float * y, float * z);

void fd_coeff(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, 
float * x, float * y, float * z, float * xp, float * yp, float * zp, 
float * dx, float *dxp, float * dy, float *dyp, float * dz, float * dzp,
float *** c1111,  float *** c1122,  float *** c1133, float *** c2222, float *** c2233, float *** c3333,
float *** c2323h, float *** c1313h, float *** c1212h, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk);

void FD_ORTHO();

void info(FILE *fp);

void model(float *** rho, 
float *** c1111, float *** c1122, float *** c1133, float *** c2222, float *** c2233, float *** c3333, 
float *** c2323, float *** c1313, float *** c1212, 
float *** taus, float *** taup, float * eta, float * x, float * y, float * z);

void readmod(float *** rho, 
float *** c1111, float *** c1122, float *** c1133, float *** c2222, float *** c2233, float *** c3333, 
float *** c2323, float *** c1313, float *** c1212, 
float *** taus, float *** taup, float * eta);

void timing(double * time_avg,  double * time_std, int ishot);

double update_s_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct tensor3d ***e, struct tensor3d ***t,
float *** c1111,  float *** c1122,  float *** c1133,  
                  float *** c2222,  float *** c2233,  
                                    float *** c3333,  
float *** c2323h, float *** c1313h, float *** c1212h, 
float *** absorb_coeff,
int infoout);

void write_par(FILE *fp);
