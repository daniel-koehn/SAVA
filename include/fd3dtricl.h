/*------------------------------------------------------------------------
 *  fdaniso.h - include file for program FD3D (elastic, anisotropic)
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* declaration of functions */

void checkfd(FILE *fp, float *** rho,
float *** c1111, float *** c1112, float *** c1113, float *** c1122, float *** c1133, float *** c1123, float *** c1212,
float *** c1213, float *** c1222, float *** c1223, float *** c1233, float *** c1313, float *** c1322, float *** c1323,
float *** c1333, float *** c2222, float *** c2223, float *** c2233, float *** c2323, float *** c2333, float *** c3333,
float *** taus, float *** taup, float * eta, float * x, float * y, float * z);

double exchange_e(struct tensor3d ***e,
float *** sbuf_lef_to_rig, float *** sbuf_rig_to_lef, 
float *** sbuf_bac_to_fro, float *** sbuf_fro_to_bac, 
float *** sbuf_top_to_bot, float *** sbuf_bot_to_top, 
float *** rbuf_lef_to_rig, float *** rbuf_rig_to_lef, 
float *** rbuf_bac_to_fro, float *** rbuf_fro_to_bac, 
float *** rbuf_top_to_bot, float *** rbuf_bot_to_top, 
MPI_Request * request,
int infoout);

void fd_coeff(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, 
float * x, float * y, float * z, float * xp, float * yp, float * zp, 
float * dx, float *dxp, float * dy, float *dyp, float * dz, float * dzp,
float * wxp1, float * wx1, float * wxp2, float * wx2, 
float * wyp1, float * wy1, float * wyp2, float * wy2,
float * wzp1, float * wz1, float * wzp2, float * wz2,
float *** c1111,  float *** c1122,  float *** c1133,  float *** c1123,  float *** c1113,   float *** c1112, 
	          float *** c2222,  float *** c2233,  float *** c2223,  float *** c1322,   float *** c1222, 
	                            float *** c3333,  float *** c2333,  float *** c1333,   float *** c1233, 
float *** c1123h, float *** c2223h, float *** c2333h, float *** c2323h, float *** c1323ha, float *** c1223ha, 
float *** c1113h, float *** c1322h, float *** c1333h, float *** c1323h, float *** c1313h,  float *** c1213ha, 
float *** c1112h, float *** c1222h, float *** c1233h, float *** c1223h, float *** c1213h,  float *** c1212h, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk);

void FD_TRICL();

void info(FILE *fp);

void model(float *** rho, 
float *** c1111, float *** c1112, float *** c1113, float *** c1122, float *** c1133, float *** c1123, float *** c1212, 
float *** c1213, float *** c1222, float *** c1223, float *** c1233, float *** c1313, float *** c1322, float *** c1323, 
float *** c1333, float *** c2222, float *** c2223, float *** c2233, float *** c2323, float *** c2333, float *** c3333,
float *** taus, float *** taup, float * eta, float * x, float * y, float * z);

void readmod(float *** rho, 
float *** c1111, float *** c1112, float *** c1113, float *** c1122, float *** c1133, float *** c1123, float *** c1212, 
float *** c1213, float *** c1222, float *** c1223, float *** c1233, float *** c1313, float *** c1322, float *** c1323, 
float *** c1333, float *** c2222, float *** c2223, float *** c2233, float *** c2323, float *** c2333, float *** c3333,
float *** taus, float *** taup, float * eta);

void timing(double * time_avg,  double * time_std, int ishot);

double update_s_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct tensor3d ***e, struct tensor3d ***t,
float *** c1111,  float *** c1122,  float *** c1133,  float *** c1123,  float *** c1113,   float *** c1112, 
                  float *** c2222,  float *** c2233,  float *** c2223,  float *** c1322,   float *** c1222, 
                                    float *** c3333,  float *** c2333,  float *** c1333,   float *** c1233, 
float *** c1123h, float *** c2223h, float *** c2333h, float *** c2323h, float *** c1323ha, float *** c1223ha, 
float *** c1113h, float *** c1322h, float *** c1333h, float *** c1323h, float *** c1313h,  float *** c1213ha, 
float *** c1112h, float *** c1222h, float *** c1233h, float *** c1223h, float *** c1213h,  float *** c1212h, 
float *** absorb_coeff,
float * wxp1, float * wx1, float * wxp2, float * wx2, 
float * wyp1, float * wy1, float * wyp2, float * wy2, 
float * wzp1, float * wz1, float * wzp2, float * wz2,
int infoout);

void write_par(FILE *fp);
