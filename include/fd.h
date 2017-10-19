/*------------------------------------------------------------------------
 *  fd.h - general include file for all FD-programs
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    
/*#define abs(x) ((x<0.0)?(-x):x)*/

#define PI (3.141592653589793)
#define NPAR 64
#define STRING_SIZE 128
#define REQUEST_COUNT 4


/* definition of structures */
struct pml{
	float	exp;
	float	x1;
	float	x2;
	float	expp;
	float	xp1;
	float	xp2;
};
struct pmlp{
	float	exp;
	float	x1;
	float	x2;
};

struct vector2d{
	float	x;
	float	z;
};
struct tensor2dpsv{
	float	xx;
	float	yy;
	float	zz;
	float	xz;
};
struct tensor2dsh{
	float	yz;
	float	xy;
};
struct divcurl2d{
	float	div;
	float	curly;
};
struct curl2d{
	float	curlx;
	float	curlz;
};

struct vector3d{
	float	x;
	float	y;
	float	z;
};
struct tensor3d{
	float	xx;
	float	yy;
	float	zz;
	float	yz;
	float	xz;
	float	xy;
};
struct divcurl3d{
	float	div;
	float	curlx;
	float	curly;
	float	curlz;
};

struct vectorbh{
	float	r;
	float	t;
	float	z;
};
struct tensorbh{
	float	rr;
	float	tt;
	float	zz;
	float	tz;
	float	rz;
	float	rt;
};
struct divcurlbh{
	float	div;
	float	curlr;
	float	curlt;
	float	curlz;
};

/* acquisition geometry */
struct acq{
   float   **signals;
   float   **srcpos, **srcpos_loc;
   int     **recpos, **recpos_loc;
   int     *recswitch;
} acq;

/* declaration of functions */
void	catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns);
void	copydsk_array(FILE *fp_out, FILE *fp_in, int na, int format);
void	note(FILE *fp);
float	readdsk(FILE *fp_in, int format);
void	read_grid(FILE *fp, char *dx_file, int dim, float *x, float *xp, float *xg, float *xpg);
int	*irnd_to_grid(float pos, float *xg, float *xpg, int nxg);
int	irnd_to_grid_min(float pos, float *xp, int nx1, int nx2, int dx);
int	irnd_to_grid_max(float pos, float *xp, int nx1, int nx2, int dx);
void	update_timing(int nt, double * time_update, double * time_sum, double * time_avg, double * time_std, int ntime);
void	writedsk(FILE *fp_out, float amp, int format);
void	writedsk_array(FILE *fp_out, int na, float *amp, int format);


/* utility functions */
void	error(char err_text[]);
void	warning(char warn_text[]);
double	maximum(float **a, int nx, int ny);

float	 *vector(int nl, int nh);
int	*ivector(int nl, int nh);
double	*dvector(int nl, int nh);
float	 **matrix(int nrl, int nrh, int ncl, int nch);
int	**imatrix(int nrl, int nrh, int ncl, int nch);
float	***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
struct pml	*pml_vector(int nl, int nh);
struct pmlp	*pmlp_vector(int nl, int nh);
struct vector2d		   **vector2d_matrix(int nrl, int nrh, int ncl, int nch);
struct tensor2dpsv	**tensor2dpsv_matrix(int nrl, int nrh, int ncl, int nch);
struct tensor2dsh	 **tensor2dsh_matrix(int nrl, int nrh, int ncl, int nch);
struct divcurl2d	  **divcurl2d_matrix(int nrl, int nrh, int ncl, int nch);
struct curl2d	             **curl2d_matrix(int nrl, int nrh, int ncl, int nch);
struct vector3d		 ***vector3d_tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
struct tensor3d		 ***tensor3d_tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
struct divcurl3d	***divcurl3d_tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
struct vectorbh		 **vectorbh_matrix(int nrl, int nrh, int ncl, int nch);
struct tensorbh		 **tensorbh_matrix(int nrl, int nrh, int ncl, int nch);
struct divcurlbh	**divcurlbh_matrix(int nrl, int nrh, int ncl, int nch);

void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
void free_pml_vector( struct pml  *v, int nl, int nh);
void free_pmlp_vector(struct pmlp *v, int nl, int nh);
void free_vector2d_matrix(   struct vector2d    **m, int nrl, int nrh, int ncl, int nch);
void free_tensor2dpsv_matrix(struct tensor2dpsv **m, int nrl, int nrh, int ncl, int nch);
void free_tensor2dsh_matrix( struct tensor2dsh  **m, int nrl, int nrh, int ncl, int nch);
void free_divcurl2d_matrix(  struct divcurl2d   **m, int nrl, int nrh, int ncl, int nch);
void free_curl2d_matrix(     struct curl2d      **m, int nrl, int nrh, int ncl, int nch);
void free_vector3d_tensor(struct  vector3d  ***t, int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_tensor3d_tensor(struct  tensor3d  ***t, int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_divcurl3d_tensor(struct divcurl3d ***t, int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vectorbh_matrix( struct vectorbh  **m, int nrl, int nrh, int ncl, int nch);
void free_tensorbh_matrix( struct tensorbh  **m, int nrl, int nrh, int ncl, int nch);
void free_divcurlbh_matrix(struct divcurlbh **m, int nrl, int nrh, int ncl, int nch);
