/*------------------------------------------------------------------------
 *  fdaniso.h - include file for program FD3D (elastic, anisotropic)
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------*/

/* declaration of functions */

//void eigenvalue(float ** s, float * e, int n);
void eigenvalue(float * s, float * e);

double pml_update_e(struct vector3d ***v, struct tensor3d ***e, 
struct pml *pmlle, struct pml *pmlri, 
struct pml *pmlba, struct pml *pmlfr, 
struct pml *pmlto, struct pml *pmlbo, 
struct vector3d *** Pv_x_l, struct vector3d *** Pv_x_r,
struct vector3d *** Pv_y_b, struct vector3d *** Pv_y_f,
struct vector3d *** Pv_z_t, struct vector3d *** Pv_z_b,
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, int infoout);

double update_e_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, struct tensor3d ***e, struct divcurl3d ***w,
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, 
int infoout);

