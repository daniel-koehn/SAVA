/*------------------------------------------------------------------------
 *  AC_struct.h - define data structures for AC forward problem 
 *
 *  D. Koehn
 *  Kiel, 15th of October 2017 
 *  ---------------------------------------------------------------------*/

/* timing variables */
struct times{
   double   time1, time2, time3, time4;
   double   *time_update, *time_sum, *time_avg, *time_std;
} times;

/* wavefield variables */
struct wave{
   float	     ***p, ***rpp;
   struct vector3d   ***v, ***a;
   float	     ***diverg;
} wave;

/* material parameters */
struct mat{
   float    ***pi, ***rho;
   float    ***taup, *eta;
   float    ***rhoijpkp, ***rhoipjkp, ***rhoipjpk;
} mat;

/* grid geometry */
struct geom{
   float    *x, *xp, *y, *yp, *z, *zp, *xg, *xpg, *yg, *ypg, *zg, *zpg;
   float    *dx, *dxp, *dy, *dyp, *dz, *dzp;
} geom;

/* seismogram variables */
struct seis{
   float   **sectionvx,  **sectionvy, **sectionvz; 
   float   **sectionax,  **sectionay, **sectionaz; 
   float   **sectiondiv, **sectionp;
   float   **section_fulldata;
} seis;

/* PML variables */
struct pmls{
   struct pml   *pmlle, *pmlri, *pmlba, *pmlfr, *pmlto, *pmlbo;
   float	***Pp_x_l, ***Pp_x_r, ***Pp_y_b, ***Pp_y_f, ***Pp_z_t, ***Pp_z_b;
   float	***Pvx_x_l, ***Pvx_x_r, ***Pvy_y_b, ***Pvy_y_f, ***Pvz_z_t, ***Pvz_z_b;
   float	***absorb_coeff;
} pmls;

/* MPI variables and buffer for exchange between MPI processes */
struct mpi{
   float   **sbuf_lef_to_rig, **sbuf_rig_to_lef,
	   **sbuf_bac_to_fro, **sbuf_fro_to_bac, 
	   **sbuf_top_to_bot, **sbuf_bot_to_top,
	   **rbuf_lef_to_rig, **rbuf_rig_to_lef, 
	   **rbuf_bac_to_fro, **rbuf_fro_to_bac, 
	   **rbuf_top_to_bot, **rbuf_bot_to_top;
   MPI_Request   request_v[6], request_s[6];
} mpi;


