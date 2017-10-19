/*------------------------------------------------------------------------
 *  ISO_struct.h - define data structures for ISO forward problem 
 *
 *  D. Koehn
 *  Kiel, 15th of October 2017 
 *  ---------------------------------------------------------------------*/

/* timing variables */
struct times{
   double   time1, time2, time3, time4;
   double   *time_update, *time_sum, *time_avg, *time_std;
} times;

/* dynamic variables (wavefield) */
struct wave{
   struct tensor3d    ***t, ***r;
   struct vector3d    ***v, ***a;
   struct divcurl3d   ***w;
} wave;

/* static variables (material parameters and axes) */
struct mat{
float	***pi, ***u, ***rho;
float	***taus, ***taup, *eta;
float	***uipjk, ***uijpk, ***uijkp, 
	***tausipjk, ***tausijpk, ***tausijkp, 
	***rhoijpkp, ***rhoipjkp, ***rhoipjpk;
} mat;

/* grid geometry */
struct geom{
   float    *x, *xp, *y, *yp, *z, *zp, *xg, *xpg, *yg, *ypg, *zg, *zpg;
   float    *dx, *dxp, *dy, *dyp, *dz, *dzp;
} geom;

/* seismogram variables */
struct seis{
   float   **sectiontxx, **sectiontxy, **sectiontxz, **sectiontyy, **sectiontyz,  **sectiontzz; 
   float   **sectionvx,  **sectionvy,  **sectionvz; 
   float   **sectionax,  **sectionay,  **sectionaz; 
   float   **sectioncurlx, **sectioncurly, **sectioncurlz, **sectiondiv, **sectionp;
   float   **section_fulldata;
} seis;

/* PML variables */
struct pmls{
   struct pml           *pmlle, *pmlri, *pmlba, *pmlfr, *pmlto, *pmlbo;
   struct vector3d   ***Pt_x_l, ***Pt_x_r, ***Pt_y_b, ***Pt_y_f, ***Pt_z_t, ***Pt_z_b;
   struct vector3d   ***Pv_x_l, ***Pv_x_r, ***Pv_y_b, ***Pv_y_f, ***Pv_z_t, ***Pv_z_b;
   float	     ***absorb_coeff;
} pmls;

/* MPI variables and buffer for exchange between MPI processes */
struct mpi{
   float   **sbuf_v_lef_to_rig, ***sbuf_v_rig_to_lef, 
	   **sbuf_v_bac_to_fro, ***sbuf_v_fro_to_bac, 
	   **sbuf_v_top_to_bot, ***sbuf_v_bot_to_top, 
	   **rbuf_v_lef_to_rig, ***rbuf_v_rig_to_lef, 
	   **rbuf_v_bac_to_fro, ***rbuf_v_fro_to_bac, 
	   **rbuf_v_top_to_bot, ***rbuf_v_bot_to_top;

   float   ***sbuf_s_lef_to_rig, **sbuf_s_rig_to_lef, 
	   ***sbuf_s_bac_to_fro, **sbuf_s_fro_to_bac, 
	   ***sbuf_s_top_to_bot, **sbuf_s_bot_to_top, 
	   ***rbuf_s_lef_to_rig, **rbuf_s_rig_to_lef, 
	   ***rbuf_s_bac_to_fro, **rbuf_s_fro_to_bac, 
	   ***rbuf_s_top_to_bot, **rbuf_s_bot_to_top;

   MPI_Request   request_v[12], request_s[12];
} mpi;

