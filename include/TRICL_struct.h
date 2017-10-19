/*------------------------------------------------------------------------
 *  TRICL_struct.h - define data structures for TRICL forward problem 
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
   struct tensor3d    ***t, ***e, ***r;
   struct vector3d    ***v, ***a;
   struct divcurl3d   ***w;
} wave;

/* static variables (material parameters and axes) */
struct mat{
   float   ***rho;
   float   ***c1111, ***c1112, ***c1113, ***c1122, ***c1133, ***c1123, ***c1212;
   float   ***c1213, ***c1222, ***c1223, ***c1233, ***c1313, ***c1322, ***c1323;
   float   ***c1333, ***c2222, ***c2223, ***c2233, ***c2323, ***c2333, ***c3333;
   float   ***taus, ***taup, *eta;
   float   ***rhoijpkp, ***rhoipjkp, ***rhoipjpk;
   float   ***c1123h, ***c2223h, ***c2333h, ***c2323h, ***c1323ha, ***c1223ha;
   float   ***c1113h, ***c1322h, ***c1333h, ***c1323h, ***c1313h,  ***c1213ha;
   float   ***c1112h, ***c1222h, ***c1233h, ***c1223h, ***c1213h,  ***c1212h;
   float   ***tausipjk, ***tausijpk, ***tausijkp;
} mat;

/* grid geometry */
struct geom{
   float   *x, *xp, *y, *yp, *z, *zp, *xg, *xpg, *yg, *ypg, *zg, *zpg;
   float   *dxp, *dx, *dyp, *dy, *dzp, *dz;
   float   *wxp1, *wx1, *wxp2, *wx2;
   float   *wyp1, *wy1, *wyp2, *wy2;
   float   *wzp1, *wz1, *wzp2, *wz2;
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
   struct pml          *pmlle, *pmlri, *pmlba, *pmlfr, *pmlto, *pmlbo;
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

   float  ***sbuf_s_lef_to_rig, **sbuf_s_rig_to_lef, 
	  ***sbuf_s_bac_to_fro, **sbuf_s_fro_to_bac, 
	  ***sbuf_s_top_to_bot, **sbuf_s_bot_to_top, 
	  ***rbuf_s_lef_to_rig, **rbuf_s_rig_to_lef, 
	  ***rbuf_s_bac_to_fro, **rbuf_s_fro_to_bac, 
	  ***rbuf_s_top_to_bot, **rbuf_s_bot_to_top;
     
   float  ***sbuf_e_lef_to_rig, ***sbuf_e_rig_to_lef, 
	  ***sbuf_e_bac_to_fro, ***sbuf_e_fro_to_bac, 
	  ***sbuf_e_top_to_bot, ***sbuf_e_bot_to_top, 
	  ***rbuf_e_lef_to_rig, ***rbuf_e_rig_to_lef, 
	  ***rbuf_e_bac_to_fro, ***rbuf_e_fro_to_bac, 
	  ***rbuf_e_top_to_bot, ***rbuf_e_bot_to_top; 

   MPI_Request   request_e[12], request_v[12], request_s[12];

} mpi;

