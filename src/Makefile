# Makefile for fd3d

#--------------------------------------------------------
# edit here:

# source code for model generation

MODELPS = model_ps.c
MODELP  = model_p.c
MODEL_ORTHO = model_ortho.c
MODEL_TRICL = model_tricl.c
EXEC= ../bin

# On Linux NEC-cluster with Intel-MPI
# CC=mpiicc
# LFLAGS=-lm -lcseife -lfftw3 -ftz -lstdc++
# CFLAGS=-O3 -xAVX -ipo -fno-fnalias -restrict
# SFLAGS=-L./../libcseife -L$//sfs/fs2/work-sh2/sungw331/fftw-3.3.4/lib
# IFLAGS=-I./../libcseife -I$//sfs/fs2/work-sh2/sungw331/fftw-3.3.4/include -I./../include

# On Linux NEC-cluster with Intel-MPI
#CC=mpiicc
#LFLAGS=-lm -ftz -lstdc++
#CFLAGS=-O3 -xAVX -ipo -fno-fnalias -restrict
#IFLAGS=-I./../include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
CC=mpicc
LFLAGS=-lm -lstdc++
CFLAGS=-O3 -w -fno-stack-protector -D_FORTIFY_SOURCE=0
IFLAGS=-I./../include

# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) -c $< $(IFLAGS)

VPATH = MODELS:TRICL:ORTHO:ISO:AC

FD_UTIL = \
	catseis.c \
	copydsk_array.c \
	note.c \
	readdsk.c \
	read_grid.c \
	round_to_grid.c \
	round_to_grid_min.c \
	round_to_grid_max.c \
	update_timing.c \
	util.c \
	writedsk.c \
	writedsk_array.c


SNAPMERGE_UTIL = \
	merge.c \
	read_grid_glob.c \
	read_par.c \
	$(FD_UTIL)

SNAPMERGE_SCR = \
	snapmerge_el.c \
	$(SNAPMERGE_UTIL)
	
SNAPMERGEP_SCR = \
	snapmerge_p.c \
	$(SNAPMERGE_UTIL)


FD3D_UTIL = \
	absorb.c \
	av_mat.c \
	check_par.c \
	exchange_par.c \
	fsource.c \
	initproc_cart.c \
	matcopy.c \
	merge.c \
	mergemod.c \
	output_source_signal.c \
	outseis.c \
	pml_profile.c \
	read_par.c \
	read_wavelet.c \
	receiver.c \
	shotno_glob2loc.c \
	splitrec.c \
	splitsrc.c \
	wavelet.c \
	writemod.c \
	$(FD_UTIL)

FD3D_ISO_SRC= \
	timing_iso.c

FD3D_EL_SRC= \
	exchange_s_2nd_el.c \
	exchange_v_2nd_el.c \
	pml_update_v_el.c \
	psource_el.c \
	saveseis_el.c \
	savesnap_el.c \
	snap_el.c \
	seismo_el.c \
	sources_el.c \
	update_v_2nd_el.c

FD3D_ANISO_SRC= \
	eigenvalue_analytic.c \
	pml_update_e_aniso.c \
	update_e_el_2nd_aniso.c

SAVA_ISO_SRC= \
	FD_ISO.c \
	forward_ISO.c \
	forward_shot_ISO.c \
	checkfd_ps.c \
	fd_coeff_ps.c \
	info_ps.c \
	pml_update_s_el_ps.c \
	readmod_ps.c \
	reset_wavefields_iso.c \
	sava_ISO.c \
	update_s_el_2nd_ps.c \
	update_s_ve_2nd_ps.c \
	write_par_ps.c \
	$(MODELPS) \
	$(FD3D_ISO_SRC) \
	$(FD3D_EL_SRC) \
	$(FD3D_UTIL)
	 
SAVA_AC_SRC= \
	FD_AC.c \
	checkfd_p.c \
	exchange_s_2nd_p.c \
	exchange_v_2nd_p.c \
	fd_coeff_p.c \
	forward_shot_AC.c \
	info_p.c \
	pml_update_s_ac_p.c \
	pml_update_v_p.c \
	psource_p.c \
	readmod_p.c \
	sava_AC.c \
	saveseis_p.c \
	savesnap_p.c \
	seismo_p.c \
	snap_p.c \
	sources_p.c \
	update_s_ac_2nd_p.c \
	update_s_va_2nd_p.c \
	update_v_2nd_p.c \
	write_par_p.c \
	$(MODELP) \
	$(FD3D_ISO_SRC) \
	$(FD3D_UTIL)
	
SAVA_ORTHO_SRC= \
	sava_ORTHO.c \
	checkfd_ortho.c \
	fd_coeff_ortho.c \
	FD_ORTHO.c \
	forward_ORTHO.c \
	forward_shot_ORTHO.c \
	info_ortho.c \
	readmod_ortho.c \
	reset_wavefields_ortho.c \
	timing_ortho.c \
	update_s_el_2nd_ortho.c \
	write_par_ortho.c \
	$(MODEL_ORTHO) \
	$(FD3D_EL_SRC) \
	$(FD3D_ANISO_SRC) \
	$(FD3D_UTIL)
	
SAVA_TRICL_SRC= \
	sava_TRICL.c \
	checkfd_tricl.c \
	exchange_e_2nd_tricl.c \
	fd_coeff_tricl.c \
	FD_TRICL.c \
	forward_TRICL.c \
	forward_shot_TRICL.c \
	info_tricl.c \
	readmod_tricl.c \
	reset_wavefields_tricl.c \
	timing_tricl.c \
	update_s_el_2nd_tricl.c \
	write_par_tricl.c \
	$(MODEL_TRICL) \
	$(FD3D_EL_SRC) \
	$(FD3D_ANISO_SRC) \
	$(FD3D_UTIL)
	
SNAPMERGE_OBJ = $(SNAPMERGE_SCR:%.c=%.o)

SNAPMERGEP_OBJ = $(SNAPMERGEP_SCR:%.c=%.o)

SAVA_ISO_OBJ = $(SAVA_ISO_SRC:%.c=%.o)

SAVA_AC_OBJ = $(SAVA_AC_SRC:%.c=%.o)

SAVA_ORTHO_OBJ = $(SAVA_ORTHO_SRC:%.c=%.o)

SAVA_TRICL_OBJ = $(SAVA_TRICL_SRC:%.c=%.o)

snapmerge:		$(SNAPMERGE_OBJ)
	$(CC) $(SNAPMERGE_OBJ) -o $(EXEC)/snapmerge.exe $(LFLAGS)
	
snapmergep:		$(SNAPMERGEP_OBJ)
	$(CC) $(SNAPMERGEP_OBJ) -o $(EXEC)/snapmergep.exe $(LFLAGS)

sava_ISO:		$(SAVA_ISO_OBJ)
	$(CC) $(SAVA_ISO_OBJ) -o $(EXEC)/sava_ISO $(LFLAGS)
	
sava_AC:		$(SAVA_AC_OBJ)
	$(CC) $(SAVA_AC_OBJ) -o $(EXEC)/sava_AC $(LFLAGS)

fd3daniso:		$(FD3D_TRICL_OBJ)
	$(CC) $(FD3D_TRICL_OBJ) -o $(EXEC)/fd3daniso.exe $(LFLAGS)
	
sava_ORTHO:		$(SAVA_ORTHO_OBJ)
	$(CC) $(SAVA_ORTHO_OBJ) -o $(EXEC)/sava_ORTHO $(LFLAGS)
	
sava_TRICL:		$(SAVA_TRICL_OBJ)
	$(CC) $(SAVA_TRICL_OBJ) -o $(EXEC)/sava_TRICL $(LFLAGS)

clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
