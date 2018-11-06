# Makefile for running dengo 
# again it is generated from template



# print cvode statistics, i.e. number of evaluations
#OPTIONS+= -DPRINT_CVODE_STATS 


# dengo currently supports 3 different solver modes from sundials
# CVSPILS, CVDls, CVKLU
# choose one from below, default: CVDLS
# MAX_NCELLS: number of strip solved per solver call 

# Use the scaled, preconditioned, iter- ative linear solver (CVSpils) interfaces
#OPTIONS = -DCVSPILS -DMAX_NCELLS=1024

# Use the Dense Linear Solver interface (CVDLS)
#OPTIONS+= -DMAX_NCELLS=1

# Use the Sparse KLU direct solver
OPTIONS = -DCVKLU -DMAX_NCELLS=256

# Use OPENMP?
OPTIONS += -fopenmp -DNTHREADS=8


# INSTALL PATHS
CVODE_PATH    = /home/kwoksun2/cvode-3.1.0/instdir
SUITESPARSE_PATH = /home/kwoksun2/SuiteSparse
DENGO_INSTALL_PATH = /home/kwoksun2/dengo_install

CC      = gcc
LIBS     = -lhdf5 -lhdf5_hl -lm -lstdc++
INCLUDE  = -I. -I$(CVODE_PATH)/include


# REQUIRED FILES
OBJ = cvklu_solver.C initialize_cvode_solver.C
TEST_OBJ = cvklu_solver_main.C

SHELL = sh

exec_prefix  = $(CVODE_PATH)/instdir
includedir   = $(CVODE_PATH)/include
libdir       = $(CVODE_PATH)/lib

CFLAGS   = -O3 -DNDEBUG -g
LDFLAGS  = 
LIBS    +=  /usr/lib/x86_64-linux-gnu/librt.so

LINKFLAGS = -Wl,-rpath,$(libdir)


TMP_INCS  = . ${includedir} 
INCLUDES  = $(addprefix -I, ${TMP_INCS})

#------------------------------------------------------------------------
# KLU part
#

# Taken from CVODE example makefile
# LOOKS LIKE I NEED TO INCLUDE CVODE PATH IN LD_LIBRARY_PATH for it to work

INCLUDES_KLU  = $(SUITESPARSE_PATH)/include
LIBRARIES_KLU = -lsundials_sunlinsolklu   $(SUITESPARSE_PATH)/KLU/Lib/libklu.so $(SUITESPARSE_PATH)/AMD/Lib/libamd.so $(SUITESPARSE_PATH)/COLAMD/Lib/libcolamd.so $(SUITESPARSE_PATH)/BTF/Lib/libbtf.so $(SUITESPARSE_PATH)/SuiteSparse_config/libsuitesparseconfig.so
LINKFLAGS_KLU = ${LINKFLAGS}:$(SUITESPARSE_PATH)/KLU/Lib
LIBRARIES = -lsundials_cvode -lsundials_nvecserial ${LIBS}


ifneq (,$(findstring -DCVKLU, $(OPTIONS)))
LIBRARIES += $(LIBRARIES_KLU)
INCLUDES  += $(addprefix -I, $(INCLUDES_KLU) )
LINKFLAGS := $(LINKFLAGS_KLU)
endif
#------------------------------------------------------------------------




#------------------------------------------------------------------------
# MAKE LIBDENGO BY DEFAULT
# -----------------------------------------------------------------------

libdengo.la: $(addsuffix .lo, $(basename $(OBJ)) )
	mkdir -p $(DENGO_INSTALL_PATH)
	mkdir -p $(DENGO_INSTALL_PATH)/lib
	mkdir -p $(DENGO_INSTALL_PATH)/include
	libtool --mode=link --tag=CC $(CC) $(LDFLAGS) -o $@ $^ $(OPTIONS)  ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} $(LINKFLAGS) -rpath $(DENGO_INSTALL_PATH)/lib
	libtool --mode=install install -c libdengo.la $(DENGO_INSTALL_PATH)/lib
	cp cvklu_solver.h $(DENGO_INSTALL_PATH)/include/
	cp cvklu_tables.h5 $(DENGO_INSTALL_PATH)

run_dengo: cvklu_solver_main.C
	$(CC) -o $@ $^ $(OPTIONS) -I$(DENGO_INSTALL_PATH)/include -I$(CVODE_PATH)/include -I$(INCLUDES_KLU) -L$(DENGO_INSTALL_PATH)/lib -ldengo -lm

test: run_dengo 
	./run_dengo



#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------


# Inhibit removing any *.o files after compiling

.PRECIOUS: %.o %.lo %.C

%.lo: %.C
	  libtool --mode=compile --tag=CC $(CC) -c $(OPTIONS) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $*.C

clean:
	rm -f *.o
	rm -f *.so
	rm -f run_dengo
	rm -f libdengo.la
	rm -f *.lo
	libtool --mode=uninstall rm -r $(DENGO_INSTALL_PATH)/lib/libdengo.la


