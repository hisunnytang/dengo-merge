# Makefile for running dengo 
# again it is generated from template



# print cvode statistics, i.e. number of evaluations
#OPTIONS+= -DPRINT_CVODE_STATS 


# dengo currently supports 3 different solver modes from sundials
# CVSPILS, CVDls, CVKLU
# choose one from below, default: CVDLS
# MAX_NCELLS: number of strip solved per solver call 

# INSTALL PATHS
DENGO_INSTALL_PATH = /home/kwoksun2/dengo_install
HDF5_PATH = /home/kwoksun2/anaconda3


LIBTOOL = /usr/bin/libtool


CC      = gcc
LIBS    = -lhdf5 -lhdf5_hl -lm -lstdc++

# REQUIRED FILES
OBJ = predator_prey_solver.C BE_chem_solve.C
TEST_OBJ = predator_prey_solver_main.C

includedir   = $(HDF5_PATH)/include
libdir       = $(HDF5_PATH)/lib

CFLAGS   = -O3 -DNDEBUG -g -fopenmp
LDFLAGS  = 
LIBS    +=  /usr/lib/x86_64-linux-gnu/librt.so

LINKFLAGS = -Wl,-rpath,$(HDF5_PATH)/lib

TMP_INCS  = . ${includedir} 
INCLUDES  = $(addprefix -I, ${TMP_INCS})
LIBDIR    = $(addprefix -L, ${libdir})
LIBRARIES = ${LIBS}

#------------------------------------------------------------------------
# MAKE LIBDENGO BY DEFAULT
# -----------------------------------------------------------------------

libdengo.la: $(addsuffix .lo, $(basename $(OBJ)) )
	mkdir -p $(DENGO_INSTALL_PATH)
	mkdir -p $(DENGO_INSTALL_PATH)/lib
	mkdir -p $(DENGO_INSTALL_PATH)/include
	$(LIBTOOL) --mode=link --tag=CC $(CC) $(LDFLAGS) -o $@ $^ $(OPTIONS)  ${CFLAGS} ${LDFLAGS} ${INCLUDES} ${LIBDIR} ${LIBRARIES} $(LINKFLAGS) -rpath $(DENGO_INSTALL_PATH)/lib
	$(LIBTOOL) --mode=install install -c libdengo.la $(DENGO_INSTALL_PATH)/lib
	cp predator_prey_solver.h $(DENGO_INSTALL_PATH)/include/
	cp predator_prey_tables.h5 $(DENGO_INSTALL_PATH)

run_dengo: predator_prey_solver_main.C
	$(CC) -o $@ $^ $(OPTIONS) -I$(DENGO_INSTALL_PATH)/include -I$(HDF5_PATH)/include -L$(HDF5_PATH)/lib -lhdf5_hl -lhdf5 -L$(DENGO_INSTALL_PATH)/lib -ldengo -lm

test: run_dengo 
	./run_dengo



#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------


# Inhibit removing any *.o files after compiling

.PRECIOUS: %.o %.lo %.C

%.lo: %.C
	  $(LIBTOOL) --mode=compile --tag=CC $(CC) -c $(OPTIONS) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $*.C

clean:
	rm -f *.o
	rm -f *.so
	rm -f run_dengo
	rm -f libdengo.la
	rm -f *.lo
	libtool --mode=uninstall rm -r $(DENGO_INSTALL_PATH)/lib/libdengo.la


