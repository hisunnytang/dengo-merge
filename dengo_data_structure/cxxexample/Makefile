CVODE_INSTDIR= /home/kwoksun2/cvode-3.1.0/instdir



CXX = g++
CXXFLAGS=-I. -I$(CVODE_INSTDIR)/include

LDIR = -L$(CVODE_INSTDIR)/lib -lm -lhdf5 -lhdf5_hl -lstdc++
EXTRA_OBJ = $(CVODE_INSTDIR)/lib/libsundials_cvode.a $(CVODE_INSTDIR)/lib/libsundials_nvecserial.a

OBJ = cvdls_9species_solver.o cvdls_9species_solver_main.o cvodes_solver_CVDls.o

DEPS = cvdls_9species_solver.h

%.C:
	$(CXX)  -c $< $(CXXFLAGS) 

hellomake: $(OBJ)
	$(CXX) -o hellomake $(OBJ) $(EXTRA_OBJ) $(CXXFLAGS) $(LDIR)
