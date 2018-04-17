# paths
CVODE_DIR = /home/kwoksun2/cvode-3.1.0/instdir


# declare variables
CXX = g++
LIBS = -lhdf5 -lhdf5_hl -lm -lstdc++
INCLUDE = -I. -I$(CVODE_DIR)/include
CXXFLAGS =  -fPIC
OBJ = cvdls_9species_solver.o cvodes_solver_CVDls.o

# we should make a big header files 
# that has everything declared
DEPS = dengo_types.h

%.o: %.C
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE) 

libdengo.so: $(OBJ) 
	$(CXX) -shared -o $@ $^ $(CVODE_DIR)/lib/libsundials_cvode.a $(CVODE_DIR)/lib/libsundials_nvecserial.a $(LIBS) $(INCLUDE)


clean:
	rm *.o
	rm *.so

