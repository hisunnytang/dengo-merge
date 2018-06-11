# paths
CVODE_DIR = /home/kwoksun2/cvode-3.1.0/instdir


# declare variables
CXX = g++ 
LIBS = -lhdf5 -lhdf5_hl -lm -lstdc++
INCLUDE = -I. -I$(CVODE_DIR)/include
CXXFLAGS =  -fPIC -g -fopenmp
OBJ = _temp_omp_solver.o initialize_cvode_solver.o _temp_omp_solver_main.o

# we should make a big header files 
# that has everything declared
DEPS = dengo_types.h

%.o: %.C
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE) 

program: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CVODE_DIR)/lib/libsundials_cvode.a $(CVODE_DIR)/lib/libsundials_nvecserial.a $(LIBS) $(INCLUDE)

libdengo.so: $(OBJ) 
	$(CXX) -shared -o $@ $^ $(CVODE_DIR)/lib/libsundials_cvode.a $(CVODE_DIR)/lib/libsundials_nvecserial.a $(LIBS) $(INCLUDE)

dengo_mem:
	$(CXX) $(CXXFLAGS) -c cvdls_9species_solver_mem.C $(INCLUDE)
	$(CXX) $(CXXFLAGS) -c initialize_cvode_solver.C $(INCLUDE)
	$(CXX) $(CXXFLAGS) -o dengo_mem dengo_data_struct_main.C cvdls_9species_solver_mem.o initialize_cvode_solver.o $(CVODE_DIR)/lib/libsundials_cvode.a $(CVODE_DIR)/lib/libsundials_nvecserial.a  $(LIBS) $(INCLUDE)


clean:
	rm -f *.o
	rm -f *.so
	rm dengo_struct



