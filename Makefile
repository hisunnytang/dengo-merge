CC=gcc
HDF5_INCLUDES=-I${HDF5_DIR}/include
HDF5_LIBRARIES=-L${HDF5_DIR}/lib -lhdf5_hl -lhdf5
CVODE_INCLUDES=-I${HDF5_DIR}/include/cvode
CVODE_LIBRARIES=-L${HDF5_DIR}/lib -lsundials_cvode -lsundials_nvecserial
CCFLAGS=-O3

all : test_primordial_solver run_test

clean :
	rm -f simple_cvode_solver/primordial_cvode_solver.c

simple_cvode_solver/primordial_cvode_solver.c : \
        simple_cvode_solver/cvode_solver.c.template \
        simple_cvode_solver/primordial_cvode_solver.c.template \
        primordial_chemistry/*.py
	@echo "Regenerating solver from template"
	@python2.6 primordial_chemistry/write_cvode_solver.py

test_primordial_solver : simple_cvode_solver/primordial_cvode_solver.c 
	@echo "Recompiling"
	$(CC) $(CCFLAGS) \
        -Isimple_cvode_solver \
        $(HDF5_INCLUDES)   \
        $(CVODE_INCLUDES)  \
        $(HDF5_LIBRARIES)  \
        $(CVODE_LIBRARIES) \
        simple_cvode_solver/primordial_cvode_solver.c \
        -DSAMPLE_TEST_PROBLEM -g \
		-o test_primordial_solver

run_test : test_primordial_solver
	./test_primordial_solver
	python2.6 plot_cvode_output.py
