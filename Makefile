CC=gcc
HDF5_INCLUDES=-I${HDF5_DIR}/include
HDF5_LIBRARIES=-L${HDF5_DIR}/lib -lhdf5_hl -lhdf5
CVODE_INCLUDES=-I${HDF5_DIR}/include/cvode
CVODE_LIBRARIES=-L${HDF5_DIR}/lib -lsundials_cvode -lsundials_nvecserial
CCFLAGS=-O0

all : test_primordial_solver run_test

clean :
	rm -f cvode_templates/primordial_cvode_solver.c

cvode_templates/primordial_cvode_solver.c : \
        cvode_templates/cvode_solver.c.template \
        cvode_templates/primordial_cvode_solver.c.template \
        dengo/*.py
	@echo "Regenerating solver from template"
	@python2.6 dengo/write_cvode_solver.py

test_primordial_solver : cvode_templates/primordial_cvode_solver.c 
	@echo "Recompiling"
	$(CC) $(CCFLAGS) \
        -Icvode_templates \
        $(HDF5_INCLUDES)   \
        $(CVODE_INCLUDES)  \
        $(HDF5_LIBRARIES)  \
        $(CVODE_LIBRARIES) \
        cvode_templates/primordial_cvode_solver.c \
        -DSAMPLE_TEST_PROBLEM -g \
		-o test_primordial_solver

run_test : test_primordial_solver
	./test_primordial_solver
	python2.6 plot_cvode_output.py
