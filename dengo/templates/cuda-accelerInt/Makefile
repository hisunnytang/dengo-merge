CC=g++
NVCC=nvcc
CXXFLAGS= -fopenmp -Wextra -O3 -std=c++11 -g
CUDAFLAGS= -m64 -Xptxas -v -O3 -std=c++11 -Xcompiler -fopenmp --ftz=false --prec-div=true --prec-sqrt=true --fmad=false -g -lineinfo -dc
LIBS = -lstdc++ -lcuda -lcudart -L/opt/intel/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64 -lm -lmkl_rt -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lmkl_mc -lmkl_def -lhdf5 -lhdf5_hl

LIBDIRS = -L/usr/local/cuda/samples/common/lib/linux/x86_64//linux -L/usr/local/cuda/lib64
INCDIRS =-I/usr/local/cuda/samples/common/inc -I/usr/local/cuda/include
INCDIRS+=-I/home/kwoksun2/sundials-3.2.1/instdir/include


OPTIONS = -DRADAU2A -DUSE_DRRATE -DUSE_DCRATE

%.o: %.cu
	$(NVCC) $(INCDIRS) $(CUDAFLAGS) -c $< -o $@ $(LIBDIRS) $(LIBS) $(OPTIONS)

%.o: %.C
	$(CC) $(INCDIRS) $(CXXFLAGS) -c $< -o $@ $(LIBDIRS) $(LIBS) $(OPTIONS)


my_dlink.o:
	$(NVCC) $(INCDIRS) dydt.o jacob.o gpu_memory.o solver_generic.o solver_interface.o radau2a.o radau2a_init.o radau2a_props.o mechanism.o complexInverse.o inverse.o dengo_solver.o check_memory_usuage.o $(LIBDIR) $(LIBS) -dlink -o my_dlink.o

test: dydt.o jacob.o gpu_memory.o solver_generic.o solver_interface.o radau2a.o radau2a_init.o radau2a_props.o mechanism.o check_memory_usuage.o inverse.o complexInverse.o dengo_solver.o my_dlink.o
	$(CC) -fPIC -fopenmp -O3 -Xlinker -rpath /usr/local/cuda/lib64 dydt.o jacob.o gpu_memory.o solver_generic.o solver_interface.o complexInverse.o inverse.o radau2a.o radau2a_init.o radau2a_props.o mechanism.o dengo_solver.o check_memory_usuage.o  my_dlink.o $(LIBDIRS) $(LIBS) -o test

test_modules_dlink.o:
	$(NVCC) $(INCDIRS) dydt.o jacob.o gpu_memory.o solver_generic.o solver_interface.o radau2a.o radau2a_init.o radau2a_props.o mechanism.o complexInverse.o inverse.o dengo_solver.o test_modules.o $(LIBDIR) $(LIBS) -dlink -o test_modules_dlink.o


test_device_modules: dydt.o jacob.o gpu_memory.o solver_generic.o solver_interface.o radau2a.o radau2a_init.o radau2a_props.o mechanism.o test_modules.o inverse.o complexInverse.o dengo_solver.o test_modules_dlink.o
	$(CC) -fPIC -fopenmp -O3 -Xlinker -rpath /usr/local/cuda/lib64 dydt.o jacob.o gpu_memory.o solver_generic.o solver_interface.o complexInverse.o inverse.o radau2a.o radau2a_init.o radau2a_props.o mechanism.o dengo_solver.o test_modules.o  test_modules_dlink.o $(LIBDIRS) $(LIBS) -o test_device_modules


clean:
	rm -f *o
	rm -f test
	rm -f test_device_modules
