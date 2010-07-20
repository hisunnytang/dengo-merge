all : test_primordial_solver run_test

clean :
	rm -f simple_cvode_solver/primordial_cvode_solver.c

simple_cvode_solver/primordial_cvode_solver.c : \
        simple_cvode_solver/cvode_solver.c.template \
        primordial_chemistry/*.py
	@echo "Regenerating solver from template"
	@python2.6 primordial_chemistry/write_cvode_solver.py

test_primordial_solver : simple_cvode_solver/primordial_cvode_solver.c 
	@echo "Recompiling"
	@gcc -Isimple_cvode_solver -I/usr/local/include -L/usr/local/include \
		-lhdf5_hl -lhdf5 simple_cvode_solver/primordial_cvode_solver.c \
		-o test_primordial_solver

run_test : test_primordial_solver
	./test_primordial_solver
