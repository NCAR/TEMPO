FC = gfortran
FCFLAGS = -Dstandalone -fdefault-real-8 -fdefault-double-8

SRC=machine.F90 module_mp_radar.F90 ../module_mp_tempo_params.F90 ../module_mp_tempo_utils.F90 ../module_mp_tempo_main.F90 ../drivers/standalone/module_mp_tempo.F90 driver.F90
OBJ=machine.o module_mp_radar.o module_mp_tempo_params.o module_mp_tempo_utils.o module_mp_tempo_main.o module_mp_tempo.o driver.o

%.o: %.F90
	$(FC) $(FCFLAGS) -c $(SRC)

driver: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.mod *.o driver
