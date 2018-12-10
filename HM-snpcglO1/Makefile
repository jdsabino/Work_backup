EXTROOT =   /home/colin/Desktop/Bachelorarbeit/local

EXTCORE =   /home/colin/Desktop/Joao/cfl/core

M_INIT =           cfl_INIT
M_GRID =           cfl_grid
M_FFTW =           cfl_fftw
M_HDF5 =           cfl_hdf5
M_POT =            cfl_potential \
                   cfl_potential/potential/harmonic \
                   cfl_potential/potential/box \
                   cfl_potential/potential/ioffe \
		   cfl_potential/potential/double_box \
                   cfl_potential/disorder/truncated_white_noise
M_OBS =            cfl_observables
M_CSPLIT =         cfl_coherent_splitting
M_PROP_TOOLS =     cfl_propagator/tools
M_PROP =           cfl_propagator \
		   cfl_propagator/snpcglO1 \
                   cfl_propagator/snpcglO2
                   
MODULES = $(M_INIT) $(M_GRID) $(M_FFTW) $(M_HDF5) $(M_POT) $(M_OBS) $(M_CSPLIT) $(M_PROP_TOOLS) $(M_PROP)


SCR_DIR = $(addprefix $(EXTCORE)/,$(MODULES))
VPATH = $(SCR_DIR)

O_MAIN =           cfl_sim.o \
                   cfl_sysParam.o\
                   usr_update.o \
                   usr_multiSim.o \
                   usr_chData.o
O_INIT =           cfl_init.o
O_GRID =           cfl_complex_grid.o \
	           cfl_double_grid.o
O_FFTW =           cfl_fftw.o
O_HDF5 =           cfl_hdf5_complex_write.o \
                   cfl_hdf5_complex_read.o \
                   cfl_hdf5_double_write.o \
                   cfl_hdf5_double_read.o
O_POT =            cfl_potential.o \
                   cfl_potential_IO.o \
                   pot_harmonic.o  \
                   pot_harmonic_IO.o \
                   pot_box.o \
                   pot_box_IO.o \
                   pot_ioffe.o \
                   pot_ioffe_IO.o \
                   pot_double_box.o \
                   pot_double_box_IO.o \
		   dis_truncWN.o \
                   dis_truncWN_IO.o
O_OBS =            cfl_PSF.o \
                   cfl_density_phase.o \
                   cfl_energy.o
O_CSPLIT =         cfl_coherent_splitting.o
O_PROP_TOOLS =     cfl_propagator_tools.o \
                   cfl_laplace.o \
                   cfl_bogoliubov.o
O_PROP_O1 =        snpcglO1.o \
                   snpcglO1_init.o \
                   snpcglO1_mask.o \
                   snpcglO1_chempot.o \
                   snpcglO1_tsfp.o \
                   snpcglO1_ts_rprop.o \
                   snpcglO1_bogoliubov.o
O_PROP_O2 =        snpcglO2.o \
                   snpcglO2_init.o \
                   snpcglO2_mask.o \
                   snpcglO2_chempot.o \
                   snpcglO2_tsfp.o \
                   snpcglO2_ts_rprop.o \
                   snpcglO2_bogoliubov.o
                   
OBJECT = $(O_INIT) $(O_MAIN) $(O_GRID) $(O_FFTW) $(O_HDF5) $(O_POT) $(O_OBS) $(O_CSPLIT) $(O_PROP_TOOLS) $(O_PROP_O1) $(O_PROP_O2)


EXTLIB =    -L$(EXTROOT)/lib
INCLUDE =   -I$(EXTROOT)/include $(INCLUDE_CFL)
INCLUDE_CFL = $(addprefix -I$(EXTCORE)/,$(MODULES))

CC =       gcc-6 $(INCLUDE)
CFLAGS =   -Wl,-rpath=$(EXTROOT)/lib -Wl,--no-as-needed -g -Wall -std=c99 -O3 \
	   $(EXTLIB) -lgsl -lgslcblas -fopenmp -lfftw3 -lm -lfftw3_omp -lhdf5 \
                     -lconfig
                     

%.o : %.c *.h
	$(CC) $(CFLAGS) -c $<

sim : $(OBJECT)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	-rm -f *.o *~ sim

cleanall:
	-rm -f *.o *~ sim data_out/*.*

cleandata:
	-rm -f data_out/*.*

