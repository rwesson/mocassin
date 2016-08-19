source1  = source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassin.f90 

source2  = source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinWarm.f90 

source3  = source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinOutput.f90 

source4  = source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinPlot.f90 
 
F90  = mpif90
LIBS = -lm
OPT1 = -O3 
OPT2 = -check all -traceback -g
#OPT2 = -CB -g 
XTRA =  -I/software/mpich2/include -L/software/mpich2/lib 

mocassin:
	$(F90) $(OPT1) -o mocassin $(source1) $(LIBS) $(XTRA)

mocassinWarm:
	$(F90) $(OPT1) -o mocassinWarm $(source2) $(LIBS) $(XTRA)

mocassinOutput:
	$(F90) $(OPT1) -o mocassinOutput $(source3) $(LIBS) $(XTRA)

mocassinPlot:
	$(F90) $(OPT2) -o mocassinPlot $(source4) $(LIBS) $(XTRA)

clean:
	/bin/rm *.o *~ *.mod mocassin mocassinWarm mocassinOutput mocassinPlot


 







