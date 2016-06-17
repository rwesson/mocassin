source1  = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassin.f90 

source2  = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinWarm.f90 

source3  = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinOutput.f90 

source4  = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 \
	source/output_mod.f90 source/iteration_mod.f90 source/mocassinPlot.f90 
 
F90  = mpif90
LIBS =	-lm
OPT1 = -fno-range-check
OPT2 = -CB -g -traceback

mocassin:
	$(F90) $(OPT1) -o mocassin $(source1) $(LIBS)

mocassinWarm:
	$(F90) $(OPT1) -o mocassinWarm $(source2) $(LIBS) 

mocassinOutput:
	$(F90) $(OPT2) -o mocassinOutput $(source3) $(LIBS) 

mocassinPlot:
	$(F90) $(OPT2) -o mocassinPlot $(source4) $(LIBS) 

clean:
	/bin/rm *.o *~ *.mod mocassin mocassinWarm mocassinOutput mocassinPlot


 







