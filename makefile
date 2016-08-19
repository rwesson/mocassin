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


F90  = /usr/local/cluster/mpich-1.2.5/bin/mpif90-8
LIBS = 
OPT1 = -O3 -mp1 -static -static-libcxa
OPT2 = -CB
OPT3 = 

mocassin:
	$(F90) $(OPT1) -o mocassin $(source1) 

mocassinWarm:
	$(F90) $(OPT1) -o mocassinWarm $(source2) $(LIBS)

mocassinOutput:
	$(F90) $(OPT1) -o mocassinOutput $(source3) $(LIBS)

clean:
	/bin/rm *.o *~ *.mod *.d *.il mocassin mocassinWarm mocassinOutput


