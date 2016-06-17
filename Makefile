FC  = mpif90
LIBS =	-lm
FFLAGS = -fno-range-check -Jsource/ -ffree-line-length-0
MANDIR=${DESTDIR}/usr/share/man/man1
SOURCES = source/infnan.f90 source/constants_mod.f90 source/vector_mod.f90 source/common_mod.f90 source/interpolation_mod.f90 \
	source/set_input_mod.f90 source/hydro_mod.f90 source/ph_mod.f90 source/composition_mod.f90 \
	source/continuum_mod.f90 source/ionization_mod.f90 source/pathIntegration_mod.f90 \
	source/grid_mod.f90 source/dust_mod.f90 source/emission_mod.f90 source/photon_mod.f90  \
	source/update_mod.f90 source/output_mod.f90 source/iteration_mod.f90

ifeq ($(CO),debug) #to show all compiler warnings
  FFLAGS += -fbounds-check -Wall -Wuninitialized #-ffpe-trap=zero,overflow,invalid,underflow,denormal
else ifeq ($(CO),debug2) #for profiling
  FFLAGS += -g -pg
endif

.PHONY: all clean new install uninstall

all: mocassin mocassinWarm mocassinOutput mocassinPlot

new: clean all

mocassin:
	$(FC) $(FFLAGS) -o mocassin $(SOURCES) source/mocassin.f90 $(LIBS)

mocassinWarm:
	$(FC) $(FFLAGS) -o mocassinWarm $(SOURCES) source/mocassinWarm.f90 $(LIBS)

mocassinOutput:
	$(FC) $(FFLAGS) -o mocassinOutput $(SOURCES) source/mocassinOutput.f90 $(LIBS)

mocassinPlot:
	$(FC) $(FFLAGS) -o mocassinPlot $(SOURCES) source/mocassinPlot.f90 $(LIBS)

clean:
	/bin/rm -f source/*.o *~ source/*.mod mocassin mocassinWarm mocassinOutput mocassinPlot

install:
	test -e ${DESTDIR}/usr/share/mocassin || mkdir -p ${DESTDIR}/usr/share/mocassin
	test -e ${DESTDIR}/usr/bin || mkdir -p ${DESTDIR}/usr/bin
	test -e ${MANDIR} || mkdir -p ${MANDIR}
	cp -r data/ ${DESTDIR}/usr/share/mocassin
	cp -r dustData/ ${DESTDIR}/usr/share/mocassin
	cp -r benchmarks/ ${DESTDIR}/usr/share/mocassin
	cp -r examples/ ${DESTDIR}/usr/share/mocassin
	install -g 0 -o 0 -m 644 man/mocassin.1 ${MANDIR}
	gzip -f ${MANDIR}/mocassin.1
	ln -s -f ${MANDIR}/mocassin.1.gz ${MANDIR}/mocassinWarm.1.gz
	ln -s -f ${MANDIR}/mocassin.1.gz ${MANDIR}/mocassinOutpot.1.gz
	ln -s -f ${MANDIR}/mocassin.1.gz ${MANDIR}/mocassinPlot.1.gz
	install mocassin ${DESTDIR}/usr/bin
	install mocassinWarm ${DESTDIR}/usr/bin
	install mocassinPlot ${DESTDIR}/usr/bin
	install mocassinOutput ${DESTDIR}/usr/bin

uninstall:
	rm -f ${DESTDIR}/usr/bin/mocassin*
	rm -f ${MANDIR}/mocassin*.1.gz
	rm -rf ${DESTDIR}/usr/share/mocassin
