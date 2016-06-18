FC = mpif90
LD = mpif90
FFLAGS += -fno-range-check -Jsource/ -ffree-line-length-0 -lm
MANDIR=${DESTDIR}/usr/share/man/man1
SOURCES = source/constants_mod.o source/vector_mod.o source/common_mod.o source/interpolation_mod.o \
	source/set_input_mod.o source/hydro_mod.o source/ph_mod.o source/composition_mod.o \
	source/continuum_mod.o source/ionization_mod.o source/pathIntegration_mod.o \
	source/grid_mod.o source/dust_mod.o source/emission_mod.o source/photon_mod.o  \
	source/update_mod.o source/output_mod.o source/iteration_mod.o

ifeq ($(CO),debug) #to show all compiler warnings
  FFLAGS += -fbounds-check -Wall -Wuninitialized #-ffpe-trap=zero,overflow,invalid,underflow,denormal
else ifeq ($(CO),debug2) #for profiling
  FFLAGS += -g -pg
endif

.PHONY: all clean new install uninstall

all: mocassin mocassinWarm mocassinOutput mocassinPlot

new: clean all

%.o: %.f90
	$(FC) $(FFLAGS) $< -c -o $@

mocassin: $(SOURCES) source/mocassin.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

mocassinWarm: $(SOURCES) source/mocassinWarm.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

mocassinOutput: $(SOURCES) source/mocassinOutput.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

mocassinPlot: $(SOURCES) source/mocassinPlot.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

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
