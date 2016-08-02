# standard
FC = mpif90
LD = mpif90
FFLAGS += -cpp -fno-range-check -Jsource/ -ffree-line-length-0 -lm -DPREFIX=\"${PREFIX}\"

#IBM
#FC = mpxlf90_r
#LD = mpxlf90_r
#FFLAGS += -O3 -qstrict -qrealsize=4 -bmaxdata:0x40000000 -qmaxmem=-1 -cpp -DPREFIX=\"${PREFIX}\"

#SUN
#FC = mpf90
#LD = mpf90
#FFLAGS += -fast -xarch=v9b -ftrap=%none -lmpi -cpp -DPREFIX=\"${PREFIX}\"

#SGI
#FC = f90
#LD = f90
#FFLAGS += -64 -C -mpio -OPT:Olimit=3495 -O3 -lmpi -cpp -DPREFIX=\"${PREFIX}\"

PREFIX=/usr
MANDIR=${DESTDIR}${PREFIX}/share/man/man1
SOURCES = source/constants_mod.o source/vector_mod.o source/common_mod.o source/interpolation_mod.o \
	source/set_input_mod.o source/hydro_mod.o source/ph_mod.o source/composition_mod.o \
	source/continuum_mod.o source/ionization_mod.o source/pathIntegration_mod.o \
	source/grid_mod.o source/dust_mod.o source/emission_mod.o source/photon_mod.o  \
	source/update_mod.o source/output_mod.o source/iteration_mod.o

ifeq ($(CO),debug) #to show all compiler warnings
  FFLAGS += -fbounds-check -Wall -Wuninitialized -g -pg #-ffpe-trap=zero,overflow,invalid,underflow,denormal -fbacktrace -fcheck=all
else ifeq ($(CO),valgrind)
  FFLAGS += -g
else ifeq ($(CO),gprof)
  FFLAGS += -pg
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
	test -e ${DESTDIR}${PREFIX}/share/mocassin || mkdir -p ${DESTDIR}${PREFIX}/share/mocassin
	test -e ${DESTDIR}${PREFIX}/bin || mkdir -p ${DESTDIR}${PREFIX}/bin
	test -e ${MANDIR} || mkdir -p ${MANDIR}
	cp -r data/ ${DESTDIR}${PREFIX}/share/mocassin
	cp -r dustData/ ${DESTDIR}${PREFIX}/share/mocassin
	cp -r benchmarks/ ${DESTDIR}${PREFIX}/share/mocassin
	cp -r examples/ ${DESTDIR}${PREFIX}/share/mocassin
	install -m 644 man/mocassin.1 ${MANDIR}
	gzip -f ${MANDIR}/mocassin.1
	ln -s -f ${MANDIR}/mocassin.1.gz ${MANDIR}/mocassinWarm.1.gz
	ln -s -f ${MANDIR}/mocassin.1.gz ${MANDIR}/mocassinOutput.1.gz
	ln -s -f ${MANDIR}/mocassin.1.gz ${MANDIR}/mocassinPlot.1.gz
	install mocassin ${DESTDIR}${PREFIX}/bin
	install mocassinWarm ${DESTDIR}${PREFIX}/bin
	install mocassinPlot ${DESTDIR}${PREFIX}/bin
	install mocassinOutput ${DESTDIR}${PREFIX}/bin

uninstall:
	rm -f ${DESTDIR}${PREFIX}/bin/mocassin*
	rm -f ${MANDIR}/mocassin*.1.gz
	rm -rf ${DESTDIR}${PREFIX}/share/mocassin
