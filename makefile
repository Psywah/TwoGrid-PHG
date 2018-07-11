#PHG_HOME  =/home/wleng/Proj-Deform/phg-ng2/
PHG_HOME  =/home/wleng/Test/phg-0.9.3/
PHG_LIB   =${PHG_HOME}/src/libphg

FILES = Makefile functions.h poisson.c simplest.c .cvsignore *-plot.sh jobs.sh \
	maxwell*.c heat.c non-smooth.c eigen.c elastic.c navier-stokes.[ch] \
	navier-stokes-unsteady.c

#default: test-grid
default: test
#default: test-solverT

all: lib poisson simplest maxwell maxwell-complex heat non-smooth eigen \
	elastic maxwell-eigen navier-stokes

examples.zip: ${FILES}
	@zip -9 -u -y $@ $^

clean:
	-/bin/rm -f *.o

distclean: clean
	-/bin/rm -f *.o ins-flow

matclean:
	-/bin/rm -f *.m *.m.dat

lib:
	@(cd ../src; $(MAKE))


include ${PHG_HOME}/Makefile.inc
VPATH=./src
vpath %.c ./src


#CPPFLAS=${CPPFLAGS} -O0
EXTRAFLAGS=-O0

#################
#   test grid   #
#################

two-grid.o: two-grid.c  
oct-search.o: oct-search.c oct-search.h
test.o: test.c  oct-search.h
pgrid.o: pgrid.c pgrid.h
test: ${PHG_LIB}${LIB_SUFFIX}  test.o two-grid.o oct-search.o  pgrid.o
	@echo "Linking test-grid"
	${LINKER} ${LDFLAGS} -o $@ $^ ${LIBS}    


