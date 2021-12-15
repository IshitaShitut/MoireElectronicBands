SRC_DIR = ./SRC/

MODS=$(wildcard $(SRC_DIR)mod*.F90)
MODOBJS=$(patsubst %.F90, %.o,$(MODS))

SRCS=$(filter-out $(MODS), $(wildcard $(SRC_DIR)*.F90))
OBJS=$(patsubst %.F90, %.o,$(SRCS))

EXE = moire_bands.x

FC = h5pfc

BLAS_LIBS = 
LAPACK_LIBS = -L${OPENBLAS_ROOT}/lib -lopenblas 
SCALAPACK_LIBS = -L${NETLIB_SCALAPACK_ROOT}/lib -lscalapack

FLFLAGS = ${BLAS_LIBS} ${LAPACK_LIBS} ${SCALAPACK_LIBS}

FCFLAGS = -g -fcheck=bounds -fbacktrace -O3 -I$(SRC_DIR) -fPIC -Wall 

CPP = cpp -P 
CPPFLAGS = -D__DEBUG #-D__KPOOL


all: $(MODOBJS) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $(OBJS) $(MODOBJS) $(FLFLAGS)

$(MODOBJS): %.o: %.F90
# $(CPP) $(CPPFLAGS) $< | sed '/^#pragma/d' > $*.f90
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(FCFLAGS) -J $(SRC_DIR) -c $*.f90 -o $@
# $(FC) $(FCFLAGS) -module $(SRC_DIR) -c $*.f90 -o $@

$(OBJS): %.o: %.F90
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(FCFLAGS) -c $*.f90 -o $@ 

debug:
	@echo "SRCS=$(SRCS)"
	@echo "OBJS=$(OBJS)"
	@echo "MODS=$(MODS)"
	@echo "MODOBJS=$(MODOBJS)"
	@echo "EXE=$(EXE)"
	@echo "FC =$(FC)"
	@echo "FLFLAGS=$(FLFLAGS)"
	@echo "FCFLAGS=$(FCFLAGS)"
	@echo "CPP=$(CPP)"
	@echo "CPPFLAGS=$(CPPFLAGS)"

.PHONY: clean
clean:
	rm -f $(SRC_DIR)*.o $(SRC_DIR)*.mod $(EXE) $(SRC_DIR)*.f90