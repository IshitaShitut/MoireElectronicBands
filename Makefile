SRC_DIR = ./SRC/

MODS=$(wildcard $(SRC_DIR)mod*.F90)
MODOBJS=$(patsubst %.F90, %.o,$(MODS))

SRCS=$(filter-out $(MODS), $(wildcard $(SRC_DIR)*.F90))
OBJS=$(patsubst %.F90, %.o,$(SRCS))

EXE = moire_bands.x

FC = h5pfc

BLAS_LIBS = -L${AMDBLIS_ROOT}/lib -lblis-mt
LAPACK_LIBS = -L${AMDLIBFLAME_ROOT}/lib -lflame
SCALAPACK_LIBS = -L${AMDSCALAPACK_ROOT}/lib -lscalapack 

#BLAS_LIBS = 
#LAPACK_LIBS = -L${OPENBLAS_ROOT}/lib -lopenblas 
#SCALAPACK_LIBS = -L${NETLIB_SCALAPACK_ROOT}/lib -lscalapack

FLFLAGS = ${BLAS_LIBS} ${LAPACK_LIBS} ${SCALAPACK_LIBS}

FCFLAGS = -Ofast -I$(SRC_DIR) -ffree-form

CPP = cpp -P 
CPPFLAGS = -D__KPOOL -D__DEBUG


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