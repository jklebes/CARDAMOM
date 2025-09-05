EXAMPLE_BIN_DIR=./

FC=gfortran
FCFLAGS=-O3 -ffree-line-length-0
LDFLAGS= -p -pg

TARGETS=$(EXAMPLE_BIN_DIR)/cardamom.exe

SRC_DIR=LIBRARY/CARDAMOM_F
SRC_BIN_DIR=$(SRC_DIR)/executable
SRC=$(SRC_DIR)/misc/math_functions.f90 \
    $(SRC_DIR)/misc/oksofar.f90 \
    $(SRC_DIR)/misc/brent_zero/brent_zero.f90 \
    $(SRC_DIR)/model/DALEC.A1.C1.D2.F2.H2.P1.004/src/DALEC.A1.C1.D2.F2.H2.P1.004.f90 \
    $(SRC_DIR)/general/cardamom_structures.f90 \
    $(SRC_DIR)/method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 \
    $(SRC_DIR)/method/MHMCMC/MCMC_FUN/MHMCMC_StressTests.f90 \
    $(SRC_DIR)/model/DALEC.A1.C1.D2.F2.H2.P1.004/src/DALEC.A1.C1.D2.F2.H2.P1.004_PARS.f90 \
    $(SRC_DIR)/general/cardamom_io.f90 \
    $(SRC_DIR)/method/MHMCMC/MCMC_FUN/MHMCMC.f90 \
    $(SRC_DIR)/model/DALEC.A1.C1.D2.F2.H2.P1.004/likelihood/MODEL_LIKELIHOOD.f90 \
    $(SRC_DIR)/general/cardamom_main.f90

all : $(TARGETS)

$(EXAMPLE_BIN_DIR)/cardamom.exe : $(SRC_BIN_DIR)/cardamom.exe
	cp $< $@

$(SRC_BIN_DIR)/cardamom.exe : $(SRC)
	$(FC) $(FCFLAGS) $^ -o $@

.PHONY: clean
clean :
	rm -f *.o *.mod $(EXAMPLE_BIN_DIR)/cardamom.exe $(SRC_BIN_DIR)/cardamom.exe
