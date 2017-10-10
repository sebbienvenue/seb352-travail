##
# compiler to use
##
FC := gfortran-4.7
# FCFLAGS := -O3
FCFLAGS := -O0 -g -pg -Wall -Wconversion-extra -Wextra -fcheck=all -fbacktrace
##
# FC := ifort
# FCFLAGS := -O3 -xS -free -Tf
# FCFLAGS := -O0 -g -pg -traceback -warn all -check all -free -Tf

##
# preprocessor
##
CPPFLAGS := -cpp

##
# linker flags
##
LD := $(FC)
LDFLAGS := -lopenblas

##
# MPI
##
ifeq ($(MPI),1)
	FC := mpif90
	LD := $(FC)
endif
