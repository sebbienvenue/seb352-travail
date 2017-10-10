# optimization/compiler settings
# Ifort - for production
compiler_type=intelem
compiler=ifort
fflags='-xS -free -ipo'
opt='-O3'

# GFortran - just for checking
#compiler_type=gnu95
#compiler=gfortran-4.7
#fflags='-ffree-form -g -fbacktrace -fcheck=all'
#opt='-O0'

# leave as is
mpi_compiler=mpif90
