!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine writeheader()
!!
      use mpi_mod
      use fileunits
!!
      implicit none
!!
      if(mpirank.eq.0)then
        if(ounit.ne.6)then
          open(ounit,file='runner.out',form='formatted',status='replace')
        endif
        open(debugunit,file='debug.out',form='formatted',status='replace')
!!
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'---------------------- Welcome to the -----------------------'
        write(ounit,*)'Ruhr-University Neural Network Energy Representation - RuNNer'
        write(ounit,*)'-----------------   (c) Dr. Joerg Behler    -----------------'
        write(ounit,*)'------------  Lehrstuhl fuer Theoretische Chemie  -----------'
        write(ounit,*)'------------       Ruhr-Universitaet Bochum       -----------'
        write(ounit,*)'------------        44780 Bochum, Germany         -----------'
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'When using RuNNer, please cite the following papers:'
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'For general high-dimensional NNs:'
        write(ounit,*)'J. Behler and M. Parrinello, Phys. Rev. Lett. 98, 146401 (2007).'
        write(ounit,*)'J. Behler, J. Chem. Phys. 134, 074106 (2011).'
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'For high-dimensional NNs including electrostatics:'
        write(ounit,*)'N. Artrith, T. Morawietz, and J. Behler, Phys. Rev. B 83, 153101 (2011).'
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Reviews on NN potentials:'
        write(ounit,*)'C.M. Handley, and P.L.A. Popelier, J. Phys. Chem. A 114, 3371 (2010).'
        write(ounit,*)'J. Behler, Phys. Chem. Chem. Phys. 13, 17930 (2011).'
        write(ounit,*)'J. Behler, J. Phys.: Condens. Matter 26, 183001 (2014).'
        write(ounit,*)'J. Behler, Int. J. Quant. Chem. 115, 1032 (2015).'
        write(ounit,*)'-------------------------------------------------------------'
!!
        call generalinfo()
        call compileinfo()
      endif ! mpirank.eq.0
!!
      return
      end
