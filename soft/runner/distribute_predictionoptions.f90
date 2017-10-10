!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine distribute_predictionoptions()
!!
      use mpi_mod
      use predictionoptions
!!
      implicit none
!!
      call mpi_bcast(enforcetotcharge,1,mpi_integer,0,mpi_comm_world,mpierror)
!!
      call mpi_bcast(lwritepdb,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lwritexyz,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lwritepov,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lwritepw,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lpreparemd,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(ldoforces,1,mpi_logical,0,mpi_comm_world,mpierror)

      return
      end
