!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by
!! - main.f90
!!
      subroutine mpistart(ounit,nproc,lmpi)

      use mpi_mod
 
      implicit none

      integer ounit
      integer nproc

      logical lmpi(nproc)

      lmpi(:)=.false.

      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
      lmpi(mpirank+1)=.true.
!!
      if(mpirank.eq.0)then
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,'(a29,i4,a6)')&
        ' Parallel run requested with ',mpisize,' cores'
      endif ! mpirank.eq.0
!!
!! This works for all processes:
!!      write(*,'(a29,i4,a6)')&
!!       ' Parallel run requested with ',mpisize,' cores'

!!      if(mpirank.eq.0)then
!!      write(ounit,*)'Process number  ',mpirank,' ready'
!!      endif
!!
!! This works for all processes:
!!      write(*,*)'Process number  ',mpirank,' ready'


      return
      end
