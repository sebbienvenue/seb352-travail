!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! many routines 
!!
      subroutine mpifitdistribution(npoints,nstruct,n_start,n_end)
!!
      use mpi_mod
      use fileunits
!!
      implicit none
      
!!
      integer npoints            ! in
      integer nstruct            ! out
      integer, dimension(:), allocatable :: nstruct_list ! internal
      integer icount             ! internal
      integer i1                 ! internal
      integer n_start            ! out
      integer n_end              ! out 
!!
!! decide here how many structures should be calculated by which process in this block
!! nstruct is the number of structures to be handeled by the current process
      if((mpirank+1).le.(mod(npoints,mpisize)))then
        nstruct=int(npoints/mpisize)+1
      else
        nstruct=int(npoints/mpisize)
      endif
!!      write(ounit,*)'process ',mpirank,' calculates ',nstruct,' structures'
!!
      allocate(nstruct_list(mpisize))  !'
!!
      nstruct_list(:)=int(npoints/mpisize)
      icount=mod(npoints,mpisize)
      do i1=1,mpisize
        if(icount.gt.0)then
          icount=icount-1
          nstruct_list(i1)=nstruct_list(i1)+1
        endif
      enddo
!! calculate the first and last structure to be calculated by this process
      n_start=1
      do i1=1,mpirank
        n_start=n_start+nstruct_list(i1)
      enddo
      n_end=n_start+nstruct-1
!!
      deallocate(nstruct_list)
!!
      return
      end
