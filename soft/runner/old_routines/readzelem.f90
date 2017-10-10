!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Jovan Jose 2010
!######################################################################

!! called by:
!! - fittingpair.f90 (only once for pair case)
!!
      subroutine  readzelem(npoints,nblock,max_num_atoms,num_atoms_list,zelem_list)
!!
      use fileunits
!!
      implicit none
!!     
      integer npoints                                               ! in
      integer nblock                                                ! in
      integer max_num_atoms                                         ! in
      integer zelem_list(nblock,max_num_atoms)                      ! out 
      integer num_atoms_list(nblock)                                ! out
      integer i1,i2,i3                                              ! internal
      real*8  dummy
!!
      
      open(munit,file='trainstruct.data')
       do i1=1,npoints
         read(munit,*)dummy

         do i2= 1,num_atoms_list(i1)
           read(munit,*)zelem_list(i1,i2) 
           write(ounit,*)zelem_list(i1,i2),'check',i1,i2
         end do

       end do
      close(munit)

!!
!!
      return
      end
