!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - predict.f90
!!
      subroutine writexyz(num_atoms,&
           xyzstruct,elementsymbol)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer num_atoms                ! in
      integer i1,i2                    ! internal
!! 
      real*8 xyzstruct(3,max_num_atoms) ! in
!!
      character*2 elementsymbol(max_num_atoms) ! in
!!
!!
!! Note: xyz file has Angstrom length unit
!!
      open(xyzunit,file='runner.xyz',form='formatted',status='replace')
        write(xyzunit,'(i8)')num_atoms
        write(xyzunit,'(a)')'This file has been generated by RuNNer'
        do i1=1,num_atoms
          write(xyzunit,'(a2,x,3f14.8)')elementsymbol(i1),&
          (xyzstruct(i2,i1)*0.529177d0,i2=1,3)
        enddo
      close(xyzunit)
!!
      return
      end
