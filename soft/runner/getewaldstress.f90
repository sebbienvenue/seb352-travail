!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - prediction.f90
!!
      subroutine getwaldstress(max_num_atoms,num_atoms,&
         lattice,xyzstruct,nnelecforce,nnstress_ewald)
!!
      use fileunits
!!
      implicit none
!!
      integer max_num_atoms        ! in
      integer num_atoms            ! in
      integer i1,i2,i3,i4          ! internal
!!
      real*8 xyzstruct(3,max_num_atoms)     ! in
      real*8 nnelecforce(3,max_num_atoms)  ! in
      real*8 lattice(3,3)                   ! in
      real*8 nnstress_ewald(3,3)            ! out
      real*8 deltaxj                        ! internal
      real*8 deltayj                        ! internal
      real*8 deltazj                        ! internal
!!
!!
      nnstress_ewald(:,:)=0.0d0
!!
      write(ounit,*)'Ewald stress is not implemented'
      stop
!!  
      do i1=1,3
        do i2=1,3
          do i3=1,num_atoms
            do i4=1,num_atoms ! must be over neighbors (neighborlist for charge symfunctions)
!!
!! how do we truncate the R_ij? A cutoff is probably not ok for long range ewald forces, so our standard neighbor list is not applicable
!!
!!          deltaxj=-1.d0*(xyzstruct(1,i3)-lstb(i4,1))
!!          deltayj=-1.d0*(xyzstruct(2,i3)-lstb(i4,2))
!!          deltazj=-1.d0*(xyzstruct(3,i3)-lstb(i4,3))
!!
!!
!!              nnstress_ewald(i1,i2)=nnstress_ewald(i1,i2)&
!!              +

            enddo ! i4
          enddo ! i3 
        enddo ! i2
      enddo ! i1
!!
      return
      end
