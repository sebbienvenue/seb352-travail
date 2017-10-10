!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - getstructures.f90
!! - getstructurespair.f90
!!
      subroutine getonestructure(unit,num_atoms,zelem,&
        lattice,xyzstruct,totalforce,atomcharge,atomenergy,&
        lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer unit                      ! in
      integer num_atoms                 ! in 
      integer zelem(max_num_atoms)      ! out
      integer idummy                    ! internal    
      integer i1,i2
!!
      real*8 lattice(3,3)               ! out
      real*8 xyzstruct(3,max_num_atoms) ! out
      real*8 totalforce(3,max_num_atoms)! out
      real*8 atomcharge(max_num_atoms)  ! out
      real*8 atomenergy(max_num_atoms)  ! out
!!
      logical lperiodic                 ! out 
!!
      read(unit,*)idummy,lperiodic
      if(lperiodic)then
        do i1=1,3
          read(unit,*)(lattice(i1,i2),i2=1,3)
        enddo
      endif
      do i1=1,num_atoms
        read(unit,*)zelem(i1),(xyzstruct(i2,i1),i2=1,3),atomcharge(i1),&
        atomenergy(i1),(totalforce(i2,i1),i2=1,3)
      enddo
!!
      return
      end
