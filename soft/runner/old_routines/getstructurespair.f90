!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - fittingpair.f90
!! - geterrorpair.f90
!! - preconditionpair.f90
!!
      subroutine getstructurespair(unit,npoints,&
           num_atoms_list,lattice_list,&
           xyzstruct_list,atomcharge_list,zelem_list,&
           lperiodic_list)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer unit 
      integer i1
      integer npoints                               ! in
      integer num_atoms_list(nblock)                ! input
      integer num_atoms                             ! internal
      integer zelem(max_num_atoms)                  ! internal
      integer zelem_list(nblock,max_num_atoms)      ! out
!!
      real*8 atomcharge(max_num_atoms)              ! internal
      real*8 atomcharge_list(nblock,max_num_atoms)  ! out
      real*8 atomenergy(max_num_atoms)              ! internal
      real*8 atomenergy_list(nblock,max_num_atoms)  ! internal because not needed now 
      real*8 lattice_list(3,3,nblock)               ! out
      real*8 xyzstruct(3,max_num_atoms)             ! internal
      real*8 xyzstruct_list(3,max_num_atoms,nblock) ! out
      real*8 xyzforce(3,max_num_atoms)              ! internal
!!
      logical lperiodic_list(nblock)                ! out
      logical lperiodic                             ! internal
!!
!!
      do i1=1,npoints
!!
      lperiodic     =.false.
      zelem(:)      =0
      xyzstruct(:,:)=0.0d0
      xyzforce(:,:) =0.0d0
!!
      num_atoms = num_atoms_list(i1)
!!
!! read one structure from trainstruct.data
      call getonestructure(unit,num_atoms,zelem,&
        lattice_list(1,1,i1),xyzstruct,xyzforce,atomcharge,atomenergy,&
        lperiodic)
!!
      zelem_list(i1,:)      = zelem(:)
      xyzstruct_list(:,:,i1)= xyzstruct(:,:)
      lperiodic_list(i1)    = lperiodic
      atomcharge_list(i1,:) = atomcharge(:)
      atomenergy_list(i1,:) = atomenergy(:)
!!
      enddo ! i1
!!
      return
      end
