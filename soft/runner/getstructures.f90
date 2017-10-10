!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - fitting_batch.f90
!! - fittingpair.f90
!! - geterror.f90
!! - geterrorpair.f90
!! - precondition.f90
!! - preconditionpair.f90
!!
      subroutine getstructures(unit,npoints)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!
      integer unit                                  ! in
      integer npoints                               ! in
      integer zelem(max_num_atoms)                  ! internal
      integer i1
!!
      real*8 xyzstruct(3,max_num_atoms)             ! internal
      real*8 totalforce(3,max_num_atoms)            ! internal
      real*8 atomcharge(max_num_atoms)              ! internal
      real*8 atomenergy(max_num_atoms)              ! internal
!!
      logical lperiodic                             ! internal
!!
      do i1=1,npoints
!!
        lperiodic       = .false.
        zelem(:)        = 0
        xyzstruct(:,:)  = 0.0d0
        totalforce(:,:) = 0.0d0
!!
!! read one structure from trainstruct.data
        call getonestructure(unit,num_atoms_list(i1),zelem,&
          lattice_list(1,1,i1),xyzstruct,totalforce,atomcharge,atomenergy,&
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
