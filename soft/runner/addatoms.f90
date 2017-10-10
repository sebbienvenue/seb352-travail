!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!! - prediction.f90
!! - predictionpair.f90
!!
      subroutine addatoms(num_atoms,&
             zelem,num_atoms_element,&
             totalenergy,atomenergy)
!!
      use globaloptions
!!
      implicit none
!!
      integer num_atoms                            ! in
      integer num_atoms_element(nelem)             ! in  
      integer zelem(max_num_atoms)                 ! in
      integer i1                                   ! internal
!!
      real*8 totalenergy                           ! in/out
      real*8 atomenergy(max_num_atoms)             ! in/out
!!
!!
!! add atomic energies to total energy
        do i1=1,nelem
          totalenergy=totalenergy&
          +dble(num_atoms_element(i1))*atomrefenergies(i1)
        enddo
!!
!!
!! add atomic energies to atomic energy contributions
        do i1=1,num_atoms ! loop over all atoms 
          atomenergy(i1)=atomenergy(i1)&
            +atomrefenergies(elementindex(zelem(i1)))
        enddo
!!
      return
      end
