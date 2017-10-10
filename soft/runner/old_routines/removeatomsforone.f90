!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - prediction.f90 
!!
      subroutine removeatomsforone(num_atoms,&
             zelem,&
             num_atoms_element,&
             totalenergy,atomenergy)
!!
      use globaloptions
!!
      implicit none
!!
      integer num_atoms                            ! in
      integer num_atoms_element(nelem)              ! in  
      integer zelem(max_num_atoms)                 ! in
      integer i1,i2                                ! internal
!!
      real*8 totalenergy                           ! in/out
      real*8 atomenergy(max_num_atoms)             ! in/out
!!
!!
!! remove atomic energies from total energies
        do i2=1,nelem
          totalenergy=totalenergy&
          -dble(num_atoms_element(i2))*atomrefenergies(i2)
        enddo
!!
!! remove atomic energies from atomic energy contributions
      if(luseatomenergies)then 
          do i2=1,num_atoms ! loop over all atoms of each point
            atomenergy(i2)=atomenergy(i2)&
              -atomrefenergies(elementindex(zelem(i2)))
          enddo
      endif
!!
      return
      end
