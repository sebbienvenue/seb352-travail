!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine, works for a list of structures and for a single structure

!! called by:
!! - getsymmetryfunctions.f90
!! - getpairsymmetryfunctions.f90
!! - predict.f90
!!
      subroutine removeatoms(ndim,npoints, &
             num_atoms_local,zelem_local,&
             num_atoms_element_local,&
             totalenergy_local,atomenergy_local)
!!
      use globaloptions
!!
      implicit none
!!
      integer ndim                                 ! in
      integer npoints                              ! in
      integer num_atoms_local(ndim)               ! in
      integer num_atoms_element_local(ndim,nelem) ! in  
      integer zelem_local(ndim,max_num_atoms)     ! in
      integer i1,i2                                ! internal
!!
      real*8 totalenergy_local(ndim)              ! in/out
      real*8 atomenergy_local(ndim,max_num_atoms) ! in/out
!!
!!
!! remove atomic energies from total energies
      do i1=1,npoints
        do i2=1,nelem
          totalenergy_local(i1)=totalenergy_local(i1)&
          -dble(num_atoms_element_local(i1,i2))*atomrefenergies(i2)
        enddo
      enddo
!!
!! remove atomic energies from atomic energy contributions ( Not Implimented Yet !!)
      if(luseatomenergies)then 
        do i1=1,npoints  ! loop over all points
          do i2=1,num_atoms_local(i1) ! loop over all atoms of each point
            atomenergy_local(i1,i2)=atomenergy_local(i1,i2)&
              -atomrefenergies(elementindex(zelem_local(i1,i2)))
          enddo
        enddo ! i1
      endif
!!
      return
      end
