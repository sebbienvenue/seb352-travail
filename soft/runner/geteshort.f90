!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate nneshort for a set of npoints structures 

!! called by: 
!! - getshortenergies_para.f90
!! - optimize_short_combined.f90
!!
!! very similar to getcharges.f90
!!
      subroutine geteshort(ndim,npoints,&
        zelem_local,num_atoms_local,&
        symfunction_local,nneshort_local_list)
!!
      use fileunits
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer ndim                                                    ! in
      integer npoints                                                 ! in
      integer zelem_local(ndim,max_num_atoms)                         ! in
      integer zelem(max_num_atoms)                                    ! internal
      integer num_atoms_local(ndim)                                   ! in
!!
      integer i1
!!
      real*8 symfunction_local(maxnum_funcvalues_short_atomic,max_num_atoms,ndim)  ! in
      real*8 nneshort_local                                           ! internal
      real*8 nneshort_local_list(ndim)                                ! out
      real*8 nnatomenergy(max_num_atoms)                              ! internal
!!
!!
      do i1=1,npoints
!!
        zelem(:)=zelem_local(i1,:)
!! calculate the short-range contribution
        call calconeshort(num_atoms_local(i1),&
          zelem,&
          symfunction_local(1,1,i1),nneshort_local,nnatomenergy)
!!
!! normalize nneshort_local to energy per atom
        nneshort_local=nneshort_local/dble(num_atoms_local(i1))
        nneshort_local_list(i1)=nneshort_local
!!
      enddo ! i1
!!
      return
      end
