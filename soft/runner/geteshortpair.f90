!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! called by: 
!! - getshortenergies_para.f90
!! - optimize_short_combinedpair.f90
!!
!!
      subroutine geteshortpair(npoints,&
       zelemp_local,num_atoms_local,num_pairs_local,&
       symfunctionp_local,nneshort_local)
!!
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer npoints
      integer zelemp(2,max_num_pairs)
      integer i1
      integer num_atoms_local(npoints)
      integer num_pairs_local(npoints)
      integer zelemp_local(2,npoints,max_num_pairs)
!!
      real*8 nneshort                                              ! internal
      real*8 nneshort_local(npoints)                                ! out
      real*8 nnatomenergy(max_num_pairs)                           ! internal
      real*8 symfunctionp_local(maxnum_funcvalues_short_pair,max_num_pairs,npoints)
!!
!!
      do i1=1,npoints
!!
!! calculate the short-range contribution
        zelemp(1,:)=zelemp_local(1,i1,:)            
        zelemp(2,:)=zelemp_local(2,i1,:)            
!!
        call calconeshortpair(num_pairs_local(i1),&
          zelemp,symfunctionp_local(1,1,i1),nneshort,&
          nnatomenergy)
!!
!! normalize nneshort to energy per atom
        nneshort=nneshort/dble(num_atoms_local(i1))
        nneshort_local(i1)=nneshort
!!
      enddo ! i1
!!
      return
      end
