!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getshortenergies_parapair.f90
!! - scalesymfit_parapair.f90
!!
      subroutine scalesympair(ndim,npoints,&
         num_pairs_local,zelemp_local,symfunctionp_local,&
         minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)
!!
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer ndim                                           ! in
      integer npoints                                        ! in
      integer num_pairs_local(ndim)                          ! in
      integer zelemp_local(2,ndim,max_num_pairs)             ! in
      integer i1,i2,i3                                       ! internal
!!
      real*8 symfunctionp_local(maxnum_funcvalues_short_pair,max_num_pairs,ndim)  ! in/out
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)            ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)            ! in
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)             ! in
!!
!!
      do i1=1,npoints
        do i2=1,num_pairs_local(i1)
          do i3=1,num_funcvalues_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2))) 
            if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
              symfunctionp_local(i3,i2,i1)=symfunctionp_local(i3,i2,i1) &
              - avvalue_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2)),i3) 
            elseif(lscalesym.and..not.lcentersym)then
              symfunctionp_local(i3,i2,i1)=&
             (symfunctionp_local(i3,i2,i1) -  minvalue_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2)),i3))/ &
             (maxvalue_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2)),i3)-&
              minvalue_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2)),i3))&
              *(scmax_short_pair-scmin_short_pair) + scmin_short_pair
            elseif(lscalesym.and.lcentersym)then
              symfunctionp_local(i3,i2,i1)=&
             (symfunctionp_local(i3,i2,i1)-avvalue_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2)),i3))&
             / &
             (maxvalue_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2)),i3)&
              - minvalue_short_pair(pairindex(zelemp_local(1,i1,i2),zelemp_local(2,i1,i2)),i3))
            else
            endif
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
      return
      end
