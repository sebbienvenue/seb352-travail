!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - optimize_short_combinedpair.f90
!!
      subroutine scaledsfuncpair_para(num_pairs,&
        minvalue_short_pair,maxvalue_short_pair,&
        zelemp,dsfuncdxyz_pair,strs_pair)
!!
      use fileunits
      use globaloptions
      use nnshort_pair
!! 
      implicit none
!!
      integer zelemp(2,max_num_pairs) 
      integer i2,i3                                                  ! internal
      integer num_pairs
!!
      real*8 dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3)
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)     ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)     ! in

      real*8 strs_pair(3,3,maxnum_funcvalues_short_pair,max_num_pairs)
!!
      do i2= 1,num_pairs 
        do i3=1,num_funcvalues_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)))
!! scale each symmetry function derivative value for the respective element
!!           dsfuncdxyz_pair(i3,i2,:,:)=dsfuncdxyz_pair(i3,i2,:,:)/&
!!             (maxvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)&
!!             -minvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3))
!! new for variable scaling range
          dsfuncdxyz_pair(i3,i2,:,:)=dsfuncdxyz_pair(i3,i2,:,:)/&
            (maxvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)&
            -minvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3))&
            *(scmax_short_pair-scmin_short_pair)
!!
!! scale stress components, CHECK IF THIS IS RIGHT!!!
!!          strs_pair(:,:,i3,i2)=strs_pair(:,:,i3,i2)/&
!!             (maxvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)&
!!             -minvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3))
!! new for variable scaling range
          strs_pair(:,:,i3,i2)=strs_pair(:,:,i3,i2)/&
            (maxvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)&
            -minvalue_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)),i3))&
            *(scmax_short_pair-scmin_short_pair)
        enddo
      enddo

!!
      return
      end
