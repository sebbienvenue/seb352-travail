!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - predictionpair.f90
!!
      subroutine scalesympair_para(num_pairs,&
         pindex,zelemp_list,symfunctionp,&
         minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)
!!
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer num_pairs
      integer zelemp_list(2,max_num_pairs)
      integer i2,i3
      integer pindex(max_num_pairs)
!!
!!      real*8 symfunctionp(maxnum_funcvaluesp,num_pairs_para)
      real*8 symfunctionp(maxnum_funcvalues_short_pair,max_num_pairs)
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)
!!
        do i2=1,num_pairs
          do i3=1,num_funcvalues_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2)))) 
            if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
              symfunctionp(i3,pindex(i2))=symfunctionp(i3,pindex(i2))&
                  - avvalue_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2))),i3) 
            elseif(lscalesym.and..not.lcentersym)then
!! Scale each symmetry function value for the respective element
              symfunctionp(i3,pindex(i2))=&
             (symfunctionp(i3,pindex(i2)) &
             -minvalue_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2))),i3))/ &
             (maxvalue_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2))),i3)-&
              minvalue_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2))),i3))&
              *(scmax_short_pair-scmin_short_pair) + scmin_short_pair
            elseif(lscalesym.and.lcentersym)then
              symfunctionp(i3,pindex(i2))=&
             (symfunctionp(i3,pindex(i2))-avvalue_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2))),i3))&
             / &
             (maxvalue_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2))),i3)- &
              minvalue_short_pair(pairindex(zelemp_list(1,pindex(i2)),zelemp_list(2,pindex(i2))),i3))
            else
            endif
          enddo ! i3
        enddo ! i2
!!
      return
      end
