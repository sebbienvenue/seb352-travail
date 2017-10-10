!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - geteshortpair.f90
!!
      subroutine calconeshortpair(num_pairs,&
        zelemp,symfunctionp,eshort,nnatomenergy)
!!
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer num_pairs
      integer zelemp(2,max_num_pairs)
      integer i1
!!
      real*8 symfunctionp(maxnum_funcvalues_short_pair,max_num_pairs)  ! in
      real*8 symfunction_pair(maxnum_funcvalues_short_pair)            ! internal
      real*8 weightsp(maxnum_weights_short_pair)                    ! internal
!! CAUTION: nnoutput assumes just one output node here
      real*8 nnoutput                                     ! internal
      real*8 eshort                                       ! out
      real*8 nodes_values_dummy(maxnum_layers_short_pair,maxnodes_short_pair)   ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_short_pair,maxnodes_short_pair)      ! just dummy in this routine
      real*8 nnatomenergy(max_num_pairs)                  ! out
!!
!!
      eshort=0.0d0
      nnatomenergy(:)=0.0d0
!!
!!#################################################################
!! serial original:
!!
      do i1=1,num_pairs 
        symfunction_pair(:) =  symfunctionp(:,i1)
        weightsp(:)         =  weights_short_pair(:,pairindex(zelemp(1,i1),zelemp(2,i1)))  

        call calconenn(1,maxnum_funcvalues_short_pair,maxnodes_short_pair,&
          maxnum_layers_short_pair,num_layers_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1))),&
          maxnum_weights_short_pair,nodes_short_pair(0,pairindex(zelemp(1,i1),zelemp(2,i1))),&
          symfunction_pair,weightsp,nodes_values_dummy,nodes_sum_dummy,nnoutput,&
          actfunc_short_pair(1,1,pairindex(zelemp(1,i1),zelemp(2,i1))))

        nnatomenergy(i1)=nnoutput
        eshort=eshort+nnoutput
      enddo ! i1
!!
      return
      end
