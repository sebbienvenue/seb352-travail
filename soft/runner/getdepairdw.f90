!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - optimize_short_combinedpair.f90
!!
!! called once for each point
!!
      subroutine getdepairdw(&
           num_weightspairfree,&
           zelemp,num_pairs,wconstraintpidx,&
           symfunctionp,depairdw)
!!
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer num_weightspairfree(npairs)                                   ! in
      integer zelemp(2,max_num_pairs)                                       ! in
      integer num_pairs                                                     ! in
      integer wconstraintpidx(maxnum_weights_short_pair,npairs)                    ! in
      integer numelementpairs(npairs)                                       ! internal
      integer i1,i2                                                         ! internal
!!
      real*8 weightsp(maxnum_weights_short_pair)                                   ! internal
!! CAUTION: just one output node is assumed here
      real*8 depairdw(maxnum_weights_short_pair,1,npairs)                          ! out
!! CAUTION: just one output node is assumed here
      real*8 dedw(maxnum_weights_short_pair,1)                                     ! internal
      real*8 symfunction_pair(maxnum_funcvalues_short_pair)                           ! internal
      real*8 symfunctionp(maxnum_funcvalues_short_pair,max_num_pairs)                 ! in
!!
!!
!! initializations
      depairdw(:,:,:)   = 0.0d0
      numelementpairs(:)= 0 ! counts the number of pairs of each elemental combinations
!!
!! loop over all independant pairs of the structure
      do i1=1,num_pairs  
        dedw(:,:)          =0.0d0
        weightsp(:)        =weights_short_pair(:,pairindex(zelemp(1,i1),zelemp(2,i1))) 
        symfunction_pair(:)=symfunctionp(:,i1)
        numelementpairs(pairindex(zelemp(1,i1),zelemp(2,i1))) &
          = numelementpairs(pairindex(zelemp(1,i1),zelemp(2,i1)))+1 

!! calculate the derivative dedw for one specific pair 
!! even for weight constraints we calculate all dedw for simplicity

        call getonededw(1,&
          maxnum_funcvalues_short_pair,maxnum_weights_short_pair,&
          maxnodes_short_pair,maxnum_layers_short_pair,&
          num_layers_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1))),&
          windex_short_pair(1,pairindex(zelemp(1,i1),zelemp(2,i1))),nodes_short_pair(0,pairindex(zelemp(1,i1),zelemp(2,i1))),&
          symfunction_pair,weightsp,dedw,&
          actfunc_short_pair(1,1,pairindex(zelemp(1,i1),zelemp(2,i1))))

!! sum up the total derivative array for each pair
        do i2=1,num_weightspairfree(pairindex(zelemp(1,i1),zelemp(2,i1)))
          depairdw(i2,:,pairindex(zelemp(1,i1),zelemp(2,i1)))=& 
            depairdw(i2,:,pairindex(zelemp(1,i1),zelemp(2,i1))) &
            +dedw(wconstraintpidx(i2,pairindex(zelemp(1,i1),zelemp(2,i1))),:)
        enddo ! i2
      enddo ! i1
!!
!! normalization of the derivatives
      do i1=1,npairs
        if(numelementpairs(i1).gt.0)then
         depairdw(:,:,i1)=depairdw(:,:,i1)/dble(numelementpairs(i1))
!! Fit is not stable if we normalize by num_pairs:
!!          depairdw(:,:,i1)=depairdw(:,:,i1)/dble(num_pairs)
        endif
      enddo
!!
      return
      end
