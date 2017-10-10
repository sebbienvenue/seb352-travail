!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! this subroutine is derived from getonedeshortdw.f90
!! => output array is still called dedw_local
!!
!! called by:
!! - getdfshortdw.f90
!! - getdfshortdw_para.f90
!! - getdfpairdw.f90
!! - getdfpairdw_para.f90
!!
      subroutine getonededb(&
          maxnum_weights_local,&
          maxnodes_local,maxnum_layers_local,num_layers_local,&
          windexlocal,nodes_short_local,&
          weights_local,dnodes_values,dedw_local)
!!
      use fileunits
!!
      implicit none
!!
      integer maxnum_weights_local                                     ! in
      integer maxnodes_local                                           ! in
      integer maxnum_layers_local                                      ! in
      integer num_layers_local                                         ! in
      integer windexlocal(2*maxnum_layers_local)                       ! in
      integer nodes_short_local(0:maxnum_layers_local)                 ! in
      integer icount                                                   ! internal
      integer jcount                                                   ! internal
      integer kcount                                                   ! internal
      integer i1,i2,i3                                                 ! internal
!!
!! CAUTION: just one NN output node is assumed here => second index of dedw_local is 1
      real*8 dedw_local(maxnum_weights_local,1)                        ! out
      real*8 weights_local(maxnum_weights_local)                       ! in
      real*8 dnodes_values(maxnum_layers_local,maxnodes_local)         ! in
!!
!!
!! initializations
      icount=0
      jcount=0
      kcount=0
!!
!! get the bias weight derivatives
!!--------------------------------
!! output layer
      icount=windexlocal(2*num_layers_local) ! jump to correct weight
      do i2=1,nodes_short_local(num_layers_local) ! loop over noddes in output layer, typically only one
        dedw_local(icount,i2)=dnodes_values(num_layers_local,i2) ! ok
        icount=icount+1
      enddo
!! hidden layers
      if(num_layers_local.gt.1)then ! if hidden layers exist at all
        do i1=num_layers_local-1,1,-1 ! loop backwards over all hidden layers
          icount=windexlocal(2*i1) ! jump to correct starting point for bias weights in array 
!! set pointer to connecting weight a_{i2 i3}^{i1 i1+1}
          kcount=windexlocal(2*i1+1) ! ok 
          do i2=1,nodes_short_local(i1) ! loop over all nodes in hidden layer i1
            jcount=windexlocal(2*(i1+1)) ! set pointer to first bias weight in the subsequent layer
!! loop over all nodes in the subsequent layer
            do i3=1,nodes_short_local(i1+1)
              dedw_local(icount,1)=dedw_local(icount,1)&
                + dedw_local(jcount,1)*weights_local(kcount)*dnodes_values(i1,i2)
              jcount=jcount+1 ! next derivative in layer i1+1
              kcount=kcount+1
            enddo
            icount=icount+1 ! next bias weight in layer i1 (refering to node i2)
          enddo
        enddo ! i1
      endif ! num_layers_local .gt. 1
!!
      return
      end
