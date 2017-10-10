!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getdfshortdw.f90
!! - getdfshortdw_para.f90
!! - getdfpairdw.f90
!! - getdfpairdw_para.f90
!!
      subroutine get_dnodes_values_dg(&
          maxnum_layers_local,num_layers_local,maxnodes_local,&
          maxnum_funcvalues_local,&
          maxnum_weights_local,windex_local,nodes_local,&
          weights_local,dnodes_values_local,dnodes_values_dg_local)
!!
      use fileunits
!!
      implicit none
!!      
      integer i1,i2,i3,i4                                                   ! internal
      integer maxnum_layers_local                                           ! in
      integer num_layers_local                                              ! in
      integer maxnodes_local                                                ! in
      integer maxnum_funcvalues_local                                       ! in
      integer maxnum_weights_local                                          ! in
      integer windex_local(2*maxnum_layers_local)                           ! in
      integer nodes_local(0:maxnum_layers_local)                            ! in
      integer icount                                                        ! internal
!!
      real*8 weights_local(maxnum_weights_local)                            ! in
      real*8 dnodes_values_dg_local(maxnum_layers_local,maxnodes_local,maxnum_funcvalues_local)! out 
      real*8 dnodes_values_local(maxnum_layers_local,maxnodes_local)        ! in
!!
!!
!!
!! calculate dnodes_values_dg recursively starting in hidden layer 1
      dnodes_values_dg_local(:,:,:)=0.0d0
!!
      do i1=1,num_layers_local ! target layer
        do i2=1,nodes_local(i1) ! node in target layer
          do i3=1,nodes_local(0) ! input node
            if(i1.eq.1) then ! 1st hidden layer
              icount=windex_local(1)+nodes_local(i1)*(i3-1)+i2-1
              dnodes_values_dg_local(i1,i2,i3)=weights_local(icount)
            else ! i1.gt.1
              do i4=1,nodes_local(i1-1) ! sum over all nodes in the previous hidden layer
                icount=windex_local(2*i1-1)+nodes_local(i1)*(i4-1)+i2-1      
                dnodes_values_dg_local(i1,i2,i3)=dnodes_values_dg_local(i1,i2,i3)&
                  +weights_local(icount)*dnodes_values_local(i1-1,i4)*dnodes_values_dg_local(i1-1,i4,i3)
              enddo
            endif
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
      return
      end
