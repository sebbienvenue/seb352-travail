!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - predictionpair.f90
!! - optimize_short_combinedpair.f90
!!
      subroutine calconeshort_parapair(istruct,num_pairs_para,pindex,&
        zelemp,symfunctionp,nnpairenergy,nntotalenergy)
!!
      use mpi_mod
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer num_pairs_para                                        ! in
      integer pindex(num_pairs_para)                                ! in
      integer istruct                                           ! in
      integer i1,i2,i3
      integer zelemp(2,max_num_pairs)
!!
      real*8 symfunctionp(maxnum_funcvalues_short_pair,num_pairs_para)        ! in
      real*8 symfunction_pair(maxnum_funcvalues_short_pair)               ! internal
      real*8 weightsp(maxnum_weights_short_pair)                       ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                            ! internal
      real*8 nodes_values_dummy(maxnum_layers_short_pair,maxnodes_short_pair) ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_short_pair,maxnodes_short_pair)    ! just dummy in this routine
      real*8 nnpairenergy(max_num_pairs)                         ! out                                        
      real*8 nntotalenergy                                       ! out
!!
      do i1=1,num_pairs_para     
!!
        symfunction_pair(:)=symfunctionp(:,i1)
        weightsp(:)=weights_short_pair(:,pairindex(zelemp(1,i1),zelemp(2,i1))) 
!!
!! calculate nnoutput for pair i1
       call calconenn(1,maxnum_funcvalues_short_pair,maxnodes_short_pair,&
           maxnum_layers_short_pair,num_layers_short_pair(pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1)))),&
           maxnum_weights_short_pair,nodes_short_pair(0,pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1)))),&
           symfunction_pair,weightsp,nodes_values_dummy,nodes_sum_dummy,nnoutput,&
           actfunc_short_pair(1,1,pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1)))))
!!
        if(ldebug)then
!!        write(ounit,'(a,i6,a)')'nodes values for atom ',i1,' :'
          write(ounit,'(a)')'            struct  atom layer  node           nodes_sum        nodes_values'

          do i2=1,num_layers_short_pair(pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1))))
            do i3=1,nodes_short_pair(i2,pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1)))) 
              write(ounit,'(a12,4i6,2f20.10)')'NODES SHORT ',istruct,i1,i2,i3,&
                                  nodes_sum_dummy(i2,i3),nodes_values_dummy(i2,i3)
            enddo ! i3
          enddo ! i2
          write(ounit,*)'-------------------------------------------------------------'
        endif ! ldebug
!!
        nnpairenergy(pindex(i1))=nnoutput
        nntotalenergy = nntotalenergy+nnoutput
!!
      enddo ! i1
!!
      return
      end
