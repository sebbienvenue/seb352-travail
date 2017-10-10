!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!!
!! This routine is called once per atom and calculates the derivatives dedw_local 
!! of the short range energies for that atom
!!
      subroutine getonededw_atomic(nndim,&
          maxnum_funcvalues_local,maxnum_weights_local,&
          maxnodes_local,maxnum_layers_local,num_layers_local,&
          windex_local,nodes_local,&
          symfunction_atom,weights_local,dedw_local,&
          actfunc_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer maxnum_funcvalues_local                               ! in
      integer maxnum_weights_local                                  ! in
      integer maxnodes_local                                        ! in
      integer maxnum_layers_local                                   ! in
      integer num_layers_local                                      ! in
      integer windex_local(2*maxnum_layers_local)                   ! in
      integer nodes_local(0:maxnum_layers_local)                    ! in
      integer icount                                                ! internal
      integer jcount                                                ! internal
      integer kcount                                                ! internal
      integer i1,i2,i3,i4                                           ! internal
      integer nndim                                                 ! in
!!
      real*8 symfunction_atom(maxnum_funcvalues_local)              ! in
      real*8 dedw_local(maxnum_weights_local,nndim)                 ! out
      real*8 weights_local(maxnum_weights_local)                    ! in
      real*8 nnoutput                                               ! internal dummy
      real*8 nodes_values_local(maxnum_layers_local,maxnodes_local) ! internal
      real*8 nodes_sum_local(maxnum_layers_local,maxnodes_local)    ! internal
      real*8 dnodes_values_local(maxnum_layers_local,maxnodes_local)! internal
!!
      character*1 actfunc_local(maxnodes_local,maxnum_layers_local) ! in
!!
!!
!! initializations
      icount            = 0
      jcount            = 0
      kcount            = 0
      nodes_sum_local(:,:)    = 0.0d0
      nodes_values_local(:,:) = 0.0d0
      dnodes_values_local(:,:)= 0.0d0
!!
!! Calculation of the values on all nodes nodes_values_local in the NN 
!! (needed below for the derivatives dedw_local and dnodes_values_local).
!! This is done on-the-fly here to save memory, although it is also calculated before in the subroutine geteshort.
      call calconenn(1,maxnum_funcvalues_local,maxnodes_local,&
        maxnum_layers_local,num_layers_local,maxnum_weights_local,nodes_local,&
        symfunction_atom,weights_local,nodes_values_local,nodes_sum_local,&
        nnoutput,actfunc_local)
!!
!! chain rule: f(x)=u(v(x)) => f'(x)=u'(v(x))*v'(x)
!! here: y(w)=f_a(g(w)) => y'(w)=f_a'(g(w)) g'(w)
!!
!! Calculate first part of the derivative: f_a' at each node (\frac{\partial f_a(g(w))}{\partial w})
!! In order to use subroutine getdnodes_values_local here, we set ielem and nelem to 1. 
      call getdnodes_values(maxnum_layers_local,num_layers_local,maxnodes_local,&
        nodes_local,1,1,nodes_sum_local,nodes_values_local,dnodes_values_local,actfunc_local)
!!
!! get the bias weight derivatives
!!--------------------------------
!! output layer
      icount=windex_local(2*num_layers_local) ! jump to correct weight
      do i2=1,nodes_local(num_layers_local) ! loop over noddes in output layer
        dedw_local(icount,i2)=dnodes_values_local(num_layers_local,i2) ! ok
        icount=icount+1
      enddo
!! hidden layers
      if(num_layers_local.gt.1)then ! if hidden layers exist at all
        do i1=num_layers_local-1,1,-1 ! loop backwards over all hidden layers
          icount=windex_local(2*i1) ! jump to correct starting point for bias weights in array 
!! set pointer to connecting weight a_{i2 i3}^{i1 i1+1}
          kcount=windex_local(2*i1+1) ! ok 
          do i2=1,nodes_local(i1) ! loop over all nodes in hidden layer i1
            jcount=windex_local(2*(i1+1)) ! set pointer to first bias weight in the subsequent layer
!! loop over all nodes in the subsequent layer
            do i3=1,nodes_local(i1+1)
!! CHECK: is this generalisation to two output nodes ok?
              do i4=1,nodes_local(num_layers_local)
                dedw_local(icount,i4)=dedw_local(icount,i4)&
                  + dedw_local(jcount,i4)*weights_local(kcount)*dnodes_values_local(i1,i2)
                jcount=jcount+1 ! next derivative in layer i1+1
                kcount=kcount+1
              enddo ! i4
            enddo ! i3
            icount=icount+1 ! next bias weight in layer i1 (refering to node i2)
          enddo
        enddo ! i1
      endif ! num_layers_local .gt. 1
!!
!! the calculation of the bias weight derivatives is ok
!!
!! get the connection weight derivatives'
!!--------------------------------------
!! loop over all target layers
        do i1=1,num_layers_local
          icount=windex_local(2*i1-1) ! starting point in weights and dedw_local arrays
!! loop over all nodes in previous layer
          if(i1.gt.1)then ! previous layer is hidden layer 
            do i2=1,nodes_local(i1-1)
!! set the pointer for the bias weight derivative in the target layer
            jcount=windex_local(2*i1)
!! loop over all nodes in target layer
              do i3=1,nodes_local(i1)
                do i4=1,nodes_local(num_layers_local)
                  dedw_local(icount,i4)=nodes_values_local(i1-1,i2)*dedw_local(jcount,i4)
                  icount=icount+1
                  jcount=jcount+1
                enddo ! i4
              enddo ! i3
            enddo ! i2
          else ! i1.gt.1, previous layer is input layer
            do i2=1,nodes_local(0)
!! set the pointer for the bias weight derivative in the target layer
            jcount=windex_local(2*i1)
!! loop over all nodes in target layer
              do i3=1,nodes_local(i1)
                do i4=1,nodes_local(num_layers_local)
                  dedw_local(icount,i4)=symfunction_atom(i2)*dedw_local(jcount,i4)
                  icount=icount+1
                  jcount=jcount+1
                enddo ! i4
              enddo ! i3
            enddo ! i2
          endif
        enddo ! i1
!!
      return
      end
