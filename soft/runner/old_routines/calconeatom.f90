!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconecharge.f90
!! - calconecharge_para.f90
!! - calconeshort.f90
!! - calconeshort_para.f90
!! - calconeshortpair.f90
!! - calconeshort_parapair.f90
!! - getalphalm.f90
!! - getcoulombdchargedxyz_para.f90
!! - getcoulombforces.f90
!! - getdchargedxyz.f90
!! - getdchargedxyz_para.f90
!! - getdfshortdw.f90
!! - getdfshortdw_para.f90
!! - getdfpairdw.f90
!! - getdfpairdw_para.f90
!! - getonedeshortdw.f90
!! - getoneshortforce.f90
!! - getoneshortforce_para.f90
!! - getshortforcespair.f90
!! - getoneshortforcepair_para.f90
!! - getshortforces.f90
!! - getshortforces_para.f90
!! - optimize_ewald.f90
!!
      subroutine calconeatom(maxnum_funcvalues_local,maxnodes,&
          maxnum_layers,num_layers,maxnum_weights,nodes_local,&
          symfunction_atom,weights,nodes_values_local,nodes_sum_local,&
          nnoutput,actfunc_local,ldebug)
!!
      use globaloptions
      use fileunits
!!
      implicit none
!!
      integer maxnum_funcvalues_local                                       ! in
      integer maxnum_layers                                                 ! in
      integer num_layers                                                    ! in
      integer maxnum_weights                                                ! in
      integer nodes_local(0:maxnum_layers)                                  ! in
      integer vecdim1 ! number of nodes in previous layer                   ! internal
      integer vecdim2 ! number of nodes in target layer                     ! internal
      integer maxnodes                                                      ! in
      integer i1,i2,i3                                                      ! internal
      integer icount                                                        ! internal
      integer jcount                                                        ! internal
!!
      real*8 symfunction_atom(maxnum_funcvalues_local)                      ! in 
      real*8 weights(maxnum_weights)                                        ! in
      real*8 one                                                            ! internal
      real*8 nodes_sum_local(maxnum_layers,maxnodes)                        ! out
      real*8 nodes_values_local(maxnum_layers,maxnodes)                     ! out
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                                       ! out
      real*8, dimension(:,:)  , allocatable :: biasweights                  ! internal
      real*8, dimension(:,:)  , allocatable :: connectionweights            ! internal
      real*8, dimension(:,:)  , allocatable :: nodevalues_in                ! internal
      real*8, dimension(:,:)  , allocatable :: nodevalues_out               ! internal
      real*8 alphagaussian                    
!!
      character*1 actfunc_local(maxnodes,maxnum_layers)                     ! in
!!
      logical ldebug                                                        ! in
!!
      one                    = 1.0d0
      icount                 = 0
      nodes_values_local(:,:)= 0.0d0
      nodes_sum_local(:,:)   = 0.0d0
      alphagaussian          = 0.5d0
!!
!! loop over all hidden and output layers
      do i1=1,num_layers
        vecdim1=nodes_local(i1-1) 
        vecdim2=nodes_local(i1)
!!
        allocate (nodevalues_in(1,vecdim1))
        allocate (nodevalues_out(1,vecdim2))
        allocate (biasweights(1,vecdim2))
        allocate (connectionweights(vecdim1,vecdim2))
!!
!! assign connecting weights
        do i3=1,vecdim1 ! over nodes in previous layers
          do i2=1,vecdim2 ! over target nodes
            icount=icount+1
            connectionweights(i3,i2)=weights(icount)
          enddo
        enddo
!! assign bias weights
        do i2=1,vecdim2
          icount=icount+1
          biasweights(1,i2)=weights(icount)
        enddo

!! assign the node values in the previous layer(s)
        jcount=0
        if(i1.eq.1)then ! i1-1 is input layer
          do i2=1,nodes_local(0)
            jcount=jcount+1
            nodevalues_in(1,jcount)=symfunction_atom(jcount)
          enddo ! i2
        else
          do i2=1,nodes_local(i1-1)
            jcount=jcount+1
            nodevalues_in(1,jcount)=nodes_values_local(i1-1,i2)
          enddo
        endif
!!
!!
!! calculate the target node values
!!       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!      C := alpha*op( A )*op( B ) + beta*C,
!!  M      - INTEGER.
!!           On entry,  M  specifies  the number  of rows  of the  matrix
!!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!!           Unchanged on exit.
!!  N      - INTEGER.
!!           On entry,  N  specifies the number  of columns of the matrix
!!           op( B ) and the number of columns of the matrix C. N must be
!!           at least zero.
!!           Unchanged on exit.
!!  K      - INTEGER.
!!           On entry,  K  specifies  the number of columns of the matrix
!!           op( A ) and the number of rows of the matrix op( B ). K must
!!           be at least  zero.
!!           Unchanged on exit.
        call dgemm('n','n',1,vecdim2,vecdim1,one,nodevalues_in,1,&
                connectionweights,vecdim1,one,biasweights,1)
!! after dgemm the array biasweights contains the sums at the nodes
!!
!! if requested normalize by the number of nodes in previous layer
        if(lnormnodes)then
          do i2=1,vecdim2 ! loop over all nodes in current layer
            biasweights(1,i2)=biasweights(1,i2)/dble(vecdim1)
          enddo
        endif
!!
!! store node values before application of activation function (needed for derivatives) 
        do i2=1,vecdim2
          nodes_sum_local(i1,i2)=biasweights(1,i2)
        enddo
!!
!! apply the activation functions
        do i2=1,vecdim2
          if(actfunc_local(i2,i1).eq."t")then
            biasweights(1,i2)=tanh(biasweights(1,i2))
          elseif(actfunc_local(i2,i1).eq."g")then
            biasweights(1,i2)=dexp(-alphagaussian*(biasweights(1,i2))**2)
          elseif(actfunc_local(i2,i1).eq."l")then
          elseif(actfunc_local(i2,i1).eq."c")then
            biasweights(1,i2)=dcos(biasweights(1,i2))
          elseif(actfunc_local(i2,i1).eq."s")then
            biasweights(1,i2)=1.d0/(1.d0+dexp(-1.d0*biasweights(1,i2)))
          elseif(actfunc_local(i2,i1).eq."S")then
            biasweights(1,i2)=1.d0-1.d0/(1.d0+dexp(-1.d0*biasweights(1,i2)))
          elseif(actfunc_local(i2,i1).eq."e")then
            biasweights(1,i2)=dexp(-1.d0*biasweights(1,i2))
          else
            write(ounit,*)"Error: Unknown activation function ",actfunc_local(i2,i1)
            stop
          endif
        enddo
!!
!! store final node values in nodes_values array
        do i2=1,vecdim2
          nodes_values_local(i1,i2)=biasweights(1,i2)
        enddo
!!
        deallocate (connectionweights)
        deallocate (biasweights)
        deallocate (nodevalues_in)
        deallocate (nodevalues_out)
!!
      enddo ! i1 loop over all hidden layers and output layers
!!
!! prepare output of the NN
      nnoutput=nodes_values_local(num_layers,1)
!!
      return
      end
