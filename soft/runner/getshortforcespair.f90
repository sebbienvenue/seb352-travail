!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - getallshortforcespair.f90
!!
        subroutine getshortforcespair(&
          num_atoms,num_pairs,zelemp,&
          symfunctionp,dsfuncdxyz_pair,&
          depairdsfunc,nnshortforce)
!!
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer num_atoms                                                 ! in
      integer num_pairs                                                 ! in
      integer ipairs                                                    ! internal
      integer zelemp(2,max_num_pairs)
      integer i1,i2,i3,i4,i5                                            ! internal
      integer icount                                                    ! internal
!!
      real*8 dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3)  ! in

      real*8 nnshortforce(3,max_num_atoms)                                ! out 
      real*8 depairdsfunc(max_num_pairs,maxnum_funcvalues_short_pair)               ! out (for stress) 

      real*8 weightsp(maxnum_weights_short_pair)
      real*8 symfunctionp(maxnum_funcvalues_short_pair,max_num_pairs)
      real*8 symfunction_pair(maxnum_funcvalues_short_pair)                           ! internal

!! CAUTION: just one output node is assumed here
       real*8 nnoutput                                                   ! internal 
       real*8 nodes_values(maxnum_layers_short_pair,maxnodes_short_pair)              ! internal
       real*8 nodes_sum(maxnum_layers_short_pair,maxnodes_short_pair)                 ! internal
       real*8 dnodes_values(maxnum_layers_short_pair,maxnodes_short_pair)             ! internal
       real*8 tempderivativep(maxnum_layers_short_pair,maxnodes_short_pair,maxnum_funcvalues_short_pair) ! internal

       real*8 alphagaussian
!!   
!!
!! initialization
      nnshortforce(:,:)   = 0.0d0
      depairdsfunc(:,:) = 0.0d0
      alphagaussian     = 0.5d0
!!
!! get deshortdsfunc for each atom
!!--------------------------------
      do i1=1,num_pairs 
!! don't calculate forces that we don't need
!    if(lupdatebyelement.and.(zelem(i1).ne.elemupdate)) goto 98
!!
        symfunction_pair(:)= symfunctionp(:,i1)
        weightsp(:)        = weights_short_pair(:,pairindex(zelemp(1,i1),zelemp(2,i1)))
        nodes_sum(:,:)     = 0.0d0 ! initialization
        nodes_values(:,:)  = 0.0d0 ! initialization
        dnodes_values(:,:) = 0.0d0 ! initialization
!!
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!! This is done on-the-fly here although it is also calculated before
!! in order to save memory
        call calconenn(1,maxnum_funcvalues_short_pair,maxnodes_short_pair,&
          maxnum_layers_short_pair,&
          num_layers_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1))),&
          maxnum_weights_short_pair,&
          nodes_short_pair(0,pairindex(zelemp(1,i1),zelemp(2,i1))),symfunction_pair,&
          weightsp,nodes_values,nodes_sum,nnoutput,&
          actfunc_short_pair(1,1,pairindex(zelemp(1,i1),zelemp(2,i1))))
!!
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
        ipairs=pairindex(zelemp(1,i1),zelemp(2,i1))
!!
        call getdnodes_values(maxnum_layers_short_pair,num_layers_short_pair,maxnodes_short_pair,&
          nodes_short_pair,ipairs,npairs,nodes_sum,nodes_values,dnodes_values,actfunc_short_pair)
!!
!! calculate the full derivative of E_i with respect to G_j^i
!!-----------------------------------------------------------
!!
        tempderivativep(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_short_atomic(1) values for each input node

        do i2=1,nodes_short_pair(1,pairindex(zelemp(1,i1),zelemp(2,i1)))              ! over all nodes in layer 1 ("target layer") 
          do i3=1,num_funcvalues_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1)))         ! over all nodes in previous layer
            icount=(i3-1)*nodes_short_pair(1,pairindex(zelemp(1,i1),zelemp(2,i1)))+i2 ! set pointer in weights array, don't need windex_short_atomic for first weights
            tempderivativep(1,i2,i3)=dnodes_values(1,i2)*weightsp(icount)
          enddo ! i3
        enddo ! i2

!! for layers 2 and beyond (if present) 
        if(num_layers_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1))).gt.1)then
          do i2=2,num_layers_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1)))  ! over all hidden and output layers
            do i3=1,nodes_short_pair(i2,pairindex(zelemp(1,i1),zelemp(2,i1)))       ! over all nodes in the target layer

              do i5=1, num_funcvalues_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1)))
!! we have to sum over the nodes in the previous layer (i4)
                do i4=1,nodes_short_pair(i2-1,pairindex(zelemp(1,i1),zelemp(2,i1))) ! sum over all nodes in previous layer

                  icount=windex_short_pair(2*i2-1, pairindex(zelemp(1,i1),zelemp(2,i1))) +&
                     (i4-1)*nodes_short_pair(i2,pairindex(zelemp(1,i1),zelemp(2,i1)))+i3-1 ! set pointer in weight array
                  tempderivativep(i2,i3,i5)=tempderivativep(i2,i3,i5) + &
                     dnodes_values(i2,i3)*weightsp(icount)*tempderivativep(i2-1,i4,i5)

                enddo ! i4
              enddo ! i5
            enddo ! i3
          enddo ! i2
        endif
!!
        depairdsfunc(i1,:)=tempderivativep(num_layers_short_pair(pairindex(zelemp(1,i1),zelemp(2,i1)))  ,1,:)
!!
!! debug
!        do i2=1,maxnum_funcvaluesp
!         write(ounit,*)i1,' DEBUG depairdsfunc ',i2,depairdsfunc(i1,i2)
!        enddo
!!
! 98     continue
      enddo ! i1 ! num_pairs
!!
!!
!! final force calculation
!!------------------------
!! summation over i2 and i4

      do i1=1,num_atoms  ! over all atoms in structure
        do i3=1,3 ! x,y,z
          do i2=1,num_pairs  ! over all pairs in structure
            do i4=1,num_funcvalues_short_pair(pairindex(zelemp(1,i2),zelemp(2,i2)))  ! over all symmetry functions
               nnshortforce(i3,i1)=nnshortforce(i3,i1)-depairdsfunc(i2,i4)*dsfuncdxyz_pair(i4,i2,i1,i3)
            enddo ! i4
          enddo ! i2
       enddo ! i3
!!99     continue
     enddo ! i1
!!
     return
     end
