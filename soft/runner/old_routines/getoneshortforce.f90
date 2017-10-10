!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! this subroutine is derived from routine getshortforces.f90
!!
!! called by: - optimize_short_combined.f90
!!
      subroutine getoneshortforce(iatom,ixyz,&
        maxnum_funcvalues,num_funcvalues,&
        num_atoms,maxnum_weightsshort,&
        nelem,maxnum_layersshort,num_layersshort,nodes_short,maxnodes_short,&
        zelem,windex,&
        weights_short,symfunction,&
        dsfuncdxyz,nnxyzforce,actfunc_short)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer maxnum_funcvalues                                         ! in
      integer num_funcvalues(nelem)                                     ! in
      integer num_atoms                                                 ! in
      integer maxnum_weightsshort                                       ! in
      integer nelem                                                     ! in
      integer ielem                                                     ! internal
      integer maxnum_layersshort                                        ! in
      integer num_layersshort(nelem)                                    ! in
      integer nodes_short(0:maxnum_layersshort,nelem)                            ! in
      integer maxnodes_short                                            ! in
      integer zelem(max_num_atoms)                                      ! in
      integer windex(2*maxnum_layersshort,nelem)                                 ! in
      integer iatom                                                     ! in
      integer ixyz                                                      ! in
      integer i1,i2,i3,i4,i5                                            ! internal
      integer icount                                             ! internal
!!
      real*8 dsfuncdxyz(maxnum_funcvalues,max_num_atoms,max_num_atoms,3)   ! in
      real*8 nnxyzforce(3,max_num_atoms)                                ! out 
      real*8 deshortdsfunc(max_num_atoms,maxnum_funcvalues)                ! out (for stress) 
      real*8 weights_short(maxnum_weightsshort,nelem)                   ! in
      real*8 weights(maxnum_weightsshort)                               ! internal
      real*8 symfunction(maxnum_funcvalues,max_num_atoms)                  ! in
      real*8 symfunction_atom(maxnum_funcvalues)                           ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                     ! internal 
      real*8 nodes_values(maxnum_layersshort,maxnodes_short)            ! internal
      real*8 nodes_sum(maxnum_layersshort,maxnodes_short)               ! internal
      real*8 dnodes_values(maxnum_layersshort,maxnodes_short)           ! internal
      real*8 tempderivative(maxnum_layersshort,maxnodes_short,maxnum_funcvalues) ! internal
      real*8 alphagaussian
!!   
      character*1 actfunc_short(maxnodes_short,maxnum_layersshort,nelem)                        ! in
!!
!!      write(ounit,*)'serial getoneshortforce starts'
!!
!! initialization
      nnxyzforce(:,:)   = 0.0d0
      deshortdsfunc(:,:)= 0.0d0
      alphagaussian     = 0.5d0
!!
!! get deshortdsfunc for each atom
!!--------------------------------
      do i1=1,num_atoms
        symfunction_atom(:)=symfunction(:,i1)
        weights(:)=weights_short(:,elementindex(zelem(i1)))
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
!!
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!! This is done on-the-fly here although it is also calculated before
!! in order to save memory
        call calconenn(maxnum_funcvalues,maxnodes_short,&
          maxnum_layersshort,num_layersshort(elementindex(zelem(i1))),&
          maxnum_weightsshort,&
          nodes_short(0,elementindex(zelem(i1))),&
          symfunction_atom,weights,nodes_values,nodes_sum,&
          nnoutput,actfunc_short(1,1,elementindex(zelem(i1))))
!!
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
        ielem=elementindex(zelem(i1))
!!
        call getdnodes_values(maxnum_layersshort,num_layersshort,maxnodes_short,&
          nodes_short,ielem,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_short)
!!
!! calculate the full derivative of E_i with respect to G_j^i
!!-----------------------------------------------------------
!!
       tempderivative(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_short(1) values for each input node
       do i2=1,nodes_short(1,elementindex(zelem(i1))) ! over all nodes in layer 1 ("target layer") 
         do i3=1,num_funcvalues(elementindex(zelem(i1))) ! over all nodes in previous layer
           icount=(i3-1)*nodes_short(1,elementindex(zelem(i1)))+i2 ! set pointer in weights array, don't need windex for first weights
           tempderivative(1,i2,i3)=dnodes_values(1,i2)*weights(icount)
         enddo ! i3
       enddo ! i2
!! for layers 2 and beyond (if present) 
      if(num_layersshort(elementindex(zelem(i1))).gt.1)then
        do i2=2,num_layersshort(elementindex(zelem(i1))) ! over all hidden and output layers
          do i3=1,nodes_short(i2,elementindex(zelem(i1))) ! over all nodes in the target layer
            do i5=1,num_funcvalues(elementindex(zelem(i1)))
!! we have to sum over the nodes in the previous layer (i4)
              do i4=1,nodes_short(i2-1,elementindex(zelem(i1))) ! sum over all nodes in previous layer
              icount=windex(2*i2-1,elementindex(zelem(i1)))+(i4-1)*nodes_short(i2,elementindex(zelem(i1)))+i3-1 ! set pointer in weight array
              tempderivative(i2,i3,i5)=tempderivative(i2,i3,i5) + &
              dnodes_values(i2,i3)*weights(icount)*tempderivative(i2-1,i4,i5)
              enddo ! i4
            enddo ! i5
          enddo ! i3
        enddo ! i2
      endif
!!
      deshortdsfunc(i1,:)=tempderivative(num_layersshort(elementindex(zelem(i1))),1,:)
!!
      enddo ! i1 ! num_atoms
!!
!!
!! final force calculation
!!------------------------
!! summation over i2 and i4
!!      do i1=1,num_atoms  ! over all atoms in structure
!!        do i2=1,num_atoms ! over all atoms in structure
!!          do i3=1,3 ! x,y,z
!!            do i4=1,num_funcvalues ! over all symmetry functions
!!              nnxyzforce(i3,i1)=nnxyzforce(i3,i1) &
!!                -deshortdsfunc(i2,i4)*dsfuncdxyz(i4,i2,i1,i3)
!!            enddo ! i4
!!          enddo ! i3
!!        enddo ! i2
!!      enddo ! i1
!!
!!
        do i2=1,num_atoms ! over all atoms in structure
          do i4=1,num_funcvalues(elementindex(zelem(i2))) ! over all symmetry functions
            nnxyzforce(ixyz,iatom)=nnxyzforce(ixyz,iatom) &
              -deshortdsfunc(i2,i4)*dsfuncdxyz(i4,i2,iatom,ixyz)
          enddo ! i4
        enddo ! i2
!!
      return
      end
