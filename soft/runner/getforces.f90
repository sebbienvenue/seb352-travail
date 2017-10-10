!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!!
      subroutine getforces(max_num_neighbors_short_atomic,num_atoms,&
        num_neighbors_short_atomic,neighboridx_short_atomic,zelem,&
        max_num_neighbors_elec,num_neighbors_elec,neighboridx_elec,&
        symfunction,&
        dsfuncdxyz_short_atomic,nnelec,nnshortforce,nnelecforce,lperiodic)
!!
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer max_num_neighbors_short_atomic                                 ! in
      integer max_num_neighbors_elec                                         ! in
      integer num_atoms                                                      ! in
      integer num_neighbors_short_atomic(num_atoms)                          ! in
      integer num_neighbors_elec(num_atoms)                                  ! in
      integer neighboridx_short_atomic(num_atoms,0:max_num_neighbors_short_atomic) ! in
      integer neighboridx_elec(num_atoms,0:max_num_neighbors_elec)           ! in
      integer ielem                                                          ! internal
      integer zelem(max_num_atoms)                                           ! in
      integer i1,i2,i3,i4,i5                                                 ! internal
      integer icount                                                         ! internal
      integer nndim                                                          ! internal
!!
      real*8 dsfuncdxyz_short_atomic(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_short_atomic,3) ! in
      real*8 nnshortforce(3,max_num_atoms)                                   ! out 
      real*8 nnelecforce(3,max_num_atoms)                                    ! out 
      real*8 nnelec                                                          ! out 
      real*8 weights(maxnum_weights_short_atomic)                            ! internal
      real*8 symfunction(maxnum_funcvalues_short_atomic,max_num_atoms)       ! in
      real*8 symfunction_atom(maxnum_funcvalues_short_atomic)                ! internal
      real*8 lattice(3,3)                                                    ! in
      real*8 xyzstruct(3,max_num_atoms)                                      ! in
!! CAUTION: just one output node is assumed here
      real*8 nodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic)  ! internal
      real*8 nodes_sum(maxnum_layers_short_atomic,maxnodes_short_atomic)     ! internal
      real*8 dnodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic) ! internal
      real*8 tempderivative(maxnum_layers_short_atomic,maxnodes_short_atomic,maxnum_funcvalues_short_atomic) ! internal
      real*8 alphagaussian
      real*8, dimension(:), allocatable :: nnoutput                          ! internal
      real*8, dimension(:,:), allocatable :: deshortdsfunc                   ! internal 
      real*8, dimension(:,:), allocatable :: dchargedsfunc                   ! internal 
      real*8, dimension(:,:,:), allocatable :: dchargedxyz                   ! internal 
      real*8 nnatomcharge(max_num_atoms)                                     ! internal
      real*8 distance                                                        ! internal
      real*8 invrij2                                                         ! internal
      real*8 drijdxyz(max_num_atoms,3)                                       ! internal
      real*8 fscreen                                                         ! internal
      real*8 fscreenderiv                                                    ! internal
!!   
      logical lperiodic                                                      ! in
!!
!!=======================================================
!! initializations
!!=======================================================
      nnshortforce(:,:) = 0.0d0
      nnelecforce(:,:)  = 0.0d0
      alphagaussian     = 0.5d0
      nnatomcharge(:)   = 0.0d0
      nndim             = 1
      if(lelec.and.(nn_type_elec.eq.2))then
        nndim=2
      endif
      allocate(nnoutput(nndim))
      allocate(deshortdsfunc(max_num_atoms,maxnum_funcvalues_short_atomic)) 
      deshortdsfunc(:,:)= 0.0d0
      if(lelec.and.(nn_type_elec.eq.2))then
        allocate(dchargedsfunc(max_num_atoms,maxnum_funcvalues_short_atomic)) 
        dchargedsfunc(:,:)= 0.0d0
        allocate(dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3))
        dchargedxyz(:,:,:)= 0.0d0
      endif
!!
!!=======================================================
!! get deshortdsfunc and dchargedsfunc for each atom
!!=======================================================
      do i1=1,num_atoms
!! don't calculate forces that we don't need
        if(lupdatebyelement.and.(zelem(i1).ne.elemupdate)) goto 98
!!
        symfunction_atom(:)=symfunction(:,i1)
        weights(:)=weights_short_atomic(:,elementindex(zelem(i1)))
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
!!
!!=======================================================
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!! This is done on-the-fly here although it is also calculated before
!! in order to save memory
!!=======================================================
        call calconenn(nndim,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic(elementindex(zelem(i1))),&
          maxnum_weights_short_atomic,nodes_short_atomic(0,elementindex(zelem(i1))),&
          symfunction_atom,weights,nodes_values,nodes_sum,&
          nnoutput,actfunc_short_atomic(1,1,elementindex(zelem(i1))))
!!
!!=======================================================
!! keep atomic charges for electrostatic energy and force evaluation
!!=======================================================
        if(lelec.and.(nn_type_elec.eq.2))then
          nnatomcharge(i1)=nnoutput(2)
        endif
!!
        if(luseforces)then
!!=======================================================
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
!!=======================================================
          ielem=elementindex(zelem(i1))
          call getdnodes_values(maxnum_layers_short_atomic,num_layers_short_atomic,maxnodes_short_atomic,&
            nodes_short_atomic,ielem,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_short_atomic)
!!
!!=======================================================
!! calculate the full derivative of E_i with respect to G_j^i
!!=======================================================
          tempderivative(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_short_atomic(1) values for each input node
          do i2=1,nodes_short_atomic(1,elementindex(zelem(i1))) ! over all nodes in layer 1 ("target layer") 
            do i3=1,num_funcvalues_short_atomic(elementindex(zelem(i1))) ! over all nodes in previous layer
              icount=(i3-1)*nodes_short_atomic(1,elementindex(zelem(i1)))+i2 ! set pointer in weights array, don't need windex_short_atomic for first weights
              tempderivative(1,i2,i3)=dnodes_values(1,i2)*weights(icount)
            enddo ! i3
          enddo ! i2
!! for layers 2 and beyond (if present) 
          if(num_layers_short_atomic(elementindex(zelem(i1))).gt.1)then
            do i2=2,num_layers_short_atomic(elementindex(zelem(i1))) ! over all hidden and output layers
              do i3=1,nodes_short_atomic(i2,elementindex(zelem(i1))) ! over all nodes in the target layer
                do i5=1,num_funcvalues_short_atomic(elementindex(zelem(i1)))
!! we have to sum over the nodes in the previous layer (i4)
                  do i4=1,nodes_short_atomic(i2-1,elementindex(zelem(i1))) ! sum over all nodes in previous layer
                    icount=windex_short_atomic(2*i2-1,elementindex(zelem(i1)))&
                      +(i4-1)*nodes_short_atomic(i2,elementindex(zelem(i1)))+i3-1 ! set pointer in weight array
                    tempderivative(i2,i3,i5)=tempderivative(i2,i3,i5)&
                      +dnodes_values(i2,i3)*weights(icount)*tempderivative(i2-1,i4,i5)
                  enddo ! i4
                enddo ! i5
              enddo ! i3
            enddo ! i2
          endif
!!
          deshortdsfunc(i1,:)=tempderivative(num_layers_short_atomic(elementindex(zelem(i1))),1,:)
          if(lelec.and.(nn_type_elec.eq.2))then
            dchargedsfunc(i1,:)=tempderivative(num_layers_short_atomic(elementindex(zelem(i1))),2,:)
          endif
!!
        endif ! luseforces
!!
 98     continue
      enddo ! i1 ! num_atoms
!!
!!=======================================================
!! final short range force calculation
!!=======================================================
      do i2=1,num_atoms ! over all atoms in structure
        do i1=0,num_neighbors_short_atomic(i2) ! JB CAUTION 
!! don't calculate forces that we don't need
        if(lupdatebyelement.and.(zelem(neighboridx_short_atomic(i2,i1)).ne.elemupdate)) goto 99
          do i3=1,3 ! x,y,z
            do i4=1,num_funcvalues_short_atomic(elementindex(zelem(i2))) ! over all symmetry functions
              nnshortforce(i3,neighboridx_short_atomic(i2,i1))=nnshortforce(i3,neighboridx_short_atomic(i2,i1)) &
                -deshortdsfunc(i2,i4)*dsfuncdxyz_short_atomic(i4,i2,i1,i3)
!! now get dchargedxyz
              if(lelec.and.(nn_type_elec.eq.2))then
                dchargedxyz(i2,i1,i3)=dchargedxyz(i2,i1,i3) &
                + dchargedsfunc(i2,i4)*dsfuncdxyz_short_atomic(i4,i2,i1,i3) 
              endif
            enddo ! i4
          enddo ! i3
 99     continue
        enddo ! i1
      enddo ! i2
!!
!!=======================================================
!! calculation of electrostatic energies and forces 
!!=======================================================
      if(lelec.and.(nn_type_elec.eq.2))then
        if(lperiodic)then
          call getewaldenergy(max_num_neighbors_elec,&
            neighboridx_elec,num_neighbors_elec,num_atoms,zelem,&
            lattice,xyzstruct,&
            nnatomcharge,nnelec,&
            dchargedxyz,nnelecforce,.false.)
        else ! not periodic
          call electrostatic(num_atoms,&
            nnatomcharge,xyzstruct,nnelec)
          if(luseforces)then
            nnelecforce(:,:) =0.0d0 ! important!!! 
            if(num_atoms.gt.1) then
              do i3=0,num_neighbors_elec(i1) !! FIXME
                do i4=1,3
                  do i1=1,num_atoms
                    do i2=i1+1,num_atoms
                      distance=(xyzstruct(1,i1)-xyzstruct(1,i2))**2 + &
                        (xyzstruct(2,i1)-xyzstruct(2,i2))**2 + &
                        (xyzstruct(3,i1)-xyzstruct(3,i2))**2
                      invrij2 =1.d0/distance
                      distance=dsqrt(distance)
!! calculation of \frac{\partial r_{ij}}{\partial \alpha}
                      drijdxyz(i3,i4)=0.0d0 ! initialization
                      if(i1.eq.i3)then
                        drijdxyz(i3,i4)=(xyzstruct(i4,i1)-xyzstruct(i4,i2))/distance
                      elseif(i2.eq.i3)then
                        drijdxyz(i3,i4)=-1.d0*(xyzstruct(i4,i1)-xyzstruct(i4,i2))/distance
                      endif
                      if(lscreen) then
                        call getscreenfunctionforelectrostatics(&
                          distance,fscreen,&
                          fscreenderiv,drijdxyz(i3,i4))
                        nnelecforce(i4,i3)=nnelecforce(i4,i3) &
                          - invrij2*((dchargedxyz(i1,i3,i4)*nnatomcharge(i2) &
                          + nnatomcharge(i1)*dchargedxyz(i2,neighboridx_elec(i2,i3),i4))*&
                          distance*fscreen - nnatomcharge(i1)* &
                          nnatomcharge(i2)*drijdxyz(i3,i4)* &
                          fscreen+fscreenderiv*nnatomcharge(i1)* &
                          nnatomcharge(i2)*distance)
                      else
                        nnelecforce(i4,i3)=nnelecforce(i4,i3) &
                          - invrij2*((dchargedxyz(i1,neighboridx_elec(i1,i3),i4)*nnatomcharge(i2) &
                          + nnatomcharge(i1)*dchargedxyz(i2,neighboridx_elec(i2,i3),i4))*distance&
                          -nnatomcharge(i1)*nnatomcharge(i2)*drijdxyz(i3,i4))
                      endif ! lscreen
                    enddo ! i2
                  enddo ! i1
                enddo ! i4
              enddo ! i3
            endif ! num_atoms.gt.1
!!
          endif ! luseforces
        endif ! lperiodic
      endif ! lelec
!!
      deallocate(nnoutput)
      deallocate(deshortdsfunc) 
      if(lelec.and.(nn_type_elec.eq.2))then
        deallocate(dchargedsfunc) 
        deallocate(dchargedxyz)
      endif
!!
      return
      end
