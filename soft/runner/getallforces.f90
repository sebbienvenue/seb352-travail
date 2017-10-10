!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose:
!! calculate - nnshortforce
!!           - nnelecforce
!!           - nnelec
!! for npoints structures

!! called by: 
!!
      subroutine getallforces(ndim,npoints,&
       num_atoms_local,zelem_local,&
       lattice_local,xyzstruct_local,&
       symfunction_local,&
       minvalue_short_atomic,maxvalue_short_atomic,&
       nnshortforce_local,nnelecforce_local,&
       nnelec_local,&
       lperiodic_local)
!!
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use symfunctions
      use nnshort_atomic
!!
      implicit none
!!
!! input
      integer ndim                                                      ! in
      integer npoints                                                   ! in
      integer num_atoms_local(ndim)                                     ! in
      integer zelem_local(ndim,max_num_atoms)                           ! in
      integer i1,i2                                                     ! internal
!! internal
      integer zelem(max_num_atoms)                                      ! internal
      integer num_atoms                                                 ! internal
      integer, allocatable :: lsta(:,:)                                 ! internal, numbers of neighbors
      integer, allocatable :: lstc(:)                                   ! internal, identification of atom
      integer, allocatable :: lste(:)                                   ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_short_atomic(:)             ! internal
      integer, allocatable :: num_neighbors_elec(:)                     ! internal
      integer max_num_neighbors_short_atomic                            ! internal
      integer max_num_neighbors_elec                                    ! internal
      integer, allocatable :: neighboridx_short_atomic(:,:)             ! internal
      integer, allocatable :: invneighboridx_short_atomic(:,:)          ! internal
      integer, allocatable :: neighboridx_elec(:,:)                     ! internal
      integer, allocatable :: invneighboridx_elec(:,:)                  ! internal
      integer, allocatable :: atomindex_dummy(:)                        ! internal
!!
!! input
      real*8 symfunction_local(maxnum_funcvalues_short_atomic,max_num_atoms,ndim) ! in
      real*8 lattice_local(3,3,ndim)                                    ! in
      real*8 xyzstruct_local(3,max_num_atoms,ndim)                      ! in
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)! in
!! output
      real*8 nnshortforce_local(3,max_num_atoms,ndim)                   ! out 
      real*8 nnelecforce_local(3,max_num_atoms,ndim)                    ! out 
      real*8 nnelec_local(ndim)                                         ! out 
!! internal arrays
      real*8, allocatable :: lstb(:,:)                                  ! xyz and r_ij 
      real*8, allocatable :: dsfuncdxyz_short_atomic(:,:,:,:)           ! internal   
      real*8 nnelec                                                     ! internal 
      real*8 nnshortforce(3,max_num_atoms)                              ! internal 
      real*8 nnelecforce(3,max_num_atoms)                               ! internal 
      real*8 shortforce(3,max_num_atoms)                                ! internal
      real*8 strs(3,3,maxnum_funcvalues_short_atomic,max_num_atoms)     ! internal dummy
      real*8 symfunctiondummy(maxnum_funcvalues_short_atomic,max_num_atoms)! internal
!!
      logical lperiodic_local(ndim)                                     ! in
      logical lperiodic                                                 ! internal
      logical ldummy                                                    ! internal
!!
!!=======================================================
!! initializations 
!!=======================================================
      shortforce(:,:)=0.0d0 ! just dummy here
      nnshortforce_local(:,:,:)=0.0d0 
!!
!!=======================================================
!! loop over all structures 
!!=======================================================
      do i1=1,npoints
        num_atoms        =num_atoms_local(i1)
        zelem(:)         =zelem_local(i1,:)
        lperiodic        =lperiodic_local(i1)
!!
!!=======================================================
!! allocate neighbor list arrays for atomic short range part 
!!=======================================================
        allocate(lsta(2,max_num_atoms))
        allocate(lstc(listdim))
        allocate(lste(listdim))
        allocate(lstb(listdim,4))
        allocate(num_neighbors_short_atomic(num_atoms))
!!
!!=======================================================
!! determine max_num_neighbors_short_atomic and num_neighbors_short_tomic
!!=======================================================
        call getneighborsatomic(&
          num_atoms,num_neighbors_short_atomic,zelem,max_num_neighbors_short_atomic,&
          lsta,lstc,lste,&
          maxcutoff_short_atomic,lattice_local(1,1,i1),xyzstruct_local(1,1,i1),&
          lstb,lperiodic)
!!
!!=======================================================
!! allocate arrays and determine neighboridx and invneighboridx 
!!=======================================================
        allocate(neighboridx_short_atomic(num_atoms,0:max_num_neighbors_short_atomic))  
        allocate(invneighboridx_short_atomic(num_atoms,max_num_atoms))  
        allocate(dsfuncdxyz_short_atomic(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_short_atomic,3))
        call getneighboridxatomic(num_atoms,listdim,&
          max_num_atoms,max_num_neighbors_short_atomic,&
          lsta,lstc,neighboridx_short_atomic,invneighboridx_short_atomic)
!!
!!=======================================================
!! get dsfuncdxyz_short_atomic:
!! in principle dsfuncdxyz_short_atomic could be precalculated, but needs too much storage on disk
!! Caution: outgoing symmetry functions cannot be used because here they are not scaled
!!=======================================================
        allocate(atomindex_dummy(num_atoms))
        do i2=1,num_atoms
          atomindex_dummy(i2)=i2
        enddo
        call calconefunction_atomic(cutoff_type,max_num_neighbors_short_atomic,&
           max_num_atoms,1,num_atoms,atomindex_dummy,max_num_atoms,elementindex,&
           maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
           nelem,zelem,listdim,&
           lsta,lstc,lste,invneighboridx_short_atomic,&
           function_type_short_atomic,symelement_short_atomic,&
           xyzstruct_local(1,1,i1),symfunctiondummy,0.d0,&
           funccutoff_short_atomic,eta_short_atomic,rshift_short_atomic,&
           lambda_short_atomic,zeta_short_atomic,dsfuncdxyz_short_atomic,strs,lstb,&
           lperiodic,.true.,.false.,ldummy)
        deallocate(atomindex_dummy)
!!
!!=======================================================
!! deallocate neighbor list arrays for short range atomic case 
!!=======================================================
        deallocate(lsta)
        deallocate(lstc)
        deallocate(lste)
        deallocate(lstb)
!!
!!=======================================================
!! scale dsfuncdxyz_short_atomic
!!=======================================================
        if(lscalesym)then
          call scaledsfunc(max_num_neighbors_short_atomic,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            nelem,num_atoms,minvalue_short_atomic,maxvalue_short_atomic,&
            scmin_short_atomic,scmax_short_atomic,&
            zelem,dsfuncdxyz_short_atomic,strs)
        endif
!!
!!=======================================================
!! allocate neighbor list arrays for electrostatic part 
!!=======================================================
        allocate(lsta(2,max_num_atoms))
        allocate(lstc(listdim))
        allocate(lste(listdim))
        allocate(lstb(listdim,4))
        allocate(num_neighbors_elec(num_atoms))
!!
!!=======================================================
!! determine max_num_neighbors_elec and num_neighbors_elec
!!=======================================================
        call getneighborsatomic(&
          num_atoms,num_neighbors_elec,zelem,max_num_neighbors_elec,&
          lsta,lstc,lste,&
          maxcutoff_elec,lattice_local(1,1,i1),xyzstruct_local(1,1,i1),&
          lstb,lperiodic)
!!
!!=======================================================
!! allocate arrays and determine neighboridx_elec and invneighboridx_elec 
!!=======================================================
        allocate(neighboridx_elec(num_atoms,0:max_num_neighbors_elec))  
        allocate(invneighboridx_elec(num_atoms,max_num_atoms))
        call getneighboridxatomic(num_atoms,listdim,&
          max_num_atoms,max_num_neighbors_elec,&
          lsta,lstc,neighboridx_elec,invneighboridx_elec)
!!
!!=======================================================
!! deallocate neighbor list arrays for electrostatic case 
!!=======================================================
        deallocate(lsta)
        deallocate(lstc)
        deallocate(lste)
        deallocate(lstb)
!!
!!=======================================================
!! calculate nnshortforce, nnelecforce and nnelec 
!!=======================================================
!! FIXME this is not completed
        call getforces(max_num_neighbors_short_atomic,num_atoms,&
          num_neighbors_short_atomic,neighboridx_short_atomic,zelem,&
          max_num_neighbors_elec,num_neighbors_elec,neighboridx_elec,&
          symfunction_local(1,1,i1),&
          lattice_local(1,1,i1),xyzstruct_local(1,1,i1),&
          dsfuncdxyz_short_atomic,nnelec,nnshortforce,nnelecforce,lperiodic)
!!       
        nnshortforce_local(:,:,i1)=nnshortforce(:,:)
        if(lelec.and.(nn_type_elec.eq.2))then
          nnelec_local(i1)=nnelec/dble(num_atoms_local(i1))
          nnelecforce_local(:,:,i1)=nnelecforce(:,:)
        endif
!!
!!=======================================================
!! deallocate arrays 
!!=======================================================
        deallocate(dsfuncdxyz_short_atomic)
        deallocate(neighboridx_short_atomic)  
        deallocate(invneighboridx_short_atomic)  
        deallocate(num_neighbors_short_atomic)
        deallocate(neighboridx_elec)  
        deallocate(invneighboridx_elec)  
        deallocate(num_neighbors_elec)
!!
      enddo ! i1
!!
      return
      end
