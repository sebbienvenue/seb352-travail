!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate nnshortforce for npoints structures 

!! called by: 
!! - getshortenergies_para.f90
!! - optimize_short_combined.f90 
!! - optimize_atomic.f90   
!!
      subroutine getallshortforces(ndim,npoints,&
       num_atoms_local,zelem_local,&
       symfunction_local,nnshortforce_local,&
       lattice_local,xyzstruct_local,minvalue_short_atomic,maxvalue_short_atomic,&
       lperiodic_local)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
      use symfunctions
      use nnshort_atomic
!!
      implicit none
!!
      integer ndim                                                      ! in
      integer npoints                                                   ! in
      integer num_atoms                                                 ! internal
      integer num_atoms_local(ndim)                                     ! in
      integer zelem(max_num_atoms)                                      ! internal
      integer zelem_local(ndim,max_num_atoms)                           ! in
      integer i1,i2                                                     ! internal
      integer, allocatable :: lsta(:,:)                                 ! internal, numbers of neighbors
      integer, allocatable :: lstc(:)                                   ! internal, identification of atom
      integer, allocatable :: lste(:)                                   ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_short_atomic(:)             ! internal         
      integer max_num_neighbors_short_atomic                            ! internal
      integer, allocatable :: neighboridx_short_atomic(:,:)             ! internal
      integer, allocatable :: invneighboridx_short_atomic(:,:)          ! internal
      integer, allocatable :: atomindex_dummy(:)                        ! internal
!!
      real*8 nnshortforce(3,max_num_atoms)                              ! internal 
      real*8 nnshortforce_local(3,max_num_atoms,ndim)                   ! out 
      real*8 deshortdsfunc(max_num_atoms,maxnum_funcvalues_short_atomic)         ! dummy here 
      real*8 symfunctiondummy(maxnum_funcvalues_short_atomic,max_num_atoms)      ! internal
      real*8 symfunction_local(maxnum_funcvalues_short_atomic,max_num_atoms,ndim)! in
      real*8 lattice_local(3,3,ndim)                                    ! in
      real*8 xyzstruct_local(3,max_num_atoms,ndim)                      ! in
      real*8 shortforce(3,max_num_atoms)                                ! internal
      real*8 strs(3,3,maxnum_funcvalues_short_atomic,max_num_atoms)     ! internal dummy
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)! in
      real*8, allocatable :: lstb(:,:)                                  ! internal, xyz and r_ij 
      real*8, allocatable :: dsfuncdxyz(:,:,:,:)                        ! internal
!!
      logical lperiodic_local(ndim)                                     ! in
      logical lperiodic                                                 ! internal
      logical ldummy
!!
      shortforce(:,:)=0.0d0 ! just dummy here
      nnshortforce_local(:,:,:)=0.0d0 
!!
!!
      do i1=1,npoints
        num_atoms        =num_atoms_local(i1)
        zelem(:)         =zelem_local(i1,:)
        lperiodic        =lperiodic_local(i1)
!!
!! in principle dsfuncdxyz could be precalculated, but needs too much storage on disk
!!
        allocate(lsta(2,max_num_atoms))
        allocate(lstc(listdim))
        allocate(lste(listdim))
        allocate(lstb(listdim,4))
        allocate(num_neighbors_short_atomic(num_atoms))
!!
!! we must use the serial version here because we have already parallelized in a higher level routine
        call getneighborsatomic(num_atoms,&
          num_neighbors_short_atomic,zelem,max_num_neighbors_short_atomic,&
          lsta,lstc,lste,&
          maxcutoff_short_atomic,lattice_local(1,1,i1),xyzstruct_local(1,1,i1),&
          lstb,lperiodic)
!!
        allocate(dsfuncdxyz(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_short_atomic,3))
!! allocation of index arrays must be over num_atoms here because we cannot parallelize over atoms in geterror
        allocate(neighboridx_short_atomic(num_atoms,0:max_num_neighbors_short_atomic))  
        allocate(invneighboridx_short_atomic(num_atoms,max_num_atoms))  
        call getneighboridxatomic(num_atoms,listdim,&
          max_num_atoms,max_num_neighbors_short_atomic,&
          lsta,lstc,neighboridx_short_atomic,invneighboridx_short_atomic)
!!
!! get dsfuncdxyz:
!! Caution: symmetry functions cannot be used because here they are not scaled
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
           lambda_short_atomic,zeta_short_atomic,dsfuncdxyz,strs,lstb,&
           lperiodic,.true.,.false.,ldummy)
        deallocate(atomindex_dummy)
!!
        deallocate(lsta)
        deallocate(lstc)
        deallocate(lste)
        deallocate(lstb)
!!
!! scale dsfuncdxyz
        if(lscalesym)then
          call scaledsfunc(max_num_neighbors_short_atomic,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            nelem,num_atoms,minvalue_short_atomic,maxvalue_short_atomic,&
            scmin_short_atomic,scmax_short_atomic,&
            zelem,dsfuncdxyz,strs)
        endif
!!
!!
        call getshortforces(max_num_neighbors_short_atomic,num_atoms,&
          num_neighbors_short_atomic,neighboridx_short_atomic,zelem,&
          symfunction_local(1,1,i1),&
          dsfuncdxyz,deshortdsfunc,nnshortforce)
!!       
        nnshortforce_local(:,:,i1)=nnshortforce(:,:)

        deallocate(dsfuncdxyz)
        deallocate(neighboridx_short_atomic)  
        deallocate(invneighboridx_short_atomic)  
        deallocate(num_neighbors_short_atomic)
!!
      enddo ! i1
!!
      return
      end
