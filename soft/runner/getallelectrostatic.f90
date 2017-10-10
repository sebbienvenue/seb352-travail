!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - ewaldenergies_para.f90: 
!!
      subroutine getallelectrostatic(ndim1,npoints,&
        num_atoms_local,&
        zelem_local,lattice_local,minvalue_elec,maxvalue_elec,&
        nnatomcharge_local,xyzstruct_local,nnewald_local,&
        symfunctione_local,nnelecforce_local,&
        lperiodic_local)
!!
      use fileunits
      use nnflags
      use globaloptions 
      use symfunctions
      use nnewald
!!
      implicit none
!!
      integer ndim1
      integer npoints
      integer num_atoms_local(ndim1)
      integer num_atoms
      integer zelem_local(ndim1,max_num_atoms)        ! in
      integer zelem(max_num_atoms)                    ! internal
      integer i1,i2                                   ! internal
      integer, allocatable :: lsta(:,:)               ! numbers of neighbors
      integer, allocatable :: lstc(:)                 ! identification of atom
      integer, allocatable :: lste(:)                 ! nuclear charge of atom
      integer, allocatable :: num_neighbors_elec(:)   ! internal                  
      integer max_num_neighbors_elec                  ! internal
      integer, allocatable :: neighboridx_elec(:,:)   ! internal
      integer, allocatable :: invneighboridx_elec(:,:)! internal
      integer, allocatable :: atomindex_dummy(:)      ! internal
!!
      real*8 nnewald_local(ndim1)                     ! out
      real*8 elecenergy 
      real*8 xyzstruct_local(3,max_num_atoms,ndim1)
      real*8 nnatomcharge(max_num_atoms)
      real*8 nnatomcharge_local(ndim1,max_num_atoms)
      real*8 lattice_local(3,3,ndim1)                 ! in
      real*8 nnelecforce(3,max_num_atoms)             ! internal 
      real*8 nnelecforce_local(3,max_num_atoms,ndim1) ! out 
      real*8 symfunctione_local(maxnum_funcvalues_elec,max_num_atoms,ndim1) ! in
      real*8 symfunctionedummy(maxnum_funcvalues_elec,max_num_atoms)        ! internal
      real*8 strse(3,3,maxnum_funcvalues_elec,max_num_atoms)                   ! internal
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)         ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)         ! in
      real*8, allocatable :: lstb(:,:)                                  ! xyz and r_ij 
      real*8, allocatable :: dsfuncdxyze(:,:,:,:)    
      real*8, allocatable :: dchargedxyz(:,:,:)    
!!
      logical lperiodic_local(ndim1)                  ! in
      logical lperiodic                               ! internal
      logical ldummy
!!
!!
!! initializations
      elecenergy            =0.0d0
      nnelecforce(:,:)       =0.0d0  
      symfunctionedummy(:,:) =0.0d0
!!
      do i1=1,npoints
        num_atoms        =num_atoms_local(i1)
        nnatomcharge(:)  =nnatomcharge_local(i1,:)
        lperiodic        =lperiodic_local(i1)
        zelem(:)         =zelem_local(i1,:) 
        nnelecforce(:,:) =0.0d0  
!!
        allocate(lsta(2,max_num_atoms))
        allocate(lstc(listdim))
        allocate(lste(listdim))
        allocate(lstb(listdim,4))
        allocate(num_neighbors_elec(num_atoms))
!!
        call getneighborsatomic(&
          num_atoms,num_neighbors_elec,zelem,max_num_neighbors_elec,&
          lsta,lstc,lste,&
          maxcutoff_elec,lattice_local(1,1,i1),xyzstruct_local(1,1,i1),&
          lstb,lperiodic)
!!
        allocate(neighboridx_elec(num_atoms,0:max_num_neighbors_elec))  
        allocate(invneighboridx_elec(num_atoms,max_num_atoms))  
        allocate(dsfuncdxyze(maxnum_funcvalues_elec,max_num_atoms,0:max_num_neighbors_elec,3))
        call getneighboridxatomic(num_atoms,listdim,&
          max_num_atoms,max_num_neighbors_elec,&
          lsta,lstc,neighboridx_elec,invneighboridx_elec)
!!
        if(luseforces)then
!! get dsfuncdxyze:
!! in principle dsfuncdxyze could be precalculated, but needs too much storage on disk
!! Caution: symmetry functions cannot be used because here they are not scaled
          allocate(atomindex_dummy(num_atoms))
          do i2=1,num_atoms
            atomindex_dummy(i2)=i2
          enddo
          call calconefunction_atomic(cutoff_type,max_num_neighbors_elec,&
            max_num_atoms,1,num_atoms,atomindex_dummy,num_atoms,max_num_atoms,elementindex,&
            maxnum_funcvalues_elec,num_funcvalues_elec, &
            nelem,zelem,listdim,&
            lsta,lstc,lste,invneighboridx_elec,&
            function_type_elec,symelement_elec,&
            xyzstruct_local(1,1,i1),symfunctionedummy,0.d0,&
            funccutoff_elec,eta_elec,rshift_elec,lambda_elec,zeta_elec,dsfuncdxyze,strse,lstb,&
            lperiodic,.true.,.false.,ldummy)
          deallocate(atomindex_dummy)
!!
!! scale derivatives dsfuncdxyze and strse
          if(lscalesym)then
            call scaledsfunc(max_num_neighbors_elec,&
              maxnum_funcvalues_elec,num_funcvalues_elec,&
              nelem,num_atoms,minvalue_elec,maxvalue_elec,&
              scmin_elec,scmax_elec,&
              zelem,dsfuncdxyze,strse)
          endif ! lscalesym
        endif ! luseforces
!!
        if(lperiodic)then
          allocate(dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)) 
          dchargedxyz(:,:,:)     =0.0d0  
          if(luseforces.and.(nn_type_elec.ne.3))then
!! calculate dchargedxyz for forces
          call getdchargedxyz(max_num_neighbors_elec,&
               num_neighbors_elec,neighboridx_elec,num_atoms,zelem,dsfuncdxyze,&
               dchargedxyz,symfunctione_local(1,1,i1))
          endif
          call getewaldenergy(max_num_neighbors_elec,&
               neighboridx_elec,num_neighbors_elec,num_atoms,zelem,&
               lattice_local(1,1,i1),xyzstruct_local(1,1,i1),&
               nnatomcharge,elecenergy,&
               dchargedxyz,nnelecforce,.false.)
!!
          deallocate(dchargedxyz) 
        else ! not periodic
          call electrostatic(num_atoms,&
               nnatomcharge,xyzstruct_local(1,1,i1),elecenergy)
!!
!! calculate electrostatic forces
          if(luseforces)then 
            nnelecforce(:,:) =0.0d0 ! important!!! 
            call getcoulombforces(max_num_neighbors_elec,&
              num_neighbors_elec,neighboridx_elec,num_atoms,zelem,&
              dsfuncdxyze,symfunctione_local(1,1,i1),&
              nnatomcharge,xyzstruct_local(1,1,i1),nnelecforce)
          endif ! luseforces
!!
        endif ! lperiodic
!!
!! normalize the Ewald energy per atom to be consistent with the DFT data
        nnewald_local(i1)=elecenergy/dble(num_atoms_local(i1))
!!        nnelecforce_local(:,:,i1)=nnelecforce(:,:)/dble(num_atoms_local(i1))
        nnelecforce_local(:,:,i1)=nnelecforce(:,:)
!!
        deallocate(lsta)
        deallocate(lstc)
        deallocate(lste)
        deallocate(lstb)
        deallocate(invneighboridx_elec)  
        deallocate(num_neighbors_elec)
        deallocate(neighboridx_elec)  
        deallocate(dsfuncdxyze)
!!
      enddo ! i1=1,npoints
!!
      return
      end
