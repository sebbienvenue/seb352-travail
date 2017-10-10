!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate symfunction, dsfuncsdxyz and strs for natoms atoms
!! multipurpose subroutine

!! called by:
!! - calcfunctions.f90
!! - getallelectrostatic.f90
!! - getallforces.f90
!! - getallshortforces.f90
!! - optimize_atomic.f90
!! - optimize_short_combined.f90
!! - predictelec.f90
!!
      subroutine calconefunction_atomic(cutoff_type,max_num_neighbors_local,&
        max_num_atoms,n_start,natoms,atomindex,natomsdim,elementindex,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        nelem,zelem,listdim,&
        lsta,lstc,lste,invneighboridx_local,&
        function_type_local,symelement_local,&
        xyzstruct,symfunction,rmin,&
        funccutoff_local,eta_local,rshift_local,lambda_local,&
        zeta_local,dsfuncdxyz,strs,lstb,&
        lperiodic,ldoforces,ldostress,lrmin)
!!
      use fileunits
!! don't use globaloptions here
!! don't use symfunctions here
!!
      implicit none
!!
      integer n_start                                            ! in
      integer listdim                                            ! in
      integer cutoff_type                                        ! in
      integer maxnum_funcvalues_local                            ! in
      integer num_funcvalues_local(nelem)                        ! in
      integer max_num_atoms                                      ! in
      integer max_num_neighbors_local                            ! in
      integer zelem(max_num_atoms)                               ! in
      integer function_type_local(maxnum_funcvalues_local,nelem) ! in
      integer symelement_local(maxnum_funcvalues_local,2,nelem)  ! in
      integer nelem                                              ! in
      integer lsta(2,max_num_atoms)                              ! in, numbers of neighbors
      integer lstc(listdim)                                      ! in, identification of atom
      integer lste(listdim)                                      ! in, nuclear charge of atom
      integer iindex                                             ! internal
      integer i1,i2,i3,i4,i5                                     ! internal
      integer natoms                                             ! in
      integer atomindex(natoms)                                  ! in
      integer natomsdim                                          ! in
      integer jcount                                             ! internal
      integer elementindex(102)                                  ! in
      integer invneighboridx_local(natoms,max_num_atoms)         ! in
!!
      real*8 xyzstruct(3,max_num_atoms)                          ! in
      real*8 symfunction(maxnum_funcvalues_local,natomsdim)      ! out
      real*8 symfunction_temp(maxnum_funcvalues_local)           ! internal 
      real*8 dsfuncdxyz(maxnum_funcvalues_local,natomsdim,0:max_num_neighbors_local,3)  ! out
      real*8 dsfuncdxyz_temp(0:max_num_neighbors_local,3)        ! internal 
      real*8 strs(3,3,maxnum_funcvalues_local,natomsdim)         ! out
      real*8 strs_temp(3,3,maxnum_funcvalues_local)              ! internal 
      real*8 lstb(listdim,4)                                     ! in, xyz and r_ij
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)     ! in
      real*8 eta_local(maxnum_funcvalues_local,nelem)            ! in
      real*8 rshift_local(maxnum_funcvalues_local,nelem)         ! in
      real*8 lambda_local(maxnum_funcvalues_local,nelem)         ! in
      real*8 zeta_local(maxnum_funcvalues_local,nelem)           ! in
      real*8 rmin                                                ! in

      logical lperiodic                                          ! in
      logical ldoforces                                          ! in
      logical ldostress                                          ! in
      logical lrmin                                              ! in/out
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!    initialization
      symfunction(:,:)   = 0.0d0
      dsfuncdxyz(:,:,:,:)= 0.0d0
      strs(:,:,:,:)      = 0.0d0
      lrmin              = .true. ! default for sufficiently long bonds         
!!
!! check for obviously wrong structures with too short bonds
!! do this only for the atoms of this process here
      do i1=1,natoms
        do i2=lsta(1,i1),lsta(2,i1)
          if(lstb(i2,4).le.rmin)then
            lrmin=.false.
!!            stop
          endif
        enddo
      enddo
!!
!! check for isolated atoms without neighbors in the structure

      if(lperiodic.or.(natoms.gt.1))then
        do i1=1,natoms
          if(lsta(1,i1).eq.0)then
            write(ounit,*)'ERROR: atom ',n_start+i1-1,' has no neighbors'
            stop !'
          endif
        enddo
      endif
!!
!! calculate the symmetry function values
      jcount=n_start
      do i1=1,natoms ! loop over all atoms of this process
        iindex=elementindex(zelem(jcount))
        do i2=1,num_funcvalues_local(iindex) ! loop over all symmetry functions
!!
          dsfuncdxyz_temp(:,:)=0.0d0
          symfunction_temp(:)=0.0d0
!!
          call getatomsymfunctions(i1,i2,iindex,natoms,atomindex,natomsdim,&
            max_num_atoms,max_num_neighbors_local,& 
            invneighboridx_local,jcount,listdim,lsta,lstc,lste,&
            symelement_local,maxnum_funcvalues_local,&
            cutoff_type,nelem,function_type_local,& 
            lstb,funccutoff_local,xyzstruct,symfunction_temp,dsfuncdxyz_temp,strs_temp,& 
            eta_local,zeta_local,lambda_local,rshift_local,rmin,& 
            ldoforces,ldostress)
!!
          symfunction(i2,i1)   =symfunction_temp(i2)
          do i3=0,max_num_neighbors_local,1
            do i4=1,3
              dsfuncdxyz(i2,i1,i3,i4)=dsfuncdxyz(i2,i1,i3,i4)+dsfuncdxyz_temp(i3,i4)
              do i5=1,3
                strs(i5,i4,i2,i1)    =strs(i5,i4,i2,i1)+strs_temp(i5,i4,i2)
              enddo ! i5
            enddo ! i4
          enddo ! i3
!!
        enddo ! i2
        jcount=jcount+1
!!
      enddo ! i1
!!
      return
      end
