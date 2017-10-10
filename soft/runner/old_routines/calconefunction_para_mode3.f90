!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate symfunctions, dsfuncdxyz and strs for natomsdim atoms

!! multipurpose subroutine

!! called by:
!! - prediction (3x)
!! - predictionpair.f90 (2x)
!! - optimize_short_combined.f90
!! - getdshortdw.f90
!! - calcfunctions.f90 (2x)
!! - calcpairfunctions.f90
!! - getallelectrostatic.f90
!! - getallshortforces.f90
!!
      subroutine calconefunction_para_mode3(cutoff_type,max_num_neighbors,&
           max_num_atoms,n_start,natoms,natomsdim,elementindex,&
           maxnum_funcvalues_local,num_funcvalues_local,&
           nelem,zelem,listdim,&
           lsta,lstc,lste,invneighboridx,&
           function_type_local,symelement_local,lattice,&
           xyzstruct,symfunction,rmin,&
           funccutoff_local,eta_local,rshift_local,lambda_local,zeta_local,strs,lstb,&
           lperiodic,ldoforces,ldostress,lrmin)
!!
      use fileunits
!! don't use globaloptions here
!! don't use symfunctions here
!!
      implicit none
!!
      integer listdim                                            ! in
      integer cutoff_type                                        ! in
      integer maxnum_funcvalues_local                            ! in
      integer num_funcvalues_local(nelem)                        ! in
      integer max_num_atoms                                      ! in
      integer max_num_neighbors                                  ! in
      integer zelem(max_num_atoms)                               ! in
      integer function_type_local(maxnum_funcvalues_local,nelem) ! in
      integer symelement_local(maxnum_funcvalues_local,2,nelem)  ! in
      integer nelem                                              ! in
      integer lsta(2,max_num_atoms)                              ! in, numbers of neighbors
      integer lstc(listdim)                                      ! in, identification of atom
      integer lste(listdim)                                      ! in, nuclear charge of atom
      integer iindex                                             ! internal
      integer i1,i2                                              ! internal
      integer n_start                                            ! in
      integer natoms                                             ! in
      integer natomsdim                                          ! in
      integer jcount                                             ! internal
      integer elementindex(102)                                  ! in
      integer invneighboridx(natoms,max_num_atoms)               ! in
!!
      real*8 lattice(3,3)                                        ! in
      real*8 xyzstruct(3,max_num_atoms)                          ! in
      real*8 symfunction(maxnum_funcvalues_local,natomsdim)      ! out
      real*8 dsfuncdxyz_temp(0:max_num_neighbors,3)              ! out
      real*8 lstb(listdim,4)                                     ! in, xyz and r_ij 
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)     ! in
      real*8 eta_local(maxnum_funcvalues_local,nelem)            ! in
      real*8 rshift_local(maxnum_funcvalues_local,nelem)         ! in
      real*8 lambda_local(maxnum_funcvalues_local,nelem)         ! in
      real*8 zeta_local(maxnum_funcvalues_local,nelem)           ! in
      real*8 strs(3,3,maxnum_funcvalues_local,natomsdim)         ! out
      real*8 rmin

      logical lperiodic                                          ! in
      logical ldoforces                                          ! in
      logical ldostress                                          ! in
      logical lrmin                                              ! in/out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!    initialization
!!      write(ounit,*)'calconefunction_para natoms,natomsdim ',natoms,natomsdim
      symfunction(:,:)   =0.0d0
      strs(:,:,:,:)      =0.0d0
      lrmin              =.true.
!!
!! check for obviously wrong structures with too short bonds
!! do this only for the atoms of this process here
      do i1=1,natoms
        do i2=lsta(1,i1),lsta(2,i1)
          if(lstb(i2,4).le.rmin)then
            lrmin=.false.
!!            write(ounit,*)'Error: rij lt rmin ',jcount,i2,lstb(i2,4)
!!            write(ounit,*)xyzstruct(1,n_start+i1-1),lstb(i2,1)
!!            write(ounit,*)xyzstruct(2,n_start+i1-1),lstb(i2,2)
!!            write(ounit,*)xyzstruct(3,n_start+i1-1),lstb(i2,3)
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
!! calculate the symmetry function values for natoms atoms
      jcount=n_start
      do i1=1,natoms ! loop over all atoms of this process
        iindex=elementindex(zelem(jcount))
        do i2=1,num_funcvalues_local(iindex) ! loop over all symmetry functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dsfuncdxyz is only dummy here as the needed derivatives are now calculated on the fly in getshortforces_para
          dsfuncdxyz_temp(:,:)=0.0d0
!!
          call getatomsymfunctions(i1,i2,iindex,natoms,natomsdim,max_num_atoms,max_num_neighbors,&
            invneighboridx,jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
            cutoff_type,nelem,function_type_local,&
            lstb,funccutoff_local,xyzstruct,symfunction,dsfuncdxyz_temp,strs,&
            eta_local,zeta_local,lambda_local,rshift_local,&
            ldoforces,ldostress,lrmin)
!!
        enddo ! i2
        jcount=jcount+1
!!
      enddo ! i1
!!
      return
      end
