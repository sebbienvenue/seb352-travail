!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - prediction (3x)
!! - optimize_short_combined.f90
!! - getallelectrostatic.f90
!! - getallshortforces.f90
!!
      subroutine calconefunction_paraAndi(cutoff_type,max_num_neighbors,&
           max_num_atoms,n_start,n_end,natoms,natomsdim,elementindex,&
           maxnum_funcvalues_local,num_funcvalues_local,&
           nelem,num_atoms,zelem,listdim,&
           lsta,lstc,lste,invneighboridx,&
           function_type_local,symelement_local,lattice,&
           xyzstruct,symfunction,maxcutoff_local,rmin,&
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
      integer listdim                                        ! in
      integer cutoff_type                                    ! in
      integer maxnum_funcvalues_local                        ! in
      integer num_funcvalues_local(nelem)                    ! in
      integer max_num_atoms                                  ! in
      integer max_num_neighbors                              ! in
      integer num_atoms                                      ! in
      integer zelem(max_num_atoms)                           ! in
      integer function_type_local(maxnum_funcvalues_local,nelem) ! in
      integer symelement_local(maxnum_funcvalues_local,2,nelem)  ! in
      integer nelem                                              ! in
      integer lsta(2,max_num_atoms)                              ! numbers of neighbors
      integer lstc(listdim)                                      ! identification of atom
      integer lste(listdim)                                      ! nuclear charge of atom
      integer iindex                                         ! internal
      integer i1,i2,i3,i4                                    ! internal
      integer n_start                                        ! in
      integer n_end                                          ! in
      integer natoms                                         ! in
      integer natomsdim                                      ! in
      integer jcount                                         ! internal
      integer elementindex(102)                              ! in
      integer invneighboridx(natoms,max_num_atoms)           ! in
!!
      real*8 lattice(3,3)                                                   ! in
      real*8 xyzstruct(3,max_num_atoms)                                     ! in
      real*8 symfunction(maxnum_funcvalues_local,natomsdim)                 ! out
      real*8 dsfuncdxyz(maxnum_funcvalues_local,natomsdim,0:max_num_neighbors,3)  ! out
      real*8 maxcutoff_local                                                ! in
      real*8 lstb(listdim,4)                                                ! xyz and r_ij
      real*8 pi                                                             ! internal
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)                ! in
      real*8 eta_local(maxnum_funcvalues_local,nelem)                       ! in
      real*8 rshift_local(maxnum_funcvalues_local,nelem)                    ! in
      real*8 lambda_local(maxnum_funcvalues_local,nelem)                    ! in
      real*8 zeta_local(maxnum_funcvalues_local,nelem)                      ! in
      real*8 strs(3,3,maxnum_funcvalues_local,natomsdim)                  ! out
      real*8 rmin

      logical lperiodic                      ! in
      logical ldoforces                      ! in
      logical ldostress                      ! in
      logical lrmin                          ! in/out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!    initialization
      pi                 =3.141592654d0
      symfunction(:,:)   =0.0d0
      dsfuncdxyz(:,:,:,:)=0.0d0
      strs(:,:,:,:)      =0.0d0
      lrmin              =.true.
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(function_type_local(i2,iindex).eq.1)then ! radial function
!!
            call atomsymfunction1(i1,i2,iindex,natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,&
              jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
              lstb,funccutoff_local,pi,xyzstruct,symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.2)then ! radial function
!!
            call atomsymfunction2(i1,i2,iindex,natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,&
              jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
              cutoff_type,lstb,funccutoff_local,pi,xyzstruct,eta_local,rshift_local,&
              symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.3)then
!!
            call atomsymfunction3Andi(i1,i2,iindex,&
              natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,&
              jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
              cutoff_type,lstb,funccutoff_local,pi,xyzstruct,eta_local,zeta_local,&
              lambda_local,rmin,symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress,lrmin)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.4)then ! radial function
!!
            call atomsymfunction4(i1,i2,iindex,natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,&
              jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
              cutoff_type,lstb,funccutoff_local,pi,xyzstruct,eta_local,&
              symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.5)then ! just Cartesian coordinate
!!
            call atomsymfunction5(i1,i2,iindex,natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,jcount,maxnum_funcvalues_local,&
              xyzstruct,eta_local,symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.6)then ! radial function
!!
            call atomsymfunction6(i1,i2,iindex,natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,&
              jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
              lstb,funccutoff_local,xyzstruct,&
              symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.7)then ! angular function
!!
            write(ounit,*)'Error: function type not implemented in calconefunction ',function_type_local(i1,iindex)
            stop !'
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.8)then ! angular function
!!
            call atomsymfunction8(i1,i2,iindex,natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,&
              jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
              cutoff_type,lstb,funccutoff_local,pi,xyzstruct,eta_local,rshift_local,&
              lambda_local,rmin,symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress,lrmin)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          elseif(function_type_local(i2,iindex).eq.9)then ! angular function
!!
            call atomsymfunction9(i1,i2,iindex,natoms,natomsdim,nelem,max_num_atoms,max_num_neighbors,&
              invneighboridx,&
              jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
              cutoff_type,lstb,funccutoff_local,pi,xyzstruct,eta_local,zeta_local,&
              lambda_local,rmin,symfunction,dsfuncdxyz,strs,&
              ldoforces,ldostress,lrmin)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else
            write(*,*)'Error: function_type_local not implemented ',function_type_local(i2,iindex)
          endif
        enddo ! i2
        jcount=jcount+1
!!
      enddo ! i1
!!
      return
      end
