!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: Calculate reciprocal-space part of standard Ewald summation, parallel case
!! FIXME: get this properly parallelized or replace by particle mesh Ewald
!!
!! called by:
!! - getewald.f90
!!
      subroutine ewaldrecip_mode3(max_num_neighbors_elec,num_atoms,&
        zelem,invneighboridx_elec,&
        natoms,atomindex,num_funcvalues_elec,&
        nnatomcharge,dchargedsfunc,dsfuncdxyze,&
        lattice,xyzstruct,erecipforce,erecip,ldoforces_local)
!!
      use globaloptions
      use fileunits
      use nnflags
!!
      implicit none
!!
      integer max_num_neighbors_elec                      ! in
      integer invneighboridx_elec(natoms,max_num_atoms)   ! in
      integer natoms                                      ! in
      integer atomindex(natoms)                           ! in
      integer num_funcvalues_elec(nelem)                      ! in
      integer zelem(max_num_atoms)                        ! in
      integer num_atoms                                   ! in
!!
      integer i1,i2,i3,i4,i5,i6,i7                        ! internal
!!
      real*8 nnatomcharge(max_num_atoms)                  ! in
      real*8 erecipforce(3,max_num_atoms)                 ! out
      real*8 erecip                                       ! out
      real*8 sqrtpiinv
      parameter(sqrtpiinv=0.564189583d0)
      real*8 dsfuncdxyze(maxnum_funcvalues_elec,natoms,0:max_num_neighbors_elec,3) ! in
      real*8 dchargedsfunc(natoms,maxnum_funcvalues_elec)     ! in
      real*8 dchargedxyz                                  ! internal
      real*8 xyzstruct(3,max_num_atoms)                   ! in
      real*8 kvec(3)                                      ! internal
      real*8 ksquared                                     ! internal
      real*8 twopi                                        ! internal
      real*8 lattice(3,3)                                 ! in
      real*8 reclattice(3,3)                              ! internal
      real*8 determinant                                  ! internal
      real*8 sumsindev                                    ! internal
      real*8 sumcosdev                                    ! internal
      real*8 krj                                          ! internal
      real*8 dkrj                                         ! internal
      real*8 c1,c2,c3                                     ! internal abbreviations
!!
      logical ldoforces_local                             ! in
!!
!!======================================================================
!! initializations
!!======================================================================
      erecip          =0.0d0
      erecipforce(:,:)=0.0d0
      twopi           =6.283185307d0
!!
!!======================================================================
!! calculate the reciprocal lattice vectors
!!======================================================================
      call getreclat(lattice,reclattice,determinant)
!!
      do i1=-ewaldkmax,ewaldkmax
        do i2=-ewaldkmax,ewaldkmax
          do i3=-ewaldkmax,ewaldkmax
!!======================================================================
!! generate the reciprocal lattice vectors
!! later we should implement here a more clever way that automatically adapts to the
!! lengths of the reciprocal lattice vectors to really have a sphere in reciprocal space
!!======================================================================
            kvec(1)=dble(i1)*reclattice(1,1)&
                   +dble(i2)*reclattice(2,1)&
                   +dble(i3)*reclattice(3,1)
            kvec(2)=dble(i1)*reclattice(1,2)&
                   +dble(i2)*reclattice(2,2)&
                   +dble(i3)*reclattice(3,2)
            kvec(3)=dble(i1)*reclattice(1,3)&
                   +dble(i2)*reclattice(2,3)&
                   +dble(i3)*reclattice(3,3)
            ksquared=kvec(1)**2 + kvec(2)**2 + kvec(3)**2 ! in Bohr-2
            c3=(twopi/(ksquared*determinant))*dexp(-1.d0*ksquared/(4.d0*(ewaldalpha)**2))
!!
            if((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0))then
            else
              c1=0.0d0
              c2=0.0d0
!!======================================================================
!! calculate the reciprocal space energy 
!!======================================================================
              do i4=1,num_atoms ! loop over all atoms, breaking down to natoms is not possible because of sumcos**2 and sumsin**2 
                krj=kvec(1)*xyzstruct(1,i4)&
                   +kvec(2)*xyzstruct(2,i4)&
                   +kvec(3)*xyzstruct(3,i4)
                c1=c1+nnatomcharge(i4)*dcos(krj)
                c2=c2+nnatomcharge(i4)*dsin(krj)
              enddo !i
              erecip=erecip+c3*(c1**2+c2**2)
!!
!!======================================================================
!! calculate forces
!!======================================================================
              if(ldoforces_local)then
!! These sums are different for all derivatives (derivatives are defined by i5 and i4)
                do i5=1,num_atoms
                  do i4=1,3
                    sumcosdev=0.0d0
                    sumsindev=0.0d0
!! now only loop over the block of atoms of this process
                    do i6=1,natoms
                      krj=kvec(1)*xyzstruct(1,atomindex(i6))&
                         +kvec(2)*xyzstruct(2,atomindex(i6))&
                         +kvec(3)*xyzstruct(3,atomindex(i6))
!!
!! calculate the derivative d kRj /d alpha
                      if(i5.eq.atomindex(i6))then
                        dkrj=kvec(i4)
                      else
                        dkrj=0.0d0
                      endif
!!
!! calculate on the fly the derivative of the charge with respect to direction i5,i4
                      dchargedxyz=0.0d0
                      if((nn_type_elec.ne.3).and.(nn_type_elec.ne.4))then
                        do i7=1,num_funcvalues_elec(elementindex(zelem(atomindex(i6))))
                          dchargedxyz=dchargedxyz+dchargedsfunc(i6,i7)*dsfuncdxyze(i7,i6,invneighboridx_elec(i6,i5),i4)
                        enddo ! i7
                      endif
                      sumcosdev=sumcosdev+dcos(krj)*dchargedxyz-nnatomcharge(atomindex(i6))*dsin(krj)*dkrj
                      sumsindev=sumsindev+dsin(krj)*dchargedxyz+nnatomcharge(atomindex(i6))*dcos(krj)*dkrj
                    enddo ! i6
!!
                    erecipforce(i4,i5)=erecipforce(i4,i5)&
                      -c3*(2.0d0*c1*sumcosdev +2.0d0*c2*sumsindev)
!!
                  enddo ! i4
                enddo ! i5
!!
              endif ! ldoforces_local
            endif !((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0))
!!
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
      return
      end
