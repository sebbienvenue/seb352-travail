!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! Purpose: Calculate electrostatic energy and forces with standard Ewald summation
!!
!! called by:
!! - predictelec.f90
!!
      subroutine getewald(n_start,n_end,natoms,atomindex,&
        num_funcvalues_elec,zelem,num_atoms,&
        max_num_neighbors_elec,invneighboridx_elec,&
        nnatomcharge,lattice,xyzstruct,dchargedsfunc,dsfuncdxyze,&
        maxcutoff_local,nnelecenergy,erecip,nnelecforce,ldoforces_local,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5                                            ! internal
      integer natoms                                                    ! in
      integer atomindex(natoms)                                         ! in
      integer max_num_neighbors_elec                                    ! in
      integer invneighboridx_elec(natoms,max_num_atoms)                 ! in 
      integer zelem(max_num_atoms)                                      ! in
      integer num_funcvalues_elec(nelem)                                ! in
      integer num_atoms                                                 ! in
      integer n_start                                                   ! in
      integer n_end                                                     ! in
!!
      real*8 nnatomcharge(max_num_atoms)                                ! in
      real*8 xyzstruct(3,max_num_atoms)                                 ! in
      real*8 nnelecenergy                                               ! in/out
      real*8 nnelecforce(3,max_num_atoms)                               ! in/out
      real*8 dsfuncdxyze(maxnum_funcvalues_elec,natoms,0:max_num_neighbors_elec,3) ! in
      real*8 dchargedsfunc(natoms,maxnum_funcvalues_elec)               ! in
      real*8 maxcutoff_local                                            ! in
      real*8 eself                                                      ! internal
      real*8 ereal                                                      ! internal
      real*8 erecip                                                     ! internal
      real*8 ewaldcorr                                                  ! internal
      real*8 eselfforce(3,max_num_atoms)                                ! internal
      real*8 erealforce(3,max_num_atoms)                                ! internal
      real*8 erecipforce(3,max_num_atoms)                               ! internal
      real*8 ewaldcorrforce(3,max_num_atoms)                            ! internal
      real*8 lattice(3,3)                                               ! in
!!
      logical ldoforces_local                                           ! in
      logical lperiodic                                                 ! in
!!
!!======================================================================
!! calculate the real space part for natoms atoms 
!!======================================================================
      call ewaldreal_mode3(n_start,n_end,max_num_neighbors_elec,&
        num_atoms,zelem,invneighboridx_elec,&
        natoms,atomindex,num_funcvalues_elec,&
        nnatomcharge,dchargedsfunc,dsfuncdxyze,&
        lattice,xyzstruct,erealforce,ereal,ldoforces_local,lperiodic)
!! for debugging:
!      ereal=0.0d0
!      erealforce(:,:)=0.0d0

!!======================================================================
!! calculate the reciprocal space part for natoms atoms
!!======================================================================
      call ewaldrecip_mode3(max_num_neighbors_elec,num_atoms,&
        zelem,invneighboridx_elec,&
        natoms,atomindex,num_funcvalues_elec,&
        nnatomcharge,dchargedsfunc,dsfuncdxyze,&
        lattice,xyzstruct,erecipforce,erecip,ldoforces_local)
!! for debugging:
!      erecip=0.0d0
!      erecipforce(:,:)=0.0d0

!!======================================================================
!! calculate the self energy correction part for natoms atoms
!!======================================================================
      call ewaldself_mode3(max_num_neighbors_elec,&
        num_atoms,zelem,invneighboridx_elec,&
        natoms,atomindex,num_funcvalues_elec,&
        nnatomcharge,dchargedsfunc,dsfuncdxyze,&
        eselfforce,eself,ldoforces_local)
!! for debugging
!      eself=0.0d0
!      eselfforce(:,:)=0.0d0
!!
!!======================================================================
!! sum the ewald energy and force contributions
!! elecenergy is the total energy here, not per atom!
!!======================================================================
!! Caution: erecip cannot be added here because of different double counting
       nnelecenergy = ereal+eself  ! note sign of eself!
       nnelecforce(:,:)=erealforce(:,:)+erecipforce(:,:)+eselfforce(:,:)
!!
!      write(ounit,*)'getewald: ereal energy ',ereal
!      write(ounit,*)'getewald: erecip energy ',erecip
!      write(ounit,*)'getewald: eself energy ',eself
!      write(ounit,*)'getewald: erealforce ',erealforce
!      write(ounit,*)'getewald: erecipforce ',erecipforce
!      write(ounit,*)'getewald: eselfforce ',eselfforce
      return
      end
