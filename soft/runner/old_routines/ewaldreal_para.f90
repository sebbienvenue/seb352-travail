!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getewaldenergy.f90
!!
!! note the dependence between the neighborlist cutoff and ewaldalpha!!!
!!
      subroutine ewaldreal_para(max_num_neighbors_elec,num_neighbors_elec,&
           neighboridx_elec,natoms,atomindex,&
           n_start,n_end,num_atoms,zelem,&
           atomcharge,ereal,dchargedxyz,&
           lattice,xyzstruct,erealforce,&
           lperiodic,ldoforces)
!!
      use globaloptions
!!
      implicit none
!!
      integer i,j,n
      integer i1,i2
      integer num_atoms
      integer max_num_neighbors_elec
      integer num_neighbors_elec(num_atoms)
      integer neighboridx_elec(num_atoms,0:max_num_neighbors_elec)
!!      integer lsta(2,max_num_atoms)                          ! pointer
      integer lsta(2,natoms)                          ! pointer
      integer lstc(listdim)                                  ! neighbor label
      integer lste(listdim)                                  ! nuclear charge of neighbor
      integer zelem(max_num_atoms)
      integer natoms                                         ! in
      integer atomindex(natoms)                              ! in
      integer n_start                                        ! in
      integer n_end                                          ! in
!!
      real*8 atomcharge(max_num_atoms)
      real*8 ereal
      real*8 lattice(3,3)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 rij
      real*8 comperrfct                                      ! erfc(x)
      real*8 dcomperrfct                                     ! d erfc(x)/d x
      real*8 lstb(listdim,4)                                 ! x,y,z, rij of neighbor
      real*8 erealforce(3,max_num_atoms)                     ! out
      real*8 dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)      ! in 
      real*8 deltaxj
      real*8 deltayj
      real*8 deltazj
      real*8 drijdxi
      real*8 drijdyi
      real*8 drijdzi
      real*8 drijdxj
      real*8 drijdyj
      real*8 drijdzj
      real*8 temp1
      real*8 drijd(3)
!!
      logical lperiodic                                  ! in
      logical ldoforces                                  ! in
!!
!! initializations
      ereal=0.0d0
!!
!! calculate neighbor list for natoms
      call neighbor_para(n_start,n_end,&
              num_atoms,zelem,&
              lsta,lstb,lstc,lste,&
              ewaldcutoff,lattice,xyzstruct,lperiodic)
!! using the neighbor list automatically prevents self-interaction with R_ij=0 (division!)
!!
!! loop over all atoms
      do i=1,natoms
!! loop over all neighbors
        do j=lsta(1,i),lsta(2,i)
          n=lstc(j)
          rij=lstb(j,4)
!!
          temp1=ewaldalpha*rij  ! abbreviation
!!
          ereal=ereal+atomcharge(atomindex(i))*atomcharge(n)*&
                comperrfct(temp1)/rij
!!
!! TODO: Put this in extra double loop to avoid hundreds of ifs
          if(ldoforces)then
            deltaxj=-1.d0*(xyzstruct(1,atomindex(i))-lstb(j,1))
            deltayj=-1.d0*(xyzstruct(2,atomindex(i))-lstb(j,2))
            deltazj=-1.d0*(xyzstruct(3,atomindex(i))-lstb(j,3))
            drijdxi=-deltaxj/rij
            drijdyi=-deltayj/rij
            drijdzi=-deltazj/rij
            drijdxj=-1.d0*drijdxi
            drijdyj=-1.d0*drijdyi
            drijdzj=-1.d0*drijdzi
!!
            do i1=1,num_neighbors_elec(i)
              if(i1.eq.atomindex(i)) then
                drijd(1)=drijdxi
                drijd(2)=drijdyi
                drijd(3)=drijdzi
              elseif(i1.eq.n) then
                drijd(1)=drijdxj
                drijd(2)=drijdyj
                drijd(3)=drijdzj
              else
                drijd(:)=0.0d0
              endif
!!
              do i2=1,3
            erealforce(i2,neighboridx_elec(i,i1))=erealforce(i2,neighboridx_elec(i,i1))      &
              -((atomcharge(n)*dchargedxyz(atomindex(i),i1,i2)  &
!! CHECK: dchargedxyz(n) geht so nicht
                +atomcharge(atomindex(i))*dchargedxyz(n,i1,i2)) &
              *comperrfct(temp1)/rij                 &
              +atomcharge(atomindex(i))*atomcharge(n)*          &
              (dcomperrfct(temp1)*ewaldalpha*drijd(i2)/rij   &
              -comperrfct(temp1)*drijd(i2)/rij**2 ))
              enddo ! i2
            enddo ! i1
          endif ! ldoforces
        enddo
      enddo
!!
!! remove double counting
      ereal=0.5d0*ereal
!! remove double counting also for erealforce
      erealforce(:,:)=0.5d0*erealforce(:,:)
!!
      return
      end
