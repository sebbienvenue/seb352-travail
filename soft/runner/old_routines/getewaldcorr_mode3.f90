!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getewald.f90
!!
      subroutine getewaldcorr_mode3(natoms,atomindex,n_start,n_end,invneighboridx,&
        max_num_neighbors,num_atoms,zelem,num_funcvaluese,&
        nnatomcharge,dchargedsfunc,dsfuncdxyze,lattice,xyzstruct,&
        lperiodic,ldoforces,ewaldcorr,ewaldcorrforce)
!!
      use globaloptions
!!
      implicit none
!!
      integer i,j,n
      integer i1,i2,i3
      integer max_num_neighbors
      integer num_atoms
      integer lsta(2,max_num_atoms)                          ! pointer
      integer lstc(listdim)                                  ! neighbor label
      integer lste(listdim)                                  ! nuclear charge of neighbor
      integer zelem(max_num_atoms)
      integer natoms                                         ! in
      integer n_start                                        ! in
      integer n_end                                          ! in
      integer atomindex(natoms)                              ! in
      integer num_funcvaluese(nelem)                         ! in
      integer invneighboridx(natoms,max_num_atoms)           ! in
!!
      real*8 nnatomcharge(max_num_atoms)
      real*8 ewaldcorr
      real*8 lattice(3,3)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 rij
      real*8 lstb(listdim,4)                                 ! x,y,z, rij of neighbor
      real*8 ewaldcorrforce(3,max_num_atoms)                 ! out
      real*8 deltaxj
      real*8 deltayj
      real*8 deltazj
      real*8 drijdxi
      real*8 drijdyi
      real*8 drijdzi
      real*8 drijdxj
      real*8 drijdyj
      real*8 drijdzj
      real*8 drijd(3)
      real*8 fscreen                                     ! internal
      real*8 fscreenderiv                                ! internal
      real*8 dsfuncdxyze(maxnum_funcvaluese,natoms,0:max_num_neighbors,3) ! in
      real*8 dchargedsfunc(natoms,maxnum_funcvaluese)    ! in
      real*8 tempsum                                     ! internal
!!
      logical lperiodic                                  ! in
      logical ldoforces                                  ! in
!!
!! initializations
      ewaldcorr=0.0d0
      ewaldcorrforce(:,:)=0.0d0
!!
!! calculate neighbor list
      call neighbor_para(n_start,n_end,&
        num_atoms,zelem,&
        lsta,lstb,lstc,lste,&
        rscreen_cut,lattice,xyzstruct,lperiodic)
!!
!! loop over all atoms
      do i=1,natoms
!! loop over all neighbors
        do j=lsta(1,i),lsta(2,i)
          n=lstc(j)
          rij=lstb(j,4)
!!
          call getscreenfunctionforelectrostatics(&
            rij,fscreen,fscreenderiv,0.0d0)
!!
          ewaldcorr=ewaldcorr + nnatomcharge(atomindex(i))*nnatomcharge(n) &
            /rij*(1.0d0-fscreen)
!!
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
            do i1=1,num_atoms
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
!!
                call getscreenfunctionforelectrostatics(&
                  rij,fscreen,fscreenderiv,drijd(i2))
!!
                tempsum=0.0d0
                do i3=1,num_funcvaluese(elementindex(zelem(atomindex(i))))
                  tempsum=tempsum+dchargedsfunc(i,i3)*dsfuncdxyze(i3,i,invneighboridx(i,i1),i2)
                enddo

                ewaldcorrforce(i2,i1)=ewaldcorrforce(i2,i1)   &
                 -2.d0*(nnatomcharge(n)*tempsum/rij)*(1-fscreen) &
                 -(nnatomcharge(atomindex(i))*nnatomcharge(n)*drijd(i2)/(rij**2)*(-1)) &
                 *(1-fscreen)                                           &
                 +(nnatomcharge(atomindex(i))*nnatomcharge(n)/rij)*fscreenderiv
!!
              enddo ! i2
            enddo ! i1
          endif ! ldoforces
        enddo
      enddo
!!
!! remove double counting
      ewaldcorr=0.5d0*ewaldcorr
!! remove double counting also for erealforce
      ewaldcorrforce(:,:)=0.5d0*ewaldcorrforce(:,:)
!!
!!
      return
      end
