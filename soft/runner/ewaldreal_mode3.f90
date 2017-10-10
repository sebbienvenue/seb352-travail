!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: Calculate real-space part of standard Ewald summation, parallel case
!!
!! called by:
!! - getewald.f90
!!
      subroutine ewaldreal_mode3(n_start,n_end,max_num_neighbors_elec,&
        num_atoms,zelem,invneighboridx_elec,&
        natoms,atomindex,num_funcvalues_elec,&
        nnatomcharge,dchargedsfunc,dsfuncdxyze,&
        lattice,xyzstruct,erealforce,ereal,ldoforces_local,lperiodic)
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
      integer lsta(2,max_num_atoms)                       ! internal, pointer
      integer lstc(listdim)                               ! internal, neighbor label
      integer lste(listdim)                               ! internal
      integer n                                           ! internal
      integer n_start                                     ! in
      integer n_end                                       ! in
      integer num_atoms                                   ! in
!!
      integer i1,i2,i3,i4,i5                              ! internal
!!
      real*8 nnatomcharge(max_num_atoms)                  ! in
      real*8 erealforce(3,max_num_atoms)                  ! out
      real*8 ewaldcorrforce(3,max_num_atoms)              ! out
      real*8 ereal                                        ! out
      real*8 ewaldcorr                                    ! internal 
      real*8 sqrtpiinv
      parameter(sqrtpiinv=0.564189583d0)
      real*8 dsfuncdxyze(maxnum_funcvalues_elec,natoms,0:max_num_neighbors_elec,3) ! in
      real*8 dchargedsfunc(natoms,maxnum_funcvalues_elec)     ! in
      real*8 lstb(listdim,4)                              ! internal, x,y,z, rij of neighbor
      real*8 temp1                                        ! internal
      real*8 tempsum                                      ! internal
      real*8 rij                                          ! internal
      real*8 deltaxj                                      ! internal
      real*8 deltayj                                      ! internal
      real*8 deltazj                                      ! internal
      real*8 drijdxi                                      ! internal
      real*8 drijdyi                                      ! internal
      real*8 drijdzi                                      ! internal
      real*8 drijdxj                                      ! internal
      real*8 drijdyj                                      ! internal
      real*8 drijdzj                                      ! internal
      real*8 drijd(3)                                     ! internal
      real*8 fscreen                                      ! internal
      real*8 fscreenderiv                                 ! internal
      real*8 comperrfct                                   ! erfc(x)
      real*8 dcomperrfct                                  ! d erfc(x)/d x
      real*8 xyzstruct(3,max_num_atoms)                   ! in
      real*8 lattice(3,3)                                 ! in
!!
      logical ldoforces_local                             ! in
      logical lperiodic                                   ! in
!!
      ereal=0.0d0
      erealforce(:,:)=0.0d0
      ewaldcorr=0.0d0
      ewaldcorrforce(:,:)=0.0d0
!!

!! calculate neighbor list for natoms
!! Do not use lsta, lstc, lstb from subroutine prediction, because the cutoff for symfunctions is different from ewaldcutoff
      call neighbor_para(n_start,n_end,&
        num_atoms,zelem,&
        lsta,lstb,lstc,lste,&
        ewaldcutoff,lattice,xyzstruct,lperiodic)
!! using the neighbor list automatically prevents self-interaction with R_ij=0 (division!)
!
!! loop over all atoms
      do i1=1,natoms
!! loop over all neighbors
!        do i2=lsta(1,atomindex(i1)),lsta(2,atomindex(i1))
        do i2=lsta(1,i1),lsta(2,i1)

          n=lstc(i2)
          rij=lstb(i2,4)
!!
          if(lscreen)then
            call getscreenfunctionforelectrostatics(&
              rij,fscreen,fscreenderiv,0.0d0)
!!
            ewaldcorr=ewaldcorr+nnatomcharge(atomindex(i1))*nnatomcharge(n) &
              /rij*(1.0d0-fscreen)
          endif
!!
          temp1=ewaldalpha*rij
!!
          ereal=ereal+nnatomcharge(atomindex(i1))*nnatomcharge(n)*&
            comperrfct(temp1)/rij
!          write(ounit,'(a,2i5,4f14.8)')'erealpart ',i1,i2,nnatomcharge(atomindex(i1)),&
!            nnatomcharge(n),comperrfct(temp1),rij
!!
          if(ldoforces_local)then
            deltaxj=-1.d0*(xyzstruct(1,atomindex(i1))-lstb(i2,1))
            deltayj=-1.d0*(xyzstruct(2,atomindex(i1))-lstb(i2,2))
            deltazj=-1.d0*(xyzstruct(3,atomindex(i1))-lstb(i2,3))
            drijdxi=-deltaxj/rij
            drijdyi=-deltayj/rij
            drijdzi=-deltazj/rij
            drijdxj=-1.d0*drijdxi
            drijdyj=-1.d0*drijdyi
            drijdzj=-1.d0*drijdzi
!!
!! calculate derivative with respect to coordinate i4 of atom i3
            do i3=1,num_atoms
              if(i3.eq.atomindex(i1)) then
                drijd(1)=drijdxi
                drijd(2)=drijdyi
                drijd(3)=drijdzi
              elseif(i3.eq.n) then
                drijd(1)=drijdxj
                drijd(2)=drijdyj
                drijd(3)=drijdzj
              else
                drijd(:)=0.0d0
              endif
!              drijd(:)=-1.0d0*drijd(:) ! CHECK
!!
              do i4=1,3
!!           
                if(lscreen)then
                  call getscreenfunctionforelectrostatics(&
                    rij,fscreen,fscreenderiv,drijd(i4))
                endif
!!
!! calculate dchargedxyz(atomindex(i1),i3,i4) = tempsum
                tempsum=0.0d0
                if((nn_type_elec.ne.3).and.(nn_type_elec.ne.4))then
                  do i5=1,num_funcvalues_elec(elementindex(zelem(atomindex(i1))))
                    tempsum=tempsum+dchargedsfunc(i1,i5)*dsfuncdxyze(i5,i1,invneighboridx_elec(i1,i3),i4)
                  enddo
                endif
!!
                if(lscreen)then
                  ewaldcorrforce(i4,i3)=ewaldcorrforce(i4,i3)   &
                  -2.d0*(nnatomcharge(n)*tempsum/rij)*(1-fscreen) &
                  -(nnatomcharge(atomindex(i1))*nnatomcharge(n)*drijd(i4)/(rij**2)*(-1)) &
                  *(1-fscreen)                                           &
                  +(nnatomcharge(atomindex(i1))*nnatomcharge(n)/rij)*fscreenderiv
                endif
!!
                erealforce(i4,i3)=erealforce(i4,i3)&
                 -2.0d0*nnatomcharge(n)*tempsum*comperrfct(temp1)/rij
                if(drijd(i4).ne.0.0d0)then
                  erealforce(i4,i3)=erealforce(i4,i3)&
                    -nnatomcharge(atomindex(i1))*nnatomcharge(n)&
                    *(dcomperrfct(temp1)*ewaldalpha*drijd(i4)/rij&
                     -(comperrfct(temp1)/rij)*drijd(i4)/rij)
                endif
!                write(ounit,'(a,4i4,5f14.8)')'erealforce ',i1,i2,i3,i4,nnatomcharge(atomindex(i1)),&
!                  nnatomcharge(n),comperrfct(temp1),drijd(i4),rij
!                write(ounit,'(a,4i4,2f14.8)')'erealforce ',i1,i2,i3,i4,tempsum,erealforce(i4,i3)
              enddo ! i4
            enddo ! i3
          endif ! ldoforces
        enddo ! i2 lsta
      enddo ! i1 natoms
!!
!! remove double counting
      ereal=0.5d0*ereal
!! remove double counting also for erealforce
      erealforce(:,:)=0.5d0*erealforce(:,:)
!!
      if(lscreen)then
!! remove double counting
        ewaldcorr=0.5d0*ewaldcorr
!! remove double counting also for erealforce
        ewaldcorrforce(:,:)=0.5d0*ewaldcorrforce(:,:)
!!
        ereal=ereal-ewaldcorr
        erealforce(:,:)=erealforce(:,:)-ewaldcorrforce(:,:)
      endif
!!
      return
      end
