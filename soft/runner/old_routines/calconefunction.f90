!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calcfunctions (2x)
!! - calcpairfunctions 
!! - getallelectrostatic.f90
!! - getallshortforces.f90
!!
      subroutine calconefunction(ounit,&
           maxnum_funcvalues,num_funcvalues,&
           max_num_atoms,nelem,num_atoms,zelem,listdim,&
           elementindex,&
           function_type,symelement,lattice,&
           xyzstruct,symfunction,maxcutoff,rmin,&
           funccutoff,eta,rshift,lambda,zeta,dsfuncdxyz,strs,&
           lperiodic,ldoforces,ldostress,lrmin,ldebug)
!!
      implicit none
!!
      integer maxnum_funcvalues                              ! in
      integer num_funcvalues(nelem)                          ! in
      integer max_num_atoms                                  ! in
      integer num_atoms                                      ! in
      integer zelem(max_num_atoms)                           ! in
      integer function_type(maxnum_funcvalues,nelem)         ! in
      integer symelement(maxnum_funcvalues,2,nelem)          ! in
      integer nelem
      integer listdim
      integer lsta(2,max_num_atoms)
      integer lstc(listdim) ! identification of atom
      integer lste(listdim) ! nuclear charge of atom
      integer iindex                                         ! internal
      integer elementindex(102)
      integer ounit
      integer n,m
      integer i1,i2,i3,i4                                    ! internal
!!
      real*8 lattice(3,3)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 symfunction(maxnum_funcvalues,max_num_atoms)                   ! out
      real*8 dsfuncdxyz(maxnum_funcvalues,max_num_atoms,max_num_atoms,3)    ! out
      real*8 maxcutoff                                                   ! in
      real*8 lstb(listdim,4) ! xyz and r_ij
      real*8 rij
      real*8 rik
      real*8 rjk
      real*8 pi
      real*8 f 
      real*8 g 
      real*8 funccutoff(maxnum_funcvalues,nelem)            ! in
      real*8 eta(maxnum_funcvalues,nelem)                   ! in
      real*8 rshift(maxnum_funcvalues,nelem)                ! in
      real*8 lambda(maxnum_funcvalues,nelem)                ! in
      real*8 zeta(maxnum_funcvalues,nelem)                  ! in
      real*8 deltaxj,deltayj,deltazj
      real*8 deltaxk,deltayk,deltazk
      real*8 drijdxi, drijdyi, drijdzi
      real*8 drijdxj, drijdyj, drijdzj
      real*8 drijdxk, drijdyk, drijdzk
      real*8 drikdxi, drikdyi, drikdzi
      real*8 drikdxj, drikdyj, drikdzj
      real*8 drikdxk, drikdyk, drikdzk
      real*8 drjkdxi, drjkdyi, drjkdzi
      real*8 drjkdxj, drjkdyj, drjkdzj
      real*8 drjkdxk, drjkdyk, drjkdzk
      real*8 strs(3,3,maxnum_funcvalues,max_num_atoms)                  ! out
      real*8 temp1,temp2,temp3,temp4
      real*8 fcutij
      real*8 fcutik
      real*8 fcutjk
      real*8 dfcutijdxi,dfcutijdyi,dfcutijdzi
      real*8 dfcutijdxj,dfcutijdyj,dfcutijdzj
      real*8 dfcutijdxk,dfcutijdyk,dfcutijdzk
      real*8 dfcutikdxi,dfcutikdyi,dfcutikdzi
      real*8 dfcutikdxj,dfcutikdyj,dfcutikdzj
      real*8 dfcutikdxk,dfcutikdyk,dfcutikdzk
      real*8 dfcutjkdxi,dfcutjkdyi,dfcutjkdzi
      real*8 dfcutjkdxj,dfcutjkdyj,dfcutjkdzj
      real*8 dfcutjkdxk,dfcutjkdyk,dfcutjkdzk
      real*8 costheta
      real*8 dcosthetadxi,dcosthetadyi,dcosthetadzi
      real*8 dcosthetadxj,dcosthetadyj,dcosthetadzj
      real*8 dcosthetadxk,dcosthetadyk,dcosthetadzk
      real*8 dthetadxi,dthetadyi,dthetadzi
      real*8 dthetadxj,dthetadyj,dthetadzj
      real*8 dthetadxk,dthetadyk,dthetadzk
      real*8 expxyz
      real*8 dexpxyzdxi
      real*8 dexpxyzdyi
      real*8 dexpxyzdzi
      real*8 dexpxyzdxj
      real*8 dexpxyzdyj
      real*8 dexpxyzdzj
      real*8 dexpxyzdxk
      real*8 dexpxyzdyk
      real*8 dexpxyzdzk
      real*8 exptemp
      real*8 dexptempdxi
      real*8 dexptempdyi
      real*8 dexptempdzi
      real*8 dexptempdxj
      real*8 dexptempdyj
      real*8 dexptempdzj
      real*8 dexptempdxk
      real*8 dexptempdyk
      real*8 dexptempdzk
      real*8 dgdxi,dgdyi,dgdzi
      real*8 dgdxj,dgdyj,dgdzj
      real*8 dgdxk,dgdyk,dgdzk
      real*8 dfdxi,dfdyi,dfdzi
      real*8 dfdxj,dfdyj,dfdzj
      real*8 dfdxk,dfdyk,dfdzk
      real*8 rmin
      real*8 theta
      real*8 rad2deg

      logical lperiodic
      logical ldoforces
      logical ldostress
      logical lrmin                          ! in/out
      logical ldebug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!    initialization
      pi=3.141592654d0
      symfunction(:,:)   =0.0d0
      dsfuncdxyz(:,:,:,:)=0.0d0
      strs(:,:,:,:)      =0.0d0
      lsta(:,:)          =0
      lstb(:,:)          =0.0d0
      lstc(:)            =0
      lste(:)            =0
      rad2deg            =180.d0/pi
      lrmin              =.true.
!!
!! get neighbor list here
      call neighbor(max_num_atoms,&
              num_atoms,zelem,listdim,&
              lsta,lstb,lstc,lste,&
              maxcutoff,lattice,xyzstruct,lperiodic,ldebug)
!!
!! check for obviously wrong structures with too short bonds
      do i1=1,num_atoms
        do i2=lsta(1,i1),lsta(2,i1)
          if(lstb(i2,4).le.rmin)then
            lrmin=.false.
!!            write(ounit,*)'Error: rij lt rmin ',i1,i2,lstb(i2,4)
!!            write(ounit,*)xyzstruct(1,i1),lstb(i2,1)
!!            write(ounit,*)xyzstruct(2,i1),lstb(i2,2)
!!            write(ounit,*)xyzstruct(3,i1),lstb(i2,3)
!!            stop
          endif
        enddo
      enddo
!!
!! check for isolated atoms without neighbors in the structure
      do i1=1,num_atoms
        if(lsta(1,i1).eq.0)then
          write(ounit,*)'ERROR: atom ',i1,' has no neighbors'
          stop
        endif
      enddo
!! calculate the symmetry function values
      do i1=1,num_atoms ! loop over all atoms of the structure
      iindex=elementindex(zelem(i1))
      do i2=1,num_funcvalues(iindex) ! loop over all symmetry function types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(function_type(i2,iindex).eq.1)then ! radial function
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement(i2,1,iindex).eq.lste(i3))then
!!
          rij=lstb(i3,4)
          if(rij.le.funccutoff(i2,iindex))then
          n=lstc(i3)
!!
          symfunction(i2,i1)=symfunction(i2,i1)&
            +0.5d0*(dcos(pi*rij/funccutoff(i2,iindex))+1.d0) 
!!
          deltaxj=-1.d0*(xyzstruct(1,i1)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,i1)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,i1)-lstb(i3,3))
          drijdxi=-deltaxj/rij
          drijdyi=-deltayj/rij
          drijdzi=-deltazj/rij
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
!!
          if(ldoforces)then
          temp1=-0.5d0*dsin(pi*rij/funccutoff(i2,iindex))*pi/funccutoff(i2,iindex)
!! Calculation of derivatives for forces 
!! dsfunc/dx
          dsfuncdxyz(i2,i1,i1,1)=dsfuncdxyz(i2,i1,i1,1)&
            +(temp1*drijdxi)
          dsfuncdxyz(i2,i1,n,1) =dsfuncdxyz(i2,i1,n,1)&
            +(temp1*drijdxj)
!! dsfunc/dy
          dsfuncdxyz(i2,i1,i1,2)=dsfuncdxyz(i2,i1,i1,2)&
            +(temp1*drijdyi)
          dsfuncdxyz(i2,i1,n,2) =dsfuncdxyz(i2,i1,n,2)&
            +(temp1*drijdyj)
!! dsfunc/dz
          dsfuncdxyz(i2,i1,i1,3)=dsfuncdxyz(i2,i1,i1,3)&
            +(temp1*drijdzi)
          dsfuncdxyz(i2,i1,n,3) =dsfuncdxyz(i2,i1,n,3)&
            +(temp1*drijdzj)
!!
          if(ldostress)then
!! Calculation of derivatives for stress
!! dsfunc/dx
            strs(1,1,i2,i1)=strs(1,1,i2,i1)&
              +deltaxj*(drijdxj*temp1)
            strs(2,1,i2,i1)=strs(2,1,i2,i1)&
              +deltayj*(drijdxj*temp1)
            strs(3,1,i2,i1)=strs(3,1,i2,i1)&
              +deltazj*(drijdxj*temp1)
!! dsfunc/dy
            strs(1,2,i2,i1)=strs(1,2,i2,i1)&
              +deltaxj*(drijdyj*temp1)
            strs(2,2,i2,i1)=strs(2,2,i2,i1)&
              +deltayj*(drijdyj*temp1)
            strs(3,2,i2,i1)=strs(3,2,i2,i1)&
              +deltazj*(drijdyj*temp1)
!! dsfunc/dz
            strs(1,3,i2,i1)=strs(1,3,i2,i1)&
              +deltaxj*(drijdzj*temp1)
            strs(2,3,i2,i1)=strs(2,3,i2,i1)&
              +deltayj*(drijdzj*temp1)
            strs(3,3,i2,i1)=strs(3,3,i2,i1)&
              +deltazj*(drijdzj*temp1)
!!
          endif ! ldostress 
          endif ! ldoforces
          endif ! rij.le.funccutoff(i2,iindex)
          endif ! (symelement(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.2)then ! radial function
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement(i2,1,iindex).eq.lste(i3))then
          rij=lstb(i3,4)
          if(rij.le.funccutoff(i2,iindex))then
          n=lstc(i3)
!!
          deltaxj=-1.d0*(xyzstruct(1,i1)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,i1)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,i1)-lstb(i3,3))
          drijdxi=-deltaxj/rij
          drijdyi=-deltayj/rij
          drijdzi=-deltazj/rij
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
!!
          fcutij=0.5d0*(dcos(pi*rij/funccutoff(i2,iindex))+1.d0)
          temp1=0.5d0*(-dsin(pi*rij/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
          dfcutijdxi=temp1*drijdxi
          dfcutijdyi=temp1*drijdyi
          dfcutijdzi=temp1*drijdzi
          dfcutijdxj=-1.d0*dfcutijdxi
          dfcutijdyj=-1.d0*dfcutijdyi
          dfcutijdzj=-1.d0*dfcutijdzi
!!
          symfunction(i2,i1)=symfunction(i2,i1)&
            +dexp(-1.d0*eta(i2,iindex)*(rij-rshift(i2,iindex))**2)*fcutij 
!!
          if(ldoforces)then
!! Calculation of derivatives for forces 
            temp1=-2.d0*eta(i2,iindex)*(rij-rshift(i2,iindex))*dexp(-1.d0*eta(i2,iindex)&
                  *(rij-rshift(i2,iindex))**2)*fcutij
            temp2= dexp(-1.d0*eta(i2,iindex)*(rij-rshift(i2,iindex))**2)
!! dsfunc/dx
            dsfuncdxyz(i2,i1,i1,1)=dsfuncdxyz(i2,i1,i1,1)+&
               (drijdxi*temp1&
              + temp2*dfcutijdxi)
            dsfuncdxyz(i2,i1,n,1) =dsfuncdxyz(i2,i1,n,1)+&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
!! dsfunc/dy
            dsfuncdxyz(i2,i1,i1,2)=dsfuncdxyz(i2,i1,i1,2)+&
               (drijdyi*temp1&
              + temp2*dfcutijdyi)
            dsfuncdxyz(i2,i1,n,2) =dsfuncdxyz(i2,i1,n,2)+&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
!! dsfunc/dz
            dsfuncdxyz(i2,i1,i1,3)=dsfuncdxyz(i2,i1,i1,3)+&
               (drijdzi*temp1&
              + temp2*dfcutijdzi)
            dsfuncdxyz(i2,i1,n,3) =dsfuncdxyz(i2,i1,n,3)+&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
!!
          if(ldostress)then
!! Calculation of derivatives for stress
!! dsfunc/dx
            strs(1,1,i2,i1)=strs(1,1,i2,i1)+deltaxj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
            strs(2,1,i2,i1)=strs(2,1,i2,i1)+deltayj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
            strs(3,1,i2,i1)=strs(3,1,i2,i1)+deltazj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
!! dsfunc/dy
            strs(1,2,i2,i1)=strs(1,2,i2,i1)+deltaxj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
            strs(2,2,i2,i1)=strs(2,2,i2,i1)+deltayj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
            strs(3,2,i2,i1)=strs(3,2,i2,i1)+deltazj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
!! dsfunc/dz
            strs(1,3,i2,i1)=strs(1,3,i2,i1)+deltaxj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
            strs(2,3,i2,i1)=strs(2,3,i2,i1)+deltayj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
            strs(3,3,i2,i1)=strs(3,3,i2,i1)+deltazj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
!!
          endif ! ldostress 
          endif ! ldoforces
          endif ! rij.le.funccutoff(i2,iindex)
          endif ! (symelement(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.3)then ! angular function
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement(i2,1,iindex).eq.lste(i3))then
          rij=lstb(i3,4)
          n=lstc(i3)
          if(rij.le.funccutoff(i2,iindex))then
!!
          deltaxj=-1.d0*(xyzstruct(1,i1)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,i1)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,i1)-lstb(i3,3))
          drijdxi=-deltaxj/rij
          drijdyi=-deltayj/rij
          drijdzi=-deltazj/rij
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
          drijdxk=0.0d0
          drijdyk=0.0d0
          drijdzk=0.0d0
!!
          fcutij=0.5d0*(dcos(pi*rij/funccutoff(i2,iindex))+1.d0)
          temp1=0.5d0*(-dsin(pi*rij/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
          dfcutijdxi=temp1*drijdxi
          dfcutijdyi=temp1*drijdyi
          dfcutijdzi=temp1*drijdzi
          dfcutijdxj=-1.d0*dfcutijdxi
          dfcutijdyj=-1.d0*dfcutijdyi
          dfcutijdzj=-1.d0*dfcutijdzi
          dfcutijdxk=0.0d0
          dfcutijdyk=0.0d0
          dfcutijdzk=0.0d0
!!
          do i4=lsta(1,i1),lsta(2,i1) ! loop over all neighbors of this atom
            if(symelement(i2,2,iindex).eq.lste(i4))then
!! new double counting criterion JB and TM 08.03.2010:
            if(((i4.gt.i3).and.(lste(i3).eq.lste(i4)))&
             .or.(lste(i3).ne.lste(i4)))then ! avoid double counting
!! This was a bug, because the order of neighbor atoms is not arbitrary
!!            if(i4.gt.i3)then ! avoid double counting
!!
            m=lstc(i4)
            rik=lstb(i4,4)
!!
            if(rik.le.funccutoff(i2,iindex))then
!!
              deltaxk=-1.d0*(xyzstruct(1,i1)-lstb(i4,1))
              deltayk=-1.d0*(xyzstruct(2,i1)-lstb(i4,2))
              deltazk=-1.d0*(xyzstruct(3,i1)-lstb(i4,3))
              drikdxi=-deltaxk/rik
              drikdyi=-deltayk/rik
              drikdzi=-deltazk/rik
              drikdxk=-1.d0*drikdxi
              drikdyk=-1.d0*drikdyi
              drikdzk=-1.d0*drikdzi
              drikdxj=0.0d0
              drikdyj=0.0d0
              drikdzj=0.0d0
!!
              fcutik=0.5d0*(dcos(pi*rik/funccutoff(i2,iindex))+1.d0)
              temp1=0.5d0*(-dsin(pi*rik/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
              dfcutikdxi=temp1*drikdxi
              dfcutikdyi=temp1*drikdyi
              dfcutikdzi=temp1*drikdzi
              dfcutikdxj=0.0d0
              dfcutikdyj=0.0d0
              dfcutikdzj=0.0d0
              dfcutikdxk=-1.d0*dfcutikdxi
              dfcutikdyk=-1.d0*dfcutikdyi
              dfcutikdzk=-1.d0*dfcutikdzi
!!
!! calculate the third bond length between atom i3 and i4
!!
              rjk=(lstb(i3,1)-lstb(i4,1))**2 +&
                  (lstb(i3,2)-lstb(i4,2))**2 +&
                  (lstb(i3,3)-lstb(i4,3))**2
              rjk=dsqrt(rjk)
              if(rjk.le.rmin) then
                lrmin=.false.
!!                write(ounit,*)'Error rjk .le. rmin ',i2,rjk
!!                stop
              endif
!!
              if(rjk.le.funccutoff(i2,iindex)) then
                drjkdxj=(lstb(i3,1)-lstb(i4,1))/rjk
                drjkdyj=(lstb(i3,2)-lstb(i4,2))/rjk
                drjkdzj=(lstb(i3,3)-lstb(i4,3))/rjk
                drjkdxk=-1.d0*drjkdxj
                drjkdyk=-1.d0*drjkdyj
                drjkdzk=-1.d0*drjkdzj
                drjkdxi=0.0d0
                drjkdyi=0.0d0
                drjkdzi=0.0d0
!!
                fcutjk=0.5d0*(dcos(pi*rjk/funccutoff(i2,iindex))+1.d0)
                temp1=0.5d0*(-dsin(pi*rjk/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
                dfcutjkdxj=temp1*drjkdxj
                dfcutjkdyj=temp1*drjkdyj
                dfcutjkdzj=temp1*drjkdzj
                dfcutjkdxk=-1.d0*dfcutjkdxj
                dfcutjkdyk=-1.d0*dfcutjkdyj
                dfcutjkdzk=-1.d0*dfcutjkdzj
                dfcutjkdxi=0.0d0
                dfcutjkdyi=0.0d0
                dfcutjkdzi=0.0d0
!!
!! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
                f=lambda(i2,iindex)*(rjk**2 - rij**2 -rik**2)
                g=-2.d0*rij*rik
                costheta=f/g
                costheta=1.d0+costheta ! avoid negative values
                temp4   =2.d0**(1.d0-zeta(i2,iindex))
                costheta=temp4*(costheta**zeta(i2,iindex))
!!
!! calculate the derivatives of costheta already here
!! (f/g)' = (f'g - fg')/g^2
                if(ldoforces) then
                 dfdxi=lambda(i2,iindex)*(-2.d0*rij*drijdxi - 2.d0*rik*drikdxi)
                 dfdyi=lambda(i2,iindex)*(-2.d0*rij*drijdyi - 2.d0*rik*drikdyi)
                 dfdzi=lambda(i2,iindex)*(-2.d0*rij*drijdzi - 2.d0*rik*drikdzi)
!!
                 dfdxj=lambda(i2,iindex)*(2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj)
                 dfdyj=lambda(i2,iindex)*(2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj)
                 dfdzj=lambda(i2,iindex)*(2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj)
!!
                 dfdxk=lambda(i2,iindex)*(2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk)
                 dfdyk=lambda(i2,iindex)*(2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk)
                 dfdzk=lambda(i2,iindex)*(2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk)
!!
                 dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                 dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                 dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)
!!
                 dgdxj=-2.d0*drijdxj*rik
                 dgdyj=-2.d0*drijdyj*rik
                 dgdzj=-2.d0*drijdzj*rik
!!
                 dgdxk=-2.d0*rij*drikdxk
                 dgdyk=-2.d0*rij*drikdyk
                 dgdzk=-2.d0*rij*drikdzk
!!
                 temp4=temp4*zeta(i2,iindex)*(1.d0+f/g)**(zeta(i2,iindex)-1)/g**2
                 dcosthetadxi=temp4*(dfdxi*g - f*dgdxi) !/g**2
                 dcosthetadyi=temp4*(dfdyi*g - f*dgdyi) !/g**2
                 dcosthetadzi=temp4*(dfdzi*g - f*dgdzi) !/g**2
                 dcosthetadxj=temp4*(dfdxj*g - f*dgdxj) !/g**2
                 dcosthetadyj=temp4*(dfdyj*g - f*dgdyj) !/g**2
                 dcosthetadzj=temp4*(dfdzj*g - f*dgdzj) !/g**2
                 dcosthetadxk=temp4*(dfdxk*g - f*dgdxk) !/g**2
                 dcosthetadyk=temp4*(dfdyk*g - f*dgdyk) !/g**2
                 dcosthetadzk=temp4*(dfdzk*g - f*dgdzk) !/g**2
                endif ! ldoforces for derivative of costheta
!!
!! calculation of exponential term
                expxyz=dexp(-eta(i2,iindex)*(rij**2+rik**2+rjk**2))
                temp1=-eta(i2,iindex)*2.0d0*expxyz
                dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
                dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
                dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
                dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
                dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
                dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
                dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
                dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
                dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1
!!
!! calculation of the symmetry function
          symfunction(i2,i1)=symfunction(i2,i1)&
                                 +costheta*expxyz*fcutij*fcutik*fcutjk
!!
              if(ldoforces)then
!! Calculation of derivatives for forces 
!! (fgh)' = f'gh + fg'h + fgh'
!! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'
!! for the derivatives fcutjk is a constant
!!
!! dsfunc/dx
               temp1=(+dcosthetadxi* expxyz   * fcutij   * fcutik   * fcutjk&
                       +costheta   *dexpxyzdxi* fcutij   * fcutik   * fcutjk&
                       +costheta   * expxyz   *dfcutijdxi* fcutik   * fcutjk&
                       +costheta   * expxyz   * fcutij   *dfcutikdxi* fcutjk&
                       +costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdxi)
               temp2=(+dcosthetadxj* expxyz   * fcutij   * fcutik   * fcutjk&
                       +costheta   *dexpxyzdxj* fcutij   * fcutik   * fcutjk&
                       +costheta   * expxyz   *dfcutijdxj* fcutik   * fcutjk&
                       +costheta   * expxyz   * fcutij   *dfcutikdxj* fcutjk&
                       +costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdxj)
               temp3=(+dcosthetadxk* expxyz   * fcutij   * fcutik   * fcutjk&
                       +costheta   *dexpxyzdxk* fcutij   * fcutik   * fcutjk&
                       +costheta   * expxyz   *dfcutijdxk* fcutik   * fcutjk&
                       +costheta   * expxyz   * fcutij   *dfcutikdxk* fcutjk&
                       +costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdxk)
            dsfuncdxyz(i2,i1,i1,1)=dsfuncdxyz(i2,i1,i1,1)+&
                                       temp1 
            dsfuncdxyz(i2,i1,n,1) =dsfuncdxyz(i2,i1,n,1)+&
                                       temp2 
            dsfuncdxyz(i2,i1,m,1) =dsfuncdxyz(i2,i1,m,1)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,1,i2,i1)=strs(1,1,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,1,i2,i1)=strs(2,1,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,1,i2,i1)=strs(3,1,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
!!
!! dsfunc/dy
               temp1=(+dcosthetadyi* expxyz   * fcutij   * fcutik   * fcutjk&
                       +costheta   *dexpxyzdyi* fcutij   * fcutik   * fcutjk&
                       +costheta   * expxyz   *dfcutijdyi* fcutik   * fcutjk&
                       +costheta   * expxyz   * fcutij   *dfcutikdyi* fcutjk&
                       +costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdyi)
               temp2=(+dcosthetadyj* expxyz   * fcutij   * fcutik   * fcutjk&
                       +costheta   *dexpxyzdyj* fcutij   * fcutik   * fcutjk&
                       +costheta   * expxyz   *dfcutijdyj* fcutik   * fcutjk&
                       +costheta   * expxyz   * fcutij   *dfcutikdyj* fcutjk&
                       +costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdyj)
               temp3=(+dcosthetadyk* expxyz   * fcutij   * fcutik   * fcutjk&
                       +costheta   *dexpxyzdyk* fcutij   * fcutik   * fcutjk&
                       +costheta   * expxyz   *dfcutijdyk* fcutik   * fcutjk&
                       +costheta   * expxyz   * fcutij   *dfcutikdyk* fcutjk&
                       +costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdyk)
            dsfuncdxyz(i2,i1,i1,2)=dsfuncdxyz(i2,i1,i1,2)+&
                                       temp1 
            dsfuncdxyz(i2,i1,n,2) =dsfuncdxyz(i2,i1,n,2)+&
                                       temp2 
            dsfuncdxyz(i2,i1,m,2) =dsfuncdxyz(i2,i1,m,2)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,2,i2,i1)=strs(1,2,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,2,i2,i1)=strs(2,2,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,2,i2,i1)=strs(3,2,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
!!
!! dsfunc/dz
               temp1=(+dcosthetadzi* expxyz   * fcutij   * fcutik   * fcutjk&
                      + costheta   *dexpxyzdzi* fcutij   * fcutik   * fcutjk&
                      + costheta   * expxyz   *dfcutijdzi* fcutik   * fcutjk&
                      + costheta   * expxyz   * fcutij   *dfcutikdzi* fcutjk&
                      + costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdzi)
               temp2=(+dcosthetadzj* expxyz   * fcutij   * fcutik   * fcutjk&
                      + costheta   *dexpxyzdzj* fcutij   * fcutik   * fcutjk&
                      + costheta   * expxyz   *dfcutijdzj* fcutik   * fcutjk&
                      + costheta   * expxyz   * fcutij   *dfcutikdzj* fcutjk&
                      + costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdzj)
               temp3=(+dcosthetadzk* expxyz   * fcutij   * fcutik   * fcutjk&
                      + costheta   *dexpxyzdzk* fcutij   * fcutik   * fcutjk&
                      + costheta   * expxyz   *dfcutijdzk* fcutik   * fcutjk&
                      + costheta   * expxyz   * fcutij   *dfcutikdzk* fcutjk&
                      + costheta   * expxyz   * fcutij   * fcutik   *dfcutjkdzk)
            dsfuncdxyz(i2,i1,i1,3)=dsfuncdxyz(i2,i1,i1,3)+&
                                       temp1 
            dsfuncdxyz(i2,i1,n,3) =dsfuncdxyz(i2,i1,n,3)+&
                                       temp2 
            dsfuncdxyz(i2,i1,m,3) =dsfuncdxyz(i2,i1,m,3)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,3,i2,i1)=strs(1,3,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,3,i2,i1)=strs(2,3,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,3,i2,i1)=strs(3,3,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
              endif ! ldoforces
!!
              endif ! rjk .le. cutoff
            endif ! rik .le. cutoff
!!
          endif ! i4.gt.i3
          endif ! (symelement(i2,2,iindex).eq.lste(i4))then
          enddo ! i4
          endif ! rij.le.funccutoff(i2,iindex)
          endif ! (symelement(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.4)then ! radial function
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement(i2,1,iindex).eq.lste(i3))then
          rij=lstb(i3,4)
          if(rij.le.funccutoff(i2,iindex))then
          n=lstc(i3)
!!
          deltaxj=-1.d0*(xyzstruct(1,i1)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,i1)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,i1)-lstb(i3,3))
          drijdxi=-deltaxj/rij
          drijdyi=-deltayj/rij
          drijdzi=-deltazj/rij
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
!!
          fcutij=0.5d0*(dcos(pi*rij/funccutoff(i2,iindex))+1.d0)
          temp1=0.5d0*(-dsin(pi*rij/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
          dfcutijdxi=temp1*drijdxi
          dfcutijdyi=temp1*drijdyi
          dfcutijdzi=temp1*drijdzi
          dfcutijdxj=-1.d0*dfcutijdxi
          dfcutijdyj=-1.d0*dfcutijdyi
          dfcutijdzj=-1.d0*dfcutijdzi
!!
          symfunction(i2,i1)=symfunction(i2,i1)&
            +dcos(eta(i2,iindex)*rij)*fcutij 
!!
          if(ldoforces)then
!! Calculation of derivatives for forces 
            temp1=-1.d0*eta(i2,iindex)*dsin(eta(i2,iindex)*rij)*fcutij
            temp2= dcos(eta(i2,iindex)*rij)
!! dsfunc/dx
            dsfuncdxyz(i2,i1,i1,1)=dsfuncdxyz(i2,i1,i1,1)+&
               (drijdxi*temp1&
              + temp2*dfcutijdxi)
            dsfuncdxyz(i2,i1,n,1) =dsfuncdxyz(i2,i1,n,1)+&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
!! dsfunc/dy
            dsfuncdxyz(i2,i1,i1,2)=dsfuncdxyz(i2,i1,i1,2)+&
               (drijdyi*temp1&
              + temp2*dfcutijdyi)
            dsfuncdxyz(i2,i1,n,2) =dsfuncdxyz(i2,i1,n,2)+&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
!! dsfunc/dz
            dsfuncdxyz(i2,i1,i1,3)=dsfuncdxyz(i2,i1,i1,3)+&
               (drijdzi*temp1&
              + temp2*dfcutijdzi)
            dsfuncdxyz(i2,i1,n,3) =dsfuncdxyz(i2,i1,n,3)+&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
!!
          if(ldostress)then
!! Calculation of derivatives for stress
!! dsfunc/dx
            strs(1,1,i2,i1)=strs(1,1,i2,i1)+deltaxj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
            strs(2,1,i2,i1)=strs(2,1,i2,i1)+deltayj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
            strs(3,1,i2,i1)=strs(3,1,i2,i1)+deltazj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
!! dsfunc/dy
            strs(1,2,i2,i1)=strs(1,2,i2,i1)+deltaxj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
            strs(2,2,i2,i1)=strs(2,2,i2,i1)+deltayj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
            strs(3,2,i2,i1)=strs(3,2,i2,i1)+deltazj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
!! dsfunc/dz
            strs(1,3,i2,i1)=strs(1,3,i2,i1)+deltaxj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
            strs(2,3,i2,i1)=strs(2,3,i2,i1)+deltayj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
            strs(3,3,i2,i1)=strs(3,3,i2,i1)+deltazj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
!!
          endif ! ldostress 
          endif ! ldoforces
          endif ! rij.le.funccutoff(i2,iindex)
          endif !(symelement(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.5)then ! just Cartesian coordinate

          if(eta(i2,iindex).eq.1.0d0)then
            symfunction(i2,i1)=xyzstruct(1,i1)
          elseif(eta(i2,iindex).eq.2.0d0)then
            symfunction(i2,i1)=xyzstruct(2,i1)
          elseif(eta(i2,iindex).eq.3.0d0)then
            symfunction(i2,i1)=xyzstruct(3,i1)
          else
            write(ounit,*)'Error: undefined value for symfunction type 5 ',eta(i2,iindex)
            stop !'
          endif

          if(ldoforces)then
!! Calculation of derivatives for forces
!! dsfunc/dx
            if(eta(i2,iindex).eq.1.0d0)then
              dsfuncdxyz(i2,i1,i1,1)=1.0d0
            else
              dsfuncdxyz(i2,i1,i1,1)=0.0d0
            endif
!! dsfunc/dy
            if(eta(i2,iindex).eq.2.0d0)then
              dsfuncdxyz(i2,i1,i1,2)=1.0d0
            else
              dsfuncdxyz(i2,i1,i1,2)=0.0d0
            endif
!! dsfunc/dz
            if(eta(i2,iindex).eq.3.0d0)then
              dsfuncdxyz(i2,i1,i1,3)=1.0d0
            else
              dsfuncdxyz(i2,i1,i1,3)=0.0d0
            endif
!!
          if(ldostress)then
!! Calculation of derivatives for stress
            write(ounit,*)'Error: stress is not implemented for symfunction type 5'
            stop !'
!!
          endif ! ldostress
          endif ! ldoforces
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.6)then ! radial function
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement(i2,1,iindex).eq.lste(i3))then
          if(lsta(1,i1).ne.lsta(2,i1))then
            write(ounit,*)'ERROR: symmetry function type 6 is only for dimers'
            stop !'
          endif
          rij=lstb(i3,4)
          if(rij.le.funccutoff(i2,iindex))then
          n=lstc(i3)
!!
          deltaxj=-1.d0*(xyzstruct(1,i1)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,i1)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,i1)-lstb(i3,3))
          drijdxi=-deltaxj/rij
          drijdyi=-deltayj/rij
          drijdzi=-deltazj/rij
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
!!
          symfunction(i2,i1)=rij
!!
          if(ldoforces)then
!! Calculation of derivatives for forces
!! dsfunc/dx
            dsfuncdxyz(i2,i1,i1,1)=dsfuncdxyz(i2,i1,i1,1)+&
               drijdxi
            dsfuncdxyz(i2,i1,n,1) =dsfuncdxyz(i2,i1,n,1)+&
               drijdxj
!! dsfunc/dy
            dsfuncdxyz(i2,i1,i1,2)=dsfuncdxyz(i2,i1,i1,2)+&
               drijdyi
            dsfuncdxyz(i2,i1,n,2) =dsfuncdxyz(i2,i1,n,2)+&
               drijdyj
!! dsfunc/dz
            dsfuncdxyz(i2,i1,i1,3)=dsfuncdxyz(i2,i1,i1,3)+&
               drijdzi
            dsfuncdxyz(i2,i1,n,3) =dsfuncdxyz(i2,i1,n,3)+&
               drijdzj
!!
          if(ldostress)then
!!
            write(ounit,*)'ERROR: no stress implemented for symmetry function type 6'
            stop !'
!!
          endif ! ldostress
          endif ! ldoforces
          endif ! rij.le.funccutoff(i2,iindex)
          endif !(symelement(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.7)then ! angular function
        write(ounit,*)'Error: function type not implemented in calconefunction ',function_type(i1,iindex)
        stop !'
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.8)then ! angular function
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement(i2,1,iindex).eq.lste(i3))then
          rij=lstb(i3,4)
          n=lstc(i3)
          if(rij.le.funccutoff(i2,iindex))then
!!
          deltaxj=-1.d0*(xyzstruct(1,i1)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,i1)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,i1)-lstb(i3,3))
          drijdxi=-deltaxj/rij
          drijdyi=-deltayj/rij
          drijdzi=-deltazj/rij
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
          drijdxk=0.0d0
          drijdyk=0.0d0
          drijdzk=0.0d0
!!
          fcutij=0.5d0*(dcos(pi*rij/funccutoff(i2,iindex))+1.d0)
          temp1=0.5d0*(-dsin(pi*rij/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
          dfcutijdxi=temp1*drijdxi
          dfcutijdyi=temp1*drijdyi
          dfcutijdzi=temp1*drijdzi
          dfcutijdxj=-1.d0*dfcutijdxi
          dfcutijdyj=-1.d0*dfcutijdyi
          dfcutijdzj=-1.d0*dfcutijdzi
          dfcutijdxk=0.0d0
          dfcutijdyk=0.0d0
          dfcutijdzk=0.0d0
!!
          do i4=lsta(1,i1),lsta(2,i1) ! loop over all neighbors of this atom
            if(symelement(i2,2,iindex).eq.lste(i4))then
!! new double counting criterion JB and TM 08.03.2010:
            if(((i4.gt.i3).and.(lste(i3).eq.lste(i4)))&
             .or.(lste(i3).ne.lste(i4)))then ! avoid double counting
!!
            m=lstc(i4)
            rik=lstb(i4,4)
!!
            if(rik.le.funccutoff(i2,iindex))then
!!
              deltaxk=-1.d0*(xyzstruct(1,i1)-lstb(i4,1))
              deltayk=-1.d0*(xyzstruct(2,i1)-lstb(i4,2))
              deltazk=-1.d0*(xyzstruct(3,i1)-lstb(i4,3))
              drikdxi=-deltaxk/rik
              drikdyi=-deltayk/rik
              drikdzi=-deltazk/rik
              drikdxk=-1.d0*drikdxi
              drikdyk=-1.d0*drikdyi
              drikdzk=-1.d0*drikdzi
              drikdxj=0.0d0
              drikdyj=0.0d0
              drikdzj=0.0d0
!!
              fcutik=0.5d0*(dcos(pi*rik/funccutoff(i2,iindex))+1.d0)
              temp1=0.5d0*(-dsin(pi*rik/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
              dfcutikdxi=temp1*drikdxi
              dfcutikdyi=temp1*drikdyi
              dfcutikdzi=temp1*drikdzi
              dfcutikdxj=0.0d0
              dfcutikdyj=0.0d0
              dfcutikdzj=0.0d0
              dfcutikdxk=-1.d0*dfcutikdxi
              dfcutikdyk=-1.d0*dfcutikdyi
              dfcutikdzk=-1.d0*dfcutikdzi
!!
!! calculate the third bond length between atom i3 and i4
!!
              rjk=(lstb(i3,1)-lstb(i4,1))**2 +&
                  (lstb(i3,2)-lstb(i4,2))**2 +&
                  (lstb(i3,3)-lstb(i4,3))**2
              rjk=dsqrt(rjk)
              if(rjk.le.rmin) then
                lrmin=.false.
              endif
!!
              if(rjk.le.funccutoff(i2,iindex)) then
                drjkdxj=(lstb(i3,1)-lstb(i4,1))/rjk
                drjkdyj=(lstb(i3,2)-lstb(i4,2))/rjk
                drjkdzj=(lstb(i3,3)-lstb(i4,3))/rjk
                drjkdxk=-1.d0*drjkdxj
                drjkdyk=-1.d0*drjkdyj
                drjkdzk=-1.d0*drjkdzj
                drjkdxi=0.0d0
                drjkdyi=0.0d0
                drjkdzi=0.0d0
!!
                fcutjk=0.5d0*(dcos(pi*rjk/funccutoff(i2,iindex))+1.d0)
                temp1=0.5d0*(-dsin(pi*rjk/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
                dfcutjkdxj=temp1*drjkdxj
                dfcutjkdyj=temp1*drjkdyj
                dfcutjkdzj=temp1*drjkdzj
                dfcutjkdxk=-1.d0*dfcutjkdxj
                dfcutjkdyk=-1.d0*dfcutjkdyj
                dfcutjkdzk=-1.d0*dfcutjkdzj
                dfcutjkdxi=0.0d0
                dfcutjkdyi=0.0d0
                dfcutjkdzi=0.0d0
!!
!! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
                f=rjk**2-rij**2-rik**2
                g=-2.d0*rij*rik
!!                costheta=f/g
!! filter out numerical noise that results in NaN for acos(x):
                temp1=f/g
                if(temp1.gt.1.0d0)then
                  if(temp1.gt.1.1d0)then
                    write(ounit,*)'Error in calconefunction type 8 ',temp1
                    stop
                  endif
                  temp1=1.0d0
                endif
                if(temp1.lt.-1.0d0)then
                  if(temp1.lt.-1.1d0)then
                    write(ounit,*)'Error in calconefunction type 8 ',temp1
                    stop
                  endif
                  temp1=-1.0d0
                endif
                theta=acos(temp1)
!! convert to degree
                theta=theta*rad2deg
!!
!! calculate the derivatives of costheta already here
!! (f/g)' = (f'g - fg')/g^2
                if(ldoforces) then
                 dfdxi=lambda(i2,iindex)*(-2.d0*rij*drijdxi - 2.d0*rik*drikdxi)
                 dfdyi=lambda(i2,iindex)*(-2.d0*rij*drijdyi - 2.d0*rik*drikdyi)
                 dfdzi=lambda(i2,iindex)*(-2.d0*rij*drijdzi - 2.d0*rik*drikdzi)
!!
                 dfdxj=lambda(i2,iindex)*(2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj)
                 dfdyj=lambda(i2,iindex)*(2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj)
                 dfdzj=lambda(i2,iindex)*(2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj)
!!
                 dfdxk=lambda(i2,iindex)*(2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk)
                 dfdyk=lambda(i2,iindex)*(2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk)
                 dfdzk=lambda(i2,iindex)*(2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk)
!!
                 dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                 dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                 dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)
!!
                 dgdxj=-2.d0*drijdxj*rik
                 dgdyj=-2.d0*drijdyj*rik
                 dgdzj=-2.d0*drijdzj*rik
!!
                 dgdxk=-2.d0*rij*drikdxk
                 dgdyk=-2.d0*rij*drikdyk
                 dgdzk=-2.d0*rij*drikdzk
!!
!! derivative of acos:
!! f(x)=acos(x)
!! f'(x)=-1/sqrt(1-x**2)
!!
!! filter out numerical noise that results in Infty for temp4:
!!                 temp4=-1.d0/dsqrt(1.d0-f**2/g**2)  ! this line has problems with numerical noise
                 temp4=f**2/g**2
!!                 if(temp4.gt.1.0d0)temp4=1.0d0
                 temp4=-1.d0/dsqrt(1.d0-temp4)
                 dthetadxi=rad2deg*temp4*(dfdxi*g - f*dgdxi)/g**2
                 dthetadyi=rad2deg*temp4*(dfdyi*g - f*dgdyi)/g**2
                 dthetadzi=rad2deg*temp4*(dfdzi*g - f*dgdzi)/g**2
                 dthetadxj=rad2deg*temp4*(dfdxj*g - f*dgdxj)/g**2
                 dthetadyj=rad2deg*temp4*(dfdyj*g - f*dgdyj)/g**2
                 dthetadzj=rad2deg*temp4*(dfdzj*g - f*dgdzj)/g**2
                 dthetadxk=rad2deg*temp4*(dfdxk*g - f*dgdxk)/g**2
                 dthetadyk=rad2deg*temp4*(dfdyk*g - f*dgdyk)/g**2
                 dthetadzk=rad2deg*temp4*(dfdzk*g - f*dgdzk)/g**2
!!
                endif ! ldoforces for derivative of costheta
!!
!! calculation of the symmetry function
                temp1=theta-rshift(i2,iindex)
                temp2=theta+rshift(i2,iindex)
                exptemp=(dexp(-eta(i2,iindex)*(temp1       )**2)&
                        +dexp(-eta(i2,iindex)*(temp2-360.d0)**2)&
                        +dexp(-eta(i2,iindex)*(temp2       )**2)&
                        +dexp(-eta(i2,iindex)*(temp1-360.d0)**2))
          symfunction(i2,i1)=symfunction(i2,i1)&
                  +exptemp*fcutij*fcutik*fcutjk
!!
              if(ldoforces)then
!! Calculation of derivatives for forces
!! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'
                temp3=-2.d0*eta(i2,iindex)*temp1*dexp(-eta(i2,iindex)*(temp1       )**2)&
                      -2.d0*eta(i2,iindex)*temp2*dexp(-eta(i2,iindex)*(temp2-360.d0)**2)&
                      -2.d0*eta(i2,iindex)*temp2*dexp(-eta(i2,iindex)*(temp2       )**2)&
                      -2.d0*eta(i2,iindex)*temp1*dexp(-eta(i2,iindex)*(temp1-360.d0)**2)
!!
                dexptempdxi=temp3*dthetadxi
                dexptempdyi=temp3*dthetadyi
                dexptempdzi=temp3*dthetadzi
                dexptempdxj=temp3*dthetadxj
                dexptempdyj=temp3*dthetadyj
                dexptempdzj=temp3*dthetadzj
                dexptempdxk=temp3*dthetadxk
                dexptempdyk=temp3*dthetadyk
                dexptempdzk=temp3*dthetadzk
!!
!! dsfunc/dx
               temp1=( dexptempdxi* fcutij   * fcutik   * fcutjk&
                      + exptemp   *dfcutijdxi* fcutik   * fcutjk&
                      + exptemp   * fcutij   *dfcutikdxi* fcutjk&
                      + exptemp   * fcutij   * fcutik   *dfcutjkdxi)
               temp2=( dexptempdxj* fcutij   * fcutik   * fcutjk&
                       +exptemp   *dfcutijdxj* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdxj* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdxj)
               temp3=( dexptempdxk* fcutij   * fcutik   * fcutjk&
                       +exptemp   *dfcutijdxk* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdxk* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdxk)
            dsfuncdxyz(i2,i1,i1,1)=dsfuncdxyz(i2,i1,i1,1)+&
                                       temp1
            dsfuncdxyz(i2,i1,n,1) =dsfuncdxyz(i2,i1,n,1)+&
                                       temp2
            dsfuncdxyz(i2,i1,m,1) =dsfuncdxyz(i2,i1,m,1)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,1,i2,i1)=strs(1,1,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,1,i2,i1)=strs(2,1,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,1,i2,i1)=strs(3,1,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress
!!
!! dsfunc/dy
               temp1=( dexptempdyi* fcutij   * fcutik   * fcutjk&
                       +exptemp   *dfcutijdyi* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdyi* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdyi)
               temp2=( dexptempdyj* fcutij   * fcutik   * fcutjk&
                       +exptemp   *dfcutijdyj* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdyj* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdyj)
               temp3=( dexptempdyk* fcutij   * fcutik   * fcutjk&
                       +exptemp   *dfcutijdyk* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdyk* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdyk)
            dsfuncdxyz(i2,i1,i1,2)=dsfuncdxyz(i2,i1,i1,2)+&
                                       temp1
            dsfuncdxyz(i2,i1,n,2) =dsfuncdxyz(i2,i1,n,2)+&
                                       temp2
            dsfuncdxyz(i2,i1,m,2) =dsfuncdxyz(i2,i1,m,2)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,2,i2,i1)=strs(1,2,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,2,i2,i1)=strs(2,2,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,2,i2,i1)=strs(3,2,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress
!!
!! dsfunc/dz
               temp1=( dexptempdzi* fcutij   * fcutik   * fcutjk&
                       +exptemp   *dfcutijdzi* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdzi* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdzi)
               temp2=( dexptempdzj   * fcutij* fcutik   * fcutjk&
                       +exptemp   *dfcutijdzj* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdzj* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdzj)
               temp3=( dexptempdzk* fcutij   * fcutik   * fcutjk&
                       +exptemp   *dfcutijdzk* fcutik   * fcutjk&
                       +exptemp   * fcutij   *dfcutikdzk* fcutjk&
                       +exptemp   * fcutij   * fcutik   *dfcutjkdzk)
            dsfuncdxyz(i2,i1,i1,3)=dsfuncdxyz(i2,i1,i1,3)+&
                                       temp1
            dsfuncdxyz(i2,i1,n,3) =dsfuncdxyz(i2,i1,n,3)+&
                                       temp2
            dsfuncdxyz(i2,i1,m,3) =dsfuncdxyz(i2,i1,m,3)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,3,i2,i1)=strs(1,3,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,3,i2,i1)=strs(2,3,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,3,i2,i1)=strs(3,3,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress
              endif ! ldoforces
!!
              endif ! rjk .le. cutoff
            endif ! rik .le. cutoff
!!
          endif ! i4.gt.i3
          endif !(symelement(i2,2,iindex).eq.lste(i4))then
          enddo ! i4
          endif ! rij.le.funccutoff(i2,iindex)
          endif ! (symelement(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type(i2,iindex).eq.9)then ! angular function
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement(i2,1,iindex).eq.lste(i3))then
          rij=lstb(i3,4)
          n=lstc(i3)
          if(rij.le.funccutoff(i2,iindex))then
!!
          deltaxj=-1.d0*(xyzstruct(1,i1)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,i1)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,i1)-lstb(i3,3))
          drijdxi=-deltaxj/rij
          drijdyi=-deltayj/rij
          drijdzi=-deltazj/rij
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
          drijdxk=0.0d0
          drijdyk=0.0d0
          drijdzk=0.0d0
!!
          fcutij=0.5d0*(dcos(pi*rij/funccutoff(i2,iindex))+1.d0)
          temp1=0.5d0*(-dsin(pi*rij/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
          dfcutijdxi=temp1*drijdxi
          dfcutijdyi=temp1*drijdyi
          dfcutijdzi=temp1*drijdzi
          dfcutijdxj=-1.d0*dfcutijdxi
          dfcutijdyj=-1.d0*dfcutijdyi
          dfcutijdzj=-1.d0*dfcutijdzi
          dfcutijdxk=0.0d0
          dfcutijdyk=0.0d0
          dfcutijdzk=0.0d0
!!
          do i4=lsta(1,i1),lsta(2,i1) ! loop over all neighbors of this atom
            if(symelement(i2,2,iindex).eq.lste(i4))then
!! new double counting criterion JB and TM 08.03.2010:
            if(((i4.gt.i3).and.(lste(i3).eq.lste(i4)))&
             .or.(lste(i3).ne.lste(i4)))then ! avoid double counting
!! This was a bug, because the order of neighbor atoms is not arbitrary
!!            if(i4.gt.i3)then ! avoid double counting
!!
            m=lstc(i4)
            rik=lstb(i4,4)
!!
            if(rik.le.funccutoff(i2,iindex))then
!!
              deltaxk=-1.d0*(xyzstruct(1,i1)-lstb(i4,1))
              deltayk=-1.d0*(xyzstruct(2,i1)-lstb(i4,2))
              deltazk=-1.d0*(xyzstruct(3,i1)-lstb(i4,3))
              drikdxi=-deltaxk/rik
              drikdyi=-deltayk/rik
              drikdzi=-deltazk/rik
              drikdxk=-1.d0*drikdxi
              drikdyk=-1.d0*drikdyi
              drikdzk=-1.d0*drikdzi
              drikdxj=0.0d0
              drikdyj=0.0d0
              drikdzj=0.0d0
!!
              fcutik=0.5d0*(dcos(pi*rik/funccutoff(i2,iindex))+1.d0)
              temp1=0.5d0*(-dsin(pi*rik/funccutoff(i2,iindex)))*(pi/funccutoff(i2,iindex))
              dfcutikdxi=temp1*drikdxi
              dfcutikdyi=temp1*drikdyi
              dfcutikdzi=temp1*drikdzi
              dfcutikdxj=0.0d0
              dfcutikdyj=0.0d0
              dfcutikdzj=0.0d0
              dfcutikdxk=-1.d0*dfcutikdxi
              dfcutikdyk=-1.d0*dfcutikdyi
              dfcutikdzk=-1.d0*dfcutikdzi
!!
!! calculate the third bond length between atom i3 and i4
!!
              rjk=(lstb(i3,1)-lstb(i4,1))**2 +&
                  (lstb(i3,2)-lstb(i4,2))**2 +&
                  (lstb(i3,3)-lstb(i4,3))**2
              rjk=dsqrt(rjk)
              if(rjk.le.rmin) then
                lrmin=.false.
!!                write(ounit,*)'Error rjk .le. rmin ',i2,rjk
!!                stop
              endif
!!
                drjkdxj=(lstb(i3,1)-lstb(i4,1))/rjk
                drjkdyj=(lstb(i3,2)-lstb(i4,2))/rjk
                drjkdzj=(lstb(i3,3)-lstb(i4,3))/rjk
                drjkdxk=-1.d0*drjkdxj
                drjkdyk=-1.d0*drjkdyj
                drjkdzk=-1.d0*drjkdzj
                drjkdxi=0.0d0
                drjkdyi=0.0d0
                drjkdzi=0.0d0
!!
!! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
                f=lambda(i2,iindex)*(rjk**2 - rij**2 -rik**2)
                g=-2.d0*rij*rik
                costheta=f/g
                costheta=1.d0+costheta ! avoid negative values
                temp4   =2.d0**(1.d0-zeta(i2,iindex))
                costheta=temp4*(costheta**zeta(i2,iindex))
!!
!! calculate the derivatives of costheta already here
!! (f/g)' = (f'g - fg')/g^2
                if(ldoforces) then
                 dfdxi=lambda(i2,iindex)*(-2.d0*rij*drijdxi - 2.d0*rik*drikdxi)
                 dfdyi=lambda(i2,iindex)*(-2.d0*rij*drijdyi - 2.d0*rik*drikdyi)
                 dfdzi=lambda(i2,iindex)*(-2.d0*rij*drijdzi - 2.d0*rik*drikdzi)
!!
                 dfdxj=lambda(i2,iindex)*(2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj)
                 dfdyj=lambda(i2,iindex)*(2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj)
                 dfdzj=lambda(i2,iindex)*(2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj)
!!
                 dfdxk=lambda(i2,iindex)*(2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk)
                 dfdyk=lambda(i2,iindex)*(2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk)
                 dfdzk=lambda(i2,iindex)*(2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk)
!!
                 dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                 dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                 dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)
!!
                 dgdxj=-2.d0*drijdxj*rik
                 dgdyj=-2.d0*drijdyj*rik
                 dgdzj=-2.d0*drijdzj*rik
!!
                 dgdxk=-2.d0*rij*drikdxk
                 dgdyk=-2.d0*rij*drikdyk
                 dgdzk=-2.d0*rij*drikdzk
!!
                 temp4=temp4*zeta(i2,iindex)*(1.d0+f/g)**(zeta(i2,iindex)-1)/g**2
                 dcosthetadxi=temp4*(dfdxi*g - f*dgdxi) !/g**2
                 dcosthetadyi=temp4*(dfdyi*g - f*dgdyi) !/g**2
                 dcosthetadzi=temp4*(dfdzi*g - f*dgdzi) !/g**2
                 dcosthetadxj=temp4*(dfdxj*g - f*dgdxj) !/g**2
                 dcosthetadyj=temp4*(dfdyj*g - f*dgdyj) !/g**2
                 dcosthetadzj=temp4*(dfdzj*g - f*dgdzj) !/g**2
                 dcosthetadxk=temp4*(dfdxk*g - f*dgdxk) !/g**2
                 dcosthetadyk=temp4*(dfdyk*g - f*dgdyk) !/g**2
                 dcosthetadzk=temp4*(dfdzk*g - f*dgdzk) !/g**2
                endif ! ldoforces for derivative of costheta
!!
!! calculation of exponential term
                expxyz=dexp(-eta(i2,iindex)*(rij**2+rik**2))
                temp1=-eta(i2,iindex)*2.0d0*expxyz
                dexpxyzdxi=(rij*drijdxi+rik*drikdxi)*temp1
                dexpxyzdyi=(rij*drijdyi+rik*drikdyi)*temp1
                dexpxyzdzi=(rij*drijdzi+rik*drikdzi)*temp1
                dexpxyzdxj=(rij*drijdxj+rik*drikdxj)*temp1
                dexpxyzdyj=(rij*drijdyj+rik*drikdyj)*temp1
                dexpxyzdzj=(rij*drijdzj+rik*drikdzj)*temp1
                dexpxyzdxk=(rij*drijdxk+rik*drikdxk)*temp1
                dexpxyzdyk=(rij*drijdyk+rik*drikdyk)*temp1
                dexpxyzdzk=(rij*drijdzk+rik*drikdzk)*temp1
!!
!! calculation of the symmetry function
          symfunction(i2,i1)=symfunction(i2,i1)&
                                 +costheta*expxyz*fcutij*fcutik
!!
              if(ldoforces)then
!! Calculation of derivatives for forces 
!! (fgh)' = f'gh + fg'h + fgh'
!! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'
!!
!! dsfunc/dx
               temp1=(+dcosthetadxi* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdxi* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdxi* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdxi)
               temp2=(+dcosthetadxj* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdxj* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdxj* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdxj)
               temp3=(+dcosthetadxk* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdxk* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdxk* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdxk)
            dsfuncdxyz(i2,i1,i1,1)=dsfuncdxyz(i2,i1,i1,1)+&
                                       temp1 
            dsfuncdxyz(i2,i1,n,1) =dsfuncdxyz(i2,i1,n,1)+&
                                       temp2 
            dsfuncdxyz(i2,i1,m,1) =dsfuncdxyz(i2,i1,m,1)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,1,i2,i1)=strs(1,1,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,1,i2,i1)=strs(2,1,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,1,i2,i1)=strs(3,1,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
!!
!! dsfunc/dy
               temp1=(+dcosthetadyi* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdyi* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdyi* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdyi)
               temp2=(+dcosthetadyj* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdyj* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdyj* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdyj)
               temp3=(+dcosthetadyk* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdyk* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdyk* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdyk)
            dsfuncdxyz(i2,i1,i1,2)=dsfuncdxyz(i2,i1,i1,2)+&
                                       temp1 
            dsfuncdxyz(i2,i1,n,2) =dsfuncdxyz(i2,i1,n,2)+&
                                       temp2 
            dsfuncdxyz(i2,i1,m,2) =dsfuncdxyz(i2,i1,m,2)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,2,i2,i1)=strs(1,2,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,2,i2,i1)=strs(2,2,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,2,i2,i1)=strs(3,2,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
!!
!! dsfunc/dz
               temp1=(+dcosthetadzi* expxyz   * fcutij   * fcutik   &
                      + costheta   *dexpxyzdzi* fcutij   * fcutik   &
                      + costheta   * expxyz   *dfcutijdzi* fcutik   &
                      + costheta   * expxyz   * fcutij   *dfcutikdzi)
               temp2=(+dcosthetadzj* expxyz   * fcutij   * fcutik   &
                      + costheta   *dexpxyzdzj* fcutij   * fcutik   &
                      + costheta   * expxyz   *dfcutijdzj* fcutik   &
                      + costheta   * expxyz   * fcutij   *dfcutikdzj)
               temp3=(+dcosthetadzk* expxyz   * fcutij   * fcutik   &
                      + costheta   *dexpxyzdzk* fcutij   * fcutik   &
                      + costheta   * expxyz   *dfcutijdzk* fcutik   &
                      + costheta   * expxyz   * fcutij   *dfcutikdzk)
            dsfuncdxyz(i2,i1,i1,3)=dsfuncdxyz(i2,i1,i1,3)+&
                                       temp1 
            dsfuncdxyz(i2,i1,n,3) =dsfuncdxyz(i2,i1,n,3)+&
                                       temp2 
            dsfuncdxyz(i2,i1,m,3) =dsfuncdxyz(i2,i1,m,3)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs(1,3,i2,i1)=strs(1,3,i2,i1)+deltaxj*temp2&
                  +deltaxk*temp3
                strs(2,3,i2,i1)=strs(2,3,i2,i1)+deltayj*temp2&
                  +deltayk*temp3
                strs(3,3,i2,i1)=strs(3,3,i2,i1)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
              endif ! ldoforces
!!
            endif ! rik .le. cutoff
!!
          endif ! i4.gt.i3
          endif ! (symelement(i2,2,iindex).eq.lste(i4))then
          enddo ! i4
          endif ! rij.le.funccutoff(i2,iindex)
          endif ! (symelement(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
       write(*,*)'Error: function_type not implemented ',function_type(i2,iindex)
      endif
      enddo ! i2
      enddo ! i1
!!
      return
      end
