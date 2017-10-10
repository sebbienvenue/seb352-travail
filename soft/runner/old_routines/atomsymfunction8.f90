!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconefunction_para.f90
!!
      subroutine atomsymfunction8(i1,i2,iindex,natoms,natomsdim,nelem,&
        max_num_atoms,max_num_neighbors,invneighboridx,&
        jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
        cutoff_type,lstb,funccutoff_local,pi,xyzstruct,eta_local,rshift_local,&
        lambda_local,rmin,symfunction,dsfuncdxyz,strs,&
        ldoforces,ldostress,lrmin)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2,i3,i4
      integer listdim
      integer nelem 
      integer natomsdim 
      integer natoms 
      integer iindex
      integer max_num_atoms 
      integer max_num_neighbors 
      integer maxnum_funcvalues_local 
      integer symelement_local(maxnum_funcvalues_local,2,nelem)          ! in
      integer lsta(2,max_num_atoms)
      integer lstc(listdim)
      integer lste(listdim)
      integer n
      integer m
      integer jcount 
      integer cutoff_type
      integer invneighboridx(natoms,max_num_atoms) 
!!
      real*8 strs(3,3,maxnum_funcvalues_local,natomsdim)
      real*8 lstb(listdim,4)
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)
      real*8 lambda_local(maxnum_funcvalues_local,nelem)                ! in
      real*8 rij
      real*8 rik
      real*8 rjk
      real*8 symfunction(maxnum_funcvalues_local,natomsdim)
      real*8 xyzstruct(3,max_num_atoms)
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
      real*8 temp1
      real*8 temp2
      real*8 temp3
      real*8 temp4
      real*8 pi
      real*8 dsfuncdxyz(maxnum_funcvalues_local,natomsdim,0:max_num_neighbors,3) 
      real*8 eta_local(maxnum_funcvalues_local,nelem)                            ! in
      real*8 rshift_local(maxnum_funcvalues_local,nelem)                ! in
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
      real*8 theta
      real*8 dthetadxi,dthetadyi,dthetadzi
      real*8 dthetadxj,dthetadyj,dthetadzj
      real*8 dthetadxk,dthetadyk,dthetadzk
      real*8 rad2deg
      real*8 dgdxi,dgdyi,dgdzi
      real*8 dgdxj,dgdyj,dgdzj
      real*8 dgdxk,dgdyk,dgdzk
      real*8 dfdxi,dfdyi,dfdzi
      real*8 dfdxj,dfdyj,dfdzj
      real*8 dfdxk,dfdyk,dfdzk
      real*8 rmin
      real*8 f 
      real*8 g 
      real*8 dexptempdxi
      real*8 dexptempdyi
      real*8 dexptempdzi
      real*8 dexptempdxj
      real*8 dexptempdyj
      real*8 dexptempdzj
      real*8 dexptempdxk
      real*8 dexptempdyk
      real*8 dexptempdzk
      real*8 exptemp
!!
      logical ldoforces
      logical ldostress
      logical lrmin
!!
      rad2deg            =180.d0/pi
!!
        do i3=lsta(1,i1),lsta(2,i1) ! loop over neighbors
          if(symelement_local(i2,1,iindex).eq.lste(i3))then
          rij=lstb(i3,4)
          n=lstc(i3)
          if(rij.le.funccutoff_local(i2,iindex))then
!!
          deltaxj=-1.d0*(xyzstruct(1,jcount)-lstb(i3,1))
          deltayj=-1.d0*(xyzstruct(2,jcount)-lstb(i3,2))
          deltazj=-1.d0*(xyzstruct(3,jcount)-lstb(i3,3))
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
          if(cutoff_type.eq.1)then
            fcutij=0.5d0*(dcos(pi*rij/funccutoff_local(i2,iindex))+1.d0)
            temp1=0.5d0*(-dsin(pi*rij/funccutoff_local(i2,iindex)))*(pi/funccutoff_local(i2,iindex))
          elseif(cutoff_type.eq.2)then
            fcutij=(tanh(1.d0-rij/funccutoff_local(i2,iindex)))**3
            temp1 =(-3.d0/funccutoff_local(i2,iindex))*&
               ((tanh(1.d0-rij/funccutoff_local(i2,iindex)))**2 &
              - (tanh(1.d0-rij/funccutoff_local(i2,iindex)))**4)
          else
            write(ounit,*)'ERROR: Unknown cutoff_type in calconefunction_para'
            stop !'
          endif
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
            if(symelement_local(i2,2,iindex).eq.lste(i4))then
!! new double counting criterion JB and TM 08.03.2010:
            if(((i4.gt.i3).and.(lste(i3).eq.lste(i4)))&
             .or.(lste(i3).ne.lste(i4)))then ! avoid double counting
!!
            m=lstc(i4)
            rik=lstb(i4,4)
!!
            if(rik.le.funccutoff_local(i2,iindex))then
!!
              deltaxk=-1.d0*(xyzstruct(1,jcount)-lstb(i4,1))
              deltayk=-1.d0*(xyzstruct(2,jcount)-lstb(i4,2))
              deltazk=-1.d0*(xyzstruct(3,jcount)-lstb(i4,3))
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
              if(cutoff_type.eq.1)then
                fcutik=0.5d0*(dcos(pi*rik/funccutoff_local(i2,iindex))+1.d0)
                temp1=0.5d0*(-dsin(pi*rik/funccutoff_local(i2,iindex)))*(pi/funccutoff_local(i2,iindex))
              elseif(cutoff_type.eq.2)then
                fcutik=(tanh(1.d0-rik/funccutoff_local(i2,iindex)))**3
                temp1 =(-3.d0/funccutoff_local(i2,iindex))*&
                   ((tanh(1.d0-rik/funccutoff_local(i2,iindex)))**2 &
                  - (tanh(1.d0-rik/funccutoff_local(i2,iindex)))**4)
              else
                write(ounit,*)'ERROR: Unknown cutoff_type in calconefunction_para'
                stop !'
              endif
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
              if(rjk.le.funccutoff_local(i2,iindex)) then
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
                if(cutoff_type.eq.1)then
                  fcutjk=0.5d0*(dcos(pi*rjk/funccutoff_local(i2,iindex))+1.d0)
                  temp1=0.5d0*(-dsin(pi*rjk/funccutoff_local(i2,iindex)))*(pi/funccutoff_local(i2,iindex))
                elseif(cutoff_type.eq.2)then
                  fcutjk=(tanh(1.d0-rjk/funccutoff_local(i2,iindex)))**3
                  temp1 =(-3.d0/funccutoff_local(i2,iindex))*&
                     ((tanh(1.d0-rjk/funccutoff_local(i2,iindex)))**2 &
                    - (tanh(1.d0-rjk/funccutoff_local(i2,iindex)))**4)
                else
                  write(ounit,*)'ERROR: Unknown cutoff_type in calconefunction_para'
                  stop
                endif
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
                 dfdxi=lambda_local(i2,iindex)*(-2.d0*rij*drijdxi - 2.d0*rik*drikdxi)
                 dfdyi=lambda_local(i2,iindex)*(-2.d0*rij*drijdyi - 2.d0*rik*drikdyi)
                 dfdzi=lambda_local(i2,iindex)*(-2.d0*rij*drijdzi - 2.d0*rik*drikdzi)
!!
                 dfdxj=lambda_local(i2,iindex)*(2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj)
                 dfdyj=lambda_local(i2,iindex)*(2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj)
                 dfdzj=lambda_local(i2,iindex)*(2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj)
!!
                 dfdxk=lambda_local(i2,iindex)*(2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk)
                 dfdyk=lambda_local(i2,iindex)*(2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk)
                 dfdzk=lambda_local(i2,iindex)*(2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk)
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
                temp1=theta-rshift_local(i2,iindex)
                temp2=theta+rshift_local(i2,iindex)
                exptemp=(dexp(-eta_local(i2,iindex)*(temp1       )**2)&
                        +dexp(-eta_local(i2,iindex)*(temp2-360.d0)**2)&
                        +dexp(-eta_local(i2,iindex)*(temp2       )**2)&
                        +dexp(-eta_local(i2,iindex)*(temp1-360.d0)**2))
          symfunction(i2,i1)=symfunction(i2,i1)&
                  +exptemp*fcutij*fcutik*fcutjk
!!
              if(ldoforces)then
!! Calculation of derivatives for forces
!! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'
                temp3=-2.d0*eta_local(i2,iindex)*temp1*dexp(-eta_local(i2,iindex)*(temp1       )**2)&
                      -2.d0*eta_local(i2,iindex)*temp2*dexp(-eta_local(i2,iindex)*(temp2-360.d0)**2)&
                      -2.d0*eta_local(i2,iindex)*temp2*dexp(-eta_local(i2,iindex)*(temp2       )**2)&
                      -2.d0*eta_local(i2,iindex)*temp1*dexp(-eta_local(i2,iindex)*(temp1-360.d0)**2)
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
            dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),1)=dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),1)+&
                                       temp1
            dsfuncdxyz(i2,i1,invneighboridx(i1,n),1) =dsfuncdxyz(i2,i1,invneighboridx(i1,n),1)+&
                                       temp2
            dsfuncdxyz(i2,i1,invneighboridx(i1,m),1) =dsfuncdxyz(i2,i1,invneighboridx(i1,m),1)+&
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
            dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),2)=dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),2)+&
                                       temp1
            dsfuncdxyz(i2,i1,invneighboridx(i1,n),2) =dsfuncdxyz(i2,i1,invneighboridx(i1,n),2)+&
                                       temp2
            dsfuncdxyz(i2,i1,invneighboridx(i1,m),2) =dsfuncdxyz(i2,i1,invneighboridx(i1,m),2)+&
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
            dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),3)=dsfuncdxyz(i2,i1,invneighboridx(i1,jcount),3)+&
                                       temp1
            dsfuncdxyz(i2,i1,invneighboridx(i1,n),3) =dsfuncdxyz(i2,i1,invneighboridx(i1,n),3)+&
                                       temp2
            dsfuncdxyz(i2,i1,invneighboridx(i1,m),3) =dsfuncdxyz(i2,i1,invneighboridx(i1,m),3)+&
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
          endif !(symelement_local(i2,2,iindex).eq.lste(i4))then
          enddo ! i4
          endif ! rij.le.funccutoff_local(i2,iindex)
          endif ! (symelement_local(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!
      return
      end
