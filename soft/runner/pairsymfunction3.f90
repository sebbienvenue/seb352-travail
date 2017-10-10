!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconepairfunction.f90
!!
      subroutine pairsymfunction3(i1,i2,jcount,num_atoms,listdim,iindex,&
        pair_lstc,maxnum_funcvaluesp,npairs,max_num_pairs,&
        pair_lste,sympelement,&
        max_num_atoms,pairs_atom,rshiftp,etap,&
        pair_lsta,pair_lstb,pairs_xyz,funccutoffp,strs_pair,&
        dsfuncdxyz_pair,symfunctionp,ldoforces,ldostress)
!!
      use fileunits
      use nnconstants
!!
      implicit none
!!
      integer l,m,n
      integer i1                                                               ! in, pair
      integer i2                                                               ! in, symfunction
      integer i3                                                               ! internal
      integer jcount                                                           ! in
      integer listdim                                                          ! in
      integer num_atoms                                                        ! in
      integer pairs_atom(listdim,2)                                            ! in
      integer pair_lsta(listdim,2)                                             ! in
      integer pair_lstc(listdim)                                               ! in
      integer atomsinpairs(num_atoms)                                          ! internal
      integer maxnum_funcvaluesp                                               ! in
      integer max_num_pairs                                                    ! in
      integer max_num_atoms                                                    ! in
      integer npairs                                                           ! in
      integer iindex                                                           ! in
      integer pair_lste(listdim)                                               ! in
      integer sympelement(maxnum_funcvaluesp,2,npairs)                         ! in
!!
      real*8 rshiftp(maxnum_funcvaluesp,npairs)                                ! in
      real*8 etap(maxnum_funcvaluesp,npairs)                                   ! in
      real*8 dexpijdxi,dexpijdyi,dexpijdzi
      real*8 dexpijdxj,dexpijdyj,dexpijdzj
      real*8 dexpijdxk,dexpijdyk,dexpijdzk
      real*8 dexpjkdxi,dexpjkdyi,dexpjkdzi
      real*8 dexpjkdxj,dexpjkdyj,dexpjkdzj
      real*8 dexpjkdxk,dexpjkdyk,dexpjkdzk
      real*8 dexpikdxi,dexpikdyi,dexpikdzi
      real*8 dexpikdxj,dexpikdyj,dexpikdzj
      real*8 dexpikdxk,dexpikdyk,dexpikdzk
      real*8 expij,expjk,expik
      real*8 fcutij
      real*8 fcutik
      real*8 fcutjk
      real*8 pair_lstb(listdim,5)                                              ! in
      real*8 pairs_xyz(listdim,2,3)                                            ! in
      real*8 deltaxij,deltayij,deltazij                                        ! internal
      real*8 deltaxik,deltayik,deltazik
      real*8 rij                                                               ! internal
      real*8 rjk                                                               ! internal
      real*8 rik                                                               ! internal
      real*8 drijdxi, drijdyi, drijdzi                                         ! internal
      real*8 drijdxj, drijdyj, drijdzj                                         ! internal
      real*8 drijdxk, drijdyk, drijdzk                                         ! internal
      real*8 drikdxi, drikdyi, drikdzi
      real*8 drikdxj, drikdyj, drikdzj
      real*8 drikdxk, drikdyk, drikdzk
      real*8 drjkdxi, drjkdyi, drjkdzi
      real*8 drjkdxj, drjkdyj, drjkdzj
      real*8 drjkdxk, drjkdyk, drjkdzk
      real*8 dfcutijdxi,dfcutijdyi,dfcutijdzi
      real*8 dfcutijdxj,dfcutijdyj,dfcutijdzj
      real*8 dfcutijdxk,dfcutijdyk,dfcutijdzk
      real*8 dfcutikdxi,dfcutikdyi,dfcutikdzi
      real*8 dfcutikdxj,dfcutikdyj,dfcutikdzj
      real*8 dfcutikdxk,dfcutikdyk,dfcutikdzk
      real*8 dfcutjkdxi,dfcutjkdyi,dfcutjkdzi
      real*8 dfcutjkdxj,dfcutjkdyj,dfcutjkdzj
      real*8 dfcutjkdxk,dfcutjkdyk,dfcutjkdzk
      real*8 funccutoffp(maxnum_funcvaluesp,npairs)                            ! in
      real*8 temp1,temp2,temp3,temp4                                           ! internal
      real*8 strs_pair(3,3,maxnum_funcvaluesp,max_num_pairs)                   ! out
      real*8 dsfuncdxyz_pair(maxnum_funcvaluesp,max_num_pairs,max_num_atoms,3) ! out
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs)                    ! out
!!
      logical ldoforces
      logical ldostress
!!
      do i3 = pair_lsta(i1,1), pair_lsta(i1,2)
!! check if nuclear charge of neighbor i3 is correct for this symfunction
        if(sympelement(i2,1,iindex).eq.pair_lste(i3))then
!!
          l = pair_lstc(i3)    ! neighbor of pair
          n = pairs_atom(i1,1) ! atom A in pair
          m = pairs_atom(i1,2) ! atom B in pair
!!
          atomsinpairs(m)=atomsinpairs(m)+1
          atomsinpairs(n)=atomsinpairs(n)+1
!!
!! calculate distance between atom A=j in pair and neighbor i
          deltaxij= (pair_lstb(i3,1) - pairs_xyz(jcount,1,1))
          deltayij= (pair_lstb(i3,2) - pairs_xyz(jcount,1,2))
          deltazij= (pair_lstb(i3,3) - pairs_xyz(jcount,1,3))
          rij = deltaxij**2 + deltayij**2 + deltazij**2
          rij = sqrt(rij)
!!
          fcutij    =0.0d0       
          expij=0.0d0
          dfcutijdxi=0.0d0  
          dfcutijdyi=0.0d0  
          dfcutijdzi=0.0d0
          dfcutijdxj=0.0d0  
          dfcutijdyj=0.0d0  
          dfcutijdzj=0.0d0
          dfcutijdxk=0.0d0  
          dfcutijdyk=0.0d0  
          dfcutijdzk=0.0d0

          dexpijdxj =0.0d0  
          dexpijdyj =0.0d0  
          dexpijdzj =0.0d0
          dexpijdxi =0.0d0  
          dexpijdyi =0.0d0  
          dexpijdzi =0.0d0
          dexpijdxk =0.0d0  
          dexpijdyk =0.0d0  
          dexpijdzk =0.0d0
!!
!! if neighbor i is within cutoff
          if(rij.le.funccutoffp(i2,iindex))then
!!
            drijdxi=deltaxij/rij
            drijdyi=deltayij/rij
            drijdzi=deltazij/rij
            drijdxj=-1.d0*drijdxi
            drijdyj=-1.d0*drijdyi
            drijdzj=-1.d0*drijdzi
            drijdxk= 0.0d0
            drijdyk= 0.0d0
            drijdzk= 0.0d0
!! 
            fcutij=(0.5d0*(dcos(pi*rij/funccutoffp(i2,iindex))+1.d0))
            temp1=0.5d0*(-dsin(pi*rij/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutijdxi=temp1*drijdxi
            dfcutijdyi=temp1*drijdyi
            dfcutijdzi=temp1*drijdzi
            dfcutijdxj=-1.d0*dfcutijdxi
            dfcutijdyj=-1.d0*dfcutijdyi
            dfcutijdzj=-1.d0*dfcutijdzi
            dfcutijdxk= 0.0d0
            dfcutijdyk= 0.0d0
            dfcutijdzk= 0.0d0
            expij=(dexp(-1.d0*etap(i2,iindex)* (rij-rshiftp(i2,iindex))**2)) 
!!
            dexpijdxj=(2.0d0*expij*deltaxij* (rij-rshiftp(i2,iindex)) *etap(i2,iindex))/rij
            dexpijdyj=(2.0d0*expij*deltayij* (rij-rshiftp(i2,iindex)) *etap(i2,iindex))/rij
            dexpijdzj=(2.0d0*expij*deltazij* (rij-rshiftp(i2,iindex)) *etap(i2,iindex))/rij
            dexpijdxi=-1.0d0*dexpijdxj
            dexpijdyi=-1.0d0*dexpijdyj
            dexpijdzi=-1.0d0*dexpijdzj
            dexpijdxk= 0.0d0
            dexpijdyk= 0.0d0
            dexpijdzk= 0.0d0
          endif ! rij
!!
!! calculate distance between atom B=k in pair and neighbor i
          deltaxik= (pair_lstb(i3,1) - pairs_xyz(jcount,2,1))
          deltayik= (pair_lstb(i3,2) - pairs_xyz(jcount,2,2))
          deltazik= (pair_lstb(i3,3) - pairs_xyz(jcount,2,3))
          rik= deltaxik**2 + deltayik**2 + deltazik**2
          rik= sqrt(rik)
!!
          fcutik    =0.0d0  
          expik     =0.0d0
          dfcutikdxi=0.0d0  
          dfcutikdyi=0.0d0  
          dfcutikdzi=0.0d0
          dfcutikdxj=0.0d0  
          dfcutikdyj=0.0d0  
          dfcutikdzj=0.0d0
          dfcutikdxk=0.0d0  
          dfcutikdyk=0.0d0  
          dfcutikdzk=0.0d0

          dexpikdxk =0.0d0  
          dexpikdyk =0.0d0  
          dexpikdzk =0.0d0
          dexpikdxi =0.0d0  
          dexpikdyi =0.0d0  
          dexpikdzi =0.0d0
          dexpikdxj =0.0d0  
          dexpikdyj =0.0d0  
          dexpikdzj =0.0d0
!!
!! check if neighbor i is within cutoff of atom B=k 
          if(rik.le.funccutoffp(i2,iindex))then
            drikdxi= deltaxik/rik
            drikdyi= deltayik/rik
            drikdzi= deltazik/rik
            drikdxk=-1.d0*drikdxi
            drikdyk=-1.d0*drikdyi
            drikdzk=-1.d0*drikdzi
            drikdxj= 0.0d0
            drikdyj= 0.0d0
            drikdzj= 0.0d0
!! 
            fcutik=(0.5d0*(dcos(pi*rik/funccutoffp(i2,iindex))+1.d0))
            temp1=0.5d0*(-dsin(pi*rik/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutikdxi= temp1*drikdxi
            dfcutikdyi= temp1*drikdyi
            dfcutikdzi= temp1*drikdzi
            dfcutikdxj= 0.0d0
            dfcutikdyj= 0.0d0
            dfcutikdzj= 0.0d0
            dfcutikdxk=-1.d0*dfcutikdxi
            dfcutikdyk=-1.d0*dfcutikdyi
            dfcutikdzk=-1.d0*dfcutikdzi
!!
            expik= (dexp(-1.d0*etap(i2,iindex)* (rik-rshiftp(i2,iindex))**2)) 
!!
            dexpikdxk= (2.0d0*expik*deltaxik* (rik-rshiftp(i2,iindex)) *etap(i2,iindex))/rik
            dexpikdyk= (2.0d0*expik*deltayik* (rik-rshiftp(i2,iindex)) *etap(i2,iindex))/rik
            dexpikdzk= (2.0d0*expik*deltazik* (rik-rshiftp(i2,iindex)) *etap(i2,iindex))/rik
            dexpikdxi= -1.0d0*dexpikdxk
            dexpikdyi= -1.0d0*dexpikdyk
            dexpikdzi= -1.0d0*dexpikdzk
            dexpikdxj=  0.0d0
            dexpikdyj=  0.0d0
            dexpikdzj=  0.0d0
          endif ! rik
!!
!! check distance between atoms A=j and B=k in pair
          rjk=(pairs_xyz(jcount,1,1)-pairs_xyz(jcount,2,1))**2 +&
              (pairs_xyz(jcount,1,2)-pairs_xyz(jcount,2,2))**2+ &
              (pairs_xyz(jcount,1,3)-pairs_xyz(jcount,2,3))**2
          rjk = sqrt(rjk)
!!
          fcutjk    =0.0d0  
          expjk     =0.0d0
          dfcutjkdxj=0.0d0  
          dfcutjkdyj=0.0d0 
          dfcutjkdzj=0.0d0
          dfcutjkdxk=0.0d0  
          dfcutjkdyk=0.0d0 
          dfcutjkdzk=0.0d0
          dfcutjkdxi=0.0d0  
          dfcutjkdyi=0.0d0 
          dfcutjkdzi=0.0d0

          dexpjkdxj =0.0d0 
          dexpjkdyj =0.0d0 
          dexpjkdzj =0.0d0
          dexpjkdxk =0.0d0 
          dexpjkdyk =0.0d0 
          dexpjkdzk =0.0d0
          dexpjkdxi =0.0d0 
          dexpjkdyi =0.0d0 
          dexpjkdzi =0.0d0
!!
!! check if distance between both atoms in pair is smaller than cutoff
          if(rjk.le.funccutoffp(i2,iindex))then
            drjkdxj=(pairs_xyz(jcount,1,1)-pairs_xyz(jcount,2,1))/rjk
            drjkdyj=(pairs_xyz(jcount,1,2)-pairs_xyz(jcount,2,2))/rjk
            drjkdzj=(pairs_xyz(jcount,1,3)-pairs_xyz(jcount,2,3))/rjk
            drjkdxk=-1.d0*drjkdxj
            drjkdyk=-1.d0*drjkdyj
            drjkdzk=-1.d0*drjkdzj
            drjkdxi= 0.0d0
            drjkdyi= 0.0d0
            drjkdzi= 0.0d0
!!
            fcutjk= (0.5d0*(dcos(pi*rjk/funccutoffp(i2,iindex))+1.d0))
            temp1=0.5d0*(-dsin(pi*rjk/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutjkdxj= temp1*drjkdxj
            dfcutjkdyj= temp1*drjkdyj
            dfcutjkdzj= temp1*drjkdzj
            dfcutjkdxk=-1.d0*dfcutjkdxj
            dfcutjkdyk=-1.d0*dfcutjkdyj
            dfcutjkdzk=-1.d0*dfcutjkdzj
            dfcutjkdxi= 0.0d0
            dfcutjkdyi= 0.0d0
            dfcutjkdzi= 0.0d0
!!
            expjk=(dexp(-1.d0*etap(i2,iindex)* (rjk)**2))
!!
            dexpjkdxj=(-2.0d0*expjk*(pairs_xyz(jcount,1,1)-pairs_xyz(jcount,2,1))*etap(i2,iindex)) 
            dexpjkdyj=(-2.0d0*expjk*(pairs_xyz(jcount,1,2)-pairs_xyz(jcount,2,2))*etap(i2,iindex))
            dexpjkdzj=(-2.0d0*expjk*(pairs_xyz(jcount,1,3)-pairs_xyz(jcount,2,3))*etap(i2,iindex))
            dexpjkdxk= -1.0d0*dexpjkdxj
            dexpjkdyk= -1.0d0*dexpjkdyj
            dexpjkdzk= -1.0d0*dexpjkdzj
            dexpjkdxi=  0.0d0
            dexpjkdyi=  0.0d0
            dexpjkdzi=  0.0d0
          endif ! rjk
!! 
          symfunctionp(i2,i1)=symfunctionp(i2,i1) + (((fcutij*expij) + (fcutik*expik)) * (fcutjk*expjk))
!!
!! dsfunc/dx
          if(ldoforces) then
            temp2 = (dfcutjkdxi   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdxi  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdxi  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdxi) +& 
                    (dfcutjkdxi   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdxi  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdxi  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdxi) 
!!
            temp3 = (dfcutjkdxj   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdxj  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdxj  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdxj) +&
                    (dfcutjkdxj   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdxj  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdxj  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdxj) 
!!
            temp4 = (dfcutjkdxk   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdxk  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdxk  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdxk) +&
                    (dfcutjkdxk   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdxk  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdxk  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdxk)
!! 
            dsfuncdxyz_pair(i2,i1,l,1) =dsfuncdxyz_pair(i2,i1,l,1) + temp2
            dsfuncdxyz_pair(i2,jcount,n,1) =dsfuncdxyz_pair(i2,jcount,n,1) + temp3
            dsfuncdxyz_pair(i2,jcount,m,1) =dsfuncdxyz_pair(i2,jcount,m,1) + temp4
          endif
!! 
          if(ldostress)then
            strs_pair(1,1,i2,i1)=strs_pair(1,1,i2,i1)+deltaxij*temp3 +deltaxik*temp4
            strs_pair(2,1,i2,i1)=strs_pair(2,1,i2,i1)+deltayij*temp3 +deltayik*temp4
            strs_pair(3,1,i2,i1)=strs_pair(3,1,i2,i1)+deltazij*temp3 +deltazik*temp4
          endif ! ldostress
!!
!! dsfunc/dy
          if(ldoforces) then
            temp2 = (dfcutjkdyi   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdyi  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdyi  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdyi) +&
                    (dfcutjkdyi   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdyi  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdyi  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdyi)
!!
            temp3 = (dfcutjkdyj   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdyj  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdyj  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdyj) +&
                    (dfcutjkdyj   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdyj  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdyj  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdyj)
!!  
            temp4 = (dfcutjkdyk   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdyk  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdyk  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdyk) +&
                    (dfcutjkdyk   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdyk  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdyk  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdyk)
!! 
            dsfuncdxyz_pair(i2,i1,l,2) =dsfuncdxyz_pair(i2,i1,l,2) + temp2
            dsfuncdxyz_pair(i2,jcount,n,2) =dsfuncdxyz_pair(i2,jcount,n,2) + temp3
            dsfuncdxyz_pair(i2,jcount,m,2) =dsfuncdxyz_pair(i2,jcount,m,2) + temp4
          endif
!!
          if(ldostress)then
            strs_pair(1,2,i2,i1)=strs_pair(1,2,i2,i1)+deltaxij*temp3 +deltaxik*temp4
            strs_pair(2,2,i2,i1)=strs_pair(2,2,i2,i1)+deltayij*temp3 +deltayik*temp4
            strs_pair(3,2,i2,i1)=strs_pair(3,2,i2,i1)+deltazij*temp3 +deltazik*temp4
          endif ! ldostress
!!
!! dsfunc/dz
          if(ldoforces) then
            temp2 = (dfcutjkdzi    *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdzi  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdzi  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdzi) +&
                    (dfcutjkdzi   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdzi  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdzi  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdzi)
!!
            temp3 = (dfcutjkdzj   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdzj  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdzj  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdzj) +&
                    (dfcutjkdzj   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdzj  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdzj  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdzj)
!!
            temp4 = (dfcutjkdzk   *  expjk      *  fcutij      *  expij) +&
                    (fcutjk       *  dexpjkdzk  *  fcutij      *  expij) +&
                    (fcutjk       *  expjk      *  dfcutijdzk  *  expij) +&
                    (fcutjk       *  expjk      *  fcutij      *  dexpijdzk) +&
                    (dfcutjkdzk   *  expjk      *  fcutik      *  expik) +&
                    (fcutjk       *  dexpjkdzk  *  fcutik      *  expik) +&
                    (fcutjk       *  expjk      *  dfcutikdzk  *  expik) +&
                    (fcutjk       *  expjk      *  fcutik      *  dexpikdzk)
!!
!!            write(ounit,*)dsfuncdxyz_pair(i2,i1,l,3),i1,i2,l,temp2
            dsfuncdxyz_pair(i2,i1,l,3) =dsfuncdxyz_pair(i2,i1,l,3) + temp2
            dsfuncdxyz_pair(i2,jcount,n,3) =dsfuncdxyz_pair(i2,jcount,n,3) + temp3
            dsfuncdxyz_pair(i2,jcount,m,3) =dsfuncdxyz_pair(i2,jcount,m,3) + temp4
          endif
!!
          if(ldostress)then
            strs_pair(1,3,i2,i1)=strs_pair(1,3,i2,i1)+deltaxij*temp3 +deltaxik*temp4
            strs_pair(2,3,i2,i1)=strs_pair(2,3,i2,i1)+deltayij*temp3 +deltayik*temp4
            strs_pair(3,3,i2,i1)=strs_pair(3,3,i2,i1)+deltazij*temp3 +deltazik*temp4
          endif ! ldostress
!!
        endif ! sympelement(i2,1,iindex) ==  pair_lste(i3)
      enddo  ! loop over neighbors
!!
      return
      end
