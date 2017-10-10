!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconepairfunction.f90
!!
      subroutine pairsymfunction4(i1,i2,jcount,num_atoms,listdim,iindex,&
        pair_lstc,maxnum_funcvaluesp,npairs,max_num_pairs,&
        pair_lste,sympelement,&
        max_num_atoms,pairs_atom,pairs_charge,etap,zetap,lambdap,&
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
      integer pair_lstc(listdim)                                               ! in
      integer atomsinpairs(num_atoms)                                          ! internal
      integer maxnum_funcvaluesp                                               ! in
      integer max_num_pairs                                                    ! in
      integer max_num_atoms                                                    ! in
      integer npairs                                                           ! in
      integer iindex                                                           ! in
      integer pairs_charge(2,listdim)                                          ! in
      integer pair_lsta(listdim,2)                                             ! in
      integer pair_lste(listdim)                                               ! in
      integer sympelement(maxnum_funcvaluesp,2,npairs)                         ! in
!!
      real*8 etap(maxnum_funcvaluesp,npairs)                                   ! in
      real*8 zetap(maxnum_funcvaluesp,npairs)                                  ! in
      real*8 lambdap(maxnum_funcvaluesp,npairs)                                ! in
      real*8 costheta
      real*8 dcosthetadxi,dcosthetadyi,dcosthetadzi
      real*8 dcosthetadxj,dcosthetadyj,dcosthetadzj
      real*8 dcosthetadxk,dcosthetadyk,dcosthetadzk
      real*8 dgdxi,dgdyi,dgdzi
      real*8 dgdxj,dgdyj,dgdzj
      real*8 dgdxk,dgdyk,dgdzk
      real*8 dfdxi,dfdyi,dfdzi
      real*8 dfdxj,dfdyj,dfdzj
      real*8 dfdxk,dfdyk,dfdzk
      real*8 f
      real*8 g
      real*8 dfcutij
      real*8 dfcutjk
      real*8 dfcutik
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
      real*8 funccutoffp(maxnum_funcvaluesp,npairs)          ! in
      real*8 temp1,temp2,temp3,temp4                                           ! internal
      real*8 strs_pair(3,3,maxnum_funcvaluesp,max_num_pairs)                   ! out
      real*8 dsfuncdxyz_pair(maxnum_funcvaluesp,max_num_pairs,max_num_atoms,3) ! out
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs)                     ! out
!!
      logical ldoforces
      logical ldostress
!!
!! loop over all neighbors of pair
      do i3 = pair_lsta(i1,1), pair_lsta(i1,2)
!! check if neighbor has the right nuclear charge for this symmetry function
        if(sympelement(i2,1,iindex).eq.pair_lste(i3))then
!!
          l = pair_lstc(i3)     ! neighbor of pair
          n = pairs_atom(i1,1)  ! atom A in pair        
          m = pairs_atom(i1,2)  ! atom B in pair        
!!
          atomsinpairs(m)=atomsinpairs(m)+1
          atomsinpairs(n)=atomsinpairs(n)+1
!!
!! calculate distance between atom A=j and neighbor i
          deltaxij= (pair_lstb(i3,1) - pairs_xyz(jcount,1,1))
          deltayij= (pair_lstb(i3,2) - pairs_xyz(jcount,1,2))
          deltazij= (pair_lstb(i3,3) - pairs_xyz(jcount,1,3))
          rij = deltaxij**2 + deltayij**2 + deltazij**2           
          rij = sqrt(rij)
!!
          fcutij    =0.0d0       
          expij     =0.0d0
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

          drijdxi=0.0d0 
          drijdyi=0.0d0 
          drijdzi=0.0d0
          drijdxj=0.0d0
          drijdyj=0.0d0 
          drijdzj=0.0d0
          drijdxk=0.0d0 
          drijdyk=0.0d0 
          drijdzk=0.0d0
!!        
!! check if distance between atom A=j in pair and neighbor i is smaller than cutoff       
          if(rij.le.funccutoffp(i2,iindex))then
            drijdxi=deltaxij/rij
            drijdyi=deltayij/rij
            drijdzi=deltazij/rij
            drijdxj=-1.d0*drijdxi
            drijdyj=-1.d0*drijdyi
            drijdzj=-1.d0*drijdzi
            drijdxk=0.0d0
            drijdyk=0.0d0
            drijdzk=0.0d0
!!                                       
            fcutij=(0.5d0*(dcos(pi*rij/funccutoffp(i2,iindex))+1.d0))
            dfcutij=0.5d0*(-dsin(pi*rij/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutijdxi=dfcutij*drijdxi
            dfcutijdyi=dfcutij*drijdyi
            dfcutijdzi=dfcutij*drijdzi
            dfcutijdxj=-1.d0*dfcutijdxi
            dfcutijdyj=-1.d0*dfcutijdyi
            dfcutijdzj=-1.d0*dfcutijdzi
            dfcutijdxk=0.0d0
            dfcutijdyk=0.0d0
            dfcutijdzk=0.0d0
!!
            expij=(dexp(-etap(i2,iindex)* (rij)**2))
!!
            dexpijdxi=-2.0d0*expij*deltaxij*etap(i2,iindex)
            dexpijdyi=-2.0d0*expij*deltayij*etap(i2,iindex)
            dexpijdzi=-2.0d0*expij*deltazij*etap(i2,iindex)
            dexpijdxj=-1.0d0*dexpijdxi
            dexpijdyj=-1.0d0*dexpijdyi
            dexpijdzj=-1.0d0*dexpijdzi
            dexpijdxk=0.0d0
            dexpijdyk=0.0d0
            dexpijdzk=0.0d0
          endif

!! calculate distance between atom B=k and neighbor i
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

          drikdxi=0.0d0 
          drikdyi=0.0d0 
          drikdzi=0.0d0
          drikdxk=0.0d0 
          drikdyk=0.0d0 
          drikdzk=0.0d0
          drikdxj=0.0d0 
          drikdyj=0.0d0 
          drikdzj=0.0d0
!!
!! check if distance between atom B=k in pair and neighbor i is smaller than cutoff       
          if(rik.le.funccutoffp(i2,iindex))then

            drikdxi=deltaxik/rik
            drikdyi=deltayik/rik
            drikdzi=deltazik/rik
            drikdxk=-1.d0*drikdxi
            drikdyk=-1.d0*drikdyi
            drikdzk=-1.d0*drikdzi
            drikdxj=0.0d0
            drikdyj=0.0d0
            drikdzj=0.0d0
!!                                                          
            fcutik=(0.5d0*(dcos(pi*rik/funccutoffp(i2,iindex))+1.d0))
            dfcutik=0.5d0*(-dsin(pi*rik/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutikdxi=dfcutik*drikdxi
            dfcutikdyi=dfcutik*drikdyi
            dfcutikdzi=dfcutik*drikdzi
            dfcutikdxj=0.0d0
            dfcutikdyj=0.0d0
            dfcutikdzj=0.0d0
            dfcutikdxk=-1.d0*dfcutikdxi
            dfcutikdyk=-1.d0*dfcutikdyi
            dfcutikdzk=-1.d0*dfcutikdzi
!!
            expik= (dexp(-etap(i2,iindex)* (rik)**2))
            dexpikdxi=-2.0d0*expik*deltaxik*etap(i2,iindex)
            dexpikdyi=-2.0d0*expik*deltayik*etap(i2,iindex)
            dexpikdzi=-2.0d0*expik*deltazik*etap(i2,iindex)
            dexpikdxk=-1.0d0*dexpikdxi
            dexpikdyk=-1.0d0*dexpikdyi
            dexpikdzk=-1.0d0*dexpikdzi
            dexpikdxj=0.0d0
            dexpikdyj=0.0d0
            dexpikdzj=0.0d0
          endif

!! calculate distance between atom A=j and B=k in pair 
          rjk=(pairs_xyz(jcount,1,1)-pairs_xyz(jcount,2,1))**2+&
              (pairs_xyz(jcount,1,2)-pairs_xyz(jcount,2,2))**2+&
              (pairs_xyz(jcount,1,3)-pairs_xyz(jcount,2,3))**2
          rjk = sqrt(rjk)
!!
          fcutjk    =0.0d0 ; expjk     =0.0d0
          dfcutjkdxj=0.0d0 ; dfcutjkdyj=0.0d0 ; dfcutjkdzj=0.0d0
          dfcutjkdxk=0.0d0 ; dfcutjkdyk=0.0d0 ; dfcutjkdzk=0.0d0
          dfcutjkdxi=0.0d0 ; dfcutjkdyi=0.0d0 ; dfcutjkdzi=0.0d0

          dexpjkdxj =0.0d0 ; dexpjkdyj =0.0d0 ; dexpjkdzj =0.0d0
          dexpjkdxk =0.0d0 ; dexpjkdyk =0.0d0 ; dexpjkdzk =0.0d0
          dexpjkdxi =0.0d0 ; dexpjkdyi =0.0d0 ; dexpjkdzi =0.0d0

          drjkdxj=0.0d0 ; drjkdyj=0.0d0 ; drjkdzj=0.0d0
          drjkdxk=0.0d0 ; drjkdyk=0.0d0 ; drjkdzk=0.0d0
          drjkdxi=0.0d0 ; drjkdyi=0.0d0 ; drjkdzi=0.0d0


!! check if distance between atom A=j and B=k in pair is smaller than cutoff       
          if(rjk.le.funccutoffp(i2,iindex))then
            drjkdxj=(pairs_xyz(jcount,1,1)-pairs_xyz(jcount,2,1))/rjk
            drjkdyj=(pairs_xyz(jcount,1,2)-pairs_xyz(jcount,2,2))/rjk
            drjkdzj=(pairs_xyz(jcount,1,3)-pairs_xyz(jcount,2,3))/rjk
            drjkdxk=-1.d0*drjkdxj
            drjkdyk=-1.d0*drjkdyj
            drjkdzk=-1.d0*drjkdzj
            drjkdxi=0.0d0
            drjkdyi=0.0d0
            drjkdzi=0.0d0
!!
            fcutjk= (0.5d0*(dcos(pi*rjk/funccutoffp(i2,iindex))+1.d0))
            dfcutjk=0.5d0*(-dsin(pi*rjk/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutjkdxj=dfcutjk*drjkdxj
            dfcutjkdyj=dfcutjk*drjkdyj
            dfcutjkdzj=dfcutjk*drjkdzj
            dfcutjkdxk=-1.d0*dfcutjkdxj
            dfcutjkdyk=-1.d0*dfcutjkdyj
            dfcutjkdzk=-1.d0*dfcutjkdzj
            dfcutjkdxi=0.0d0
            dfcutjkdyi=0.0d0
            dfcutjkdzi=0.0d0
!!
            expjk=(dexp(-etap(i2,iindex)* (rjk)**2))
            dexpjkdxj=-2.0d0*expjk*(pairs_xyz(jcount,1,1)-pairs_xyz(jcount,2,1))*etap(i2,iindex)
            dexpjkdyj=-2.0d0*expjk*(pairs_xyz(jcount,1,2)-pairs_xyz(jcount,2,2))*etap(i2,iindex)
            dexpjkdzj=-2.0d0*expjk*(pairs_xyz(jcount,1,3)-pairs_xyz(jcount,2,3))*etap(i2,iindex)
            dexpjkdxk=-1.0d0*dexpjkdxj
            dexpjkdyk=-1.0d0*dexpjkdyj
            dexpjkdzk=-1.0d0*dexpjkdzj
            dexpjkdxi= 0.0d0
            dexpjkdyi= 0.0d0
            dexpjkdzi= 0.0d0
          endif
!!
!! Calculate Costheta=(C**2-A**2-B**2)/(-2.d0*A*B)
          if(pairs_charge(1,i1).eq.pairs_charge(2,i1))then
!! homonuclear pair
            f =  lambdap(i2,iindex)*(rjk**2 - rij**2 -rik**2)
            g = -2.0d0*rij*rik
            temp4 =  2.0d0**(1.d0-zetap(i2,iindex))

!! Calculate (f/g)' = (f'g - fg')/g^2
            if(ldoforces) then
              dfdxi=lambdap(i2,iindex)*(-2.d0*rij*drijdxi - 2.d0*rik*drikdxi)
              dfdyi=lambdap(i2,iindex)*(-2.d0*rij*drijdyi - 2.d0*rik*drikdyi)
              dfdzi=lambdap(i2,iindex)*(-2.d0*rij*drijdzi - 2.d0*rik*drikdzi)
!!
              dfdxj=lambdap(i2,iindex)*(2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj)
              dfdyj=lambdap(i2,iindex)*(2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj)
              dfdzj=lambdap(i2,iindex)*(2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj)
!!
              dfdxk=lambdap(i2,iindex)*(2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk)
              dfdyk=lambdap(i2,iindex)*(2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk)
              dfdzk=lambdap(i2,iindex)*(2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk)
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
              temp4=(temp4*zetap(i2,iindex)*(1.d0+f/g)**(zetap(i2,iindex)-1))/g**2
              dcosthetadxi=temp4*(dfdxi*g - f*dgdxi)
              dcosthetadyi=temp4*(dfdyi*g - f*dgdyi)
              dcosthetadzi=temp4*(dfdzi*g - f*dgdzi)
              dcosthetadxj=temp4*(dfdxj*g - f*dgdxj)
              dcosthetadyj=temp4*(dfdyj*g - f*dgdyj)
              dcosthetadzj=temp4*(dfdzj*g - f*dgdzj)
              dcosthetadxk=temp4*(dfdxk*g - f*dgdxk)
              dcosthetadyk=temp4*(dfdyk*g - f*dgdyk)
              dcosthetadzk=temp4*(dfdzk*g - f*dgdzk)
            endif ! ldoforces 
!!
          else ! heteronuclear pair
            f= lambdap(i2,iindex)*(rij**2 - rik**2 -rjk**2)
            g= -2.0d0*rjk*rik
            temp4 =  2.0d0**(1.d0-zetap(i2,iindex))

!! Calculate (f/g)' = (f'g - fg')/g^2
            if(ldoforces) then
              dfdxi=lambdap(i2,iindex)*(2.d0*rij*drijdxi - 2.d0*rik*drikdxi)
              dfdyi=lambdap(i2,iindex)*(2.d0*rij*drijdyi - 2.d0*rik*drikdyi)
              dfdzi=lambdap(i2,iindex)*(2.d0*rij*drijdzi - 2.d0*rik*drikdzi)
!!
              dfdxj=lambdap(i2,iindex)*(2.d0*rij*drijdxj - 2.d0*rjk*drjkdxj)
              dfdyj=lambdap(i2,iindex)*(2.d0*rij*drijdyj - 2.d0*rjk*drjkdyj)
              dfdzj=lambdap(i2,iindex)*(2.d0*rij*drijdzj - 2.d0*rjk*drjkdzj)
!!
              dfdxk=lambdap(i2,iindex)*(-2.d0*rik*drikdxk - 2.d0*rjk*drjkdxk)
              dfdyk=lambdap(i2,iindex)*(-2.d0*rik*drikdyk - 2.d0*rjk*drjkdyk)
              dfdzk=lambdap(i2,iindex)*(-2.d0*rik*drikdzk - 2.d0*rjk*drjkdzk)
!!
              dgdxi=-2.d0*rjk*drikdxi
              dgdyi=-2.d0*rjk*drikdyi
              dgdzi=-2.d0*rjk*drikdzi
!!
              dgdxj=-2.d0*drjkdxj*rik
              dgdyj=-2.d0*drjkdyj*rik
              dgdzj=-2.d0*drjkdzj*rik
!!
              dgdxk=-2.d0*(drjkdxk*rik+rjk*drikdxk)
              dgdyk=-2.d0*(drjkdyk*rik+rjk*drikdyk)
              dgdzk=-2.d0*(drjkdzk*rik+rjk*drikdzk)
!!
              temp4=temp4*zetap(i2,iindex)*(1.d0+f/g)**(zetap(i2,iindex)-1)/g**2
              dcosthetadxi=temp4*(dfdxi*g - f*dgdxi)
              dcosthetadyi=temp4*(dfdyi*g - f*dgdyi)
              dcosthetadzi=temp4*(dfdzi*g - f*dgdzi)
              dcosthetadxj=temp4*(dfdxj*g - f*dgdxj)
              dcosthetadyj=temp4*(dfdyj*g - f*dgdyj)
              dcosthetadzj=temp4*(dfdzj*g - f*dgdzj)
              dcosthetadxk=temp4*(dfdxk*g - f*dgdxk)
              dcosthetadyk=temp4*(dfdyk*g - f*dgdyk)
              dcosthetadzk=temp4*(dfdzk*g - f*dgdzk)
            endif
!!
          endif ! homo or heteronuclear pair
!!                
          costheta=f/g
          costheta=1.0d0+costheta
          costheta=(2.0d0**(1.0d0-zetap(i2,iindex)))*(costheta**zetap(i2,iindex))
!!                                                                      
          symfunctionp(i2,i1)=symfunctionp(i2,i1) + ((fcutjk*expjk)* (costheta*fcutik*fcutij*expij*expik))
!!
!! Calculation of Forces
!! (fghijkl)'= f'ghijkl + fg'hijkl + fgh'ijkl + fghi'jkl + fghij'kl + fghijk'l + fghijkl'
!!
!! dsfunc/dx
          if(ldoforces) then
            temp1 = dfcutjkdxi * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdxi * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadxi * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdxi * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdxi * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdxi * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdxi   
!!
            temp2 = dfcutjkdxj * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdxj * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadxj * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdxj * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdxj * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdxj * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdxj
!!
            temp3 = dfcutjkdxk * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdxk * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadxk * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdxk * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdxk * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdxk * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdxk
!!
            dsfuncdxyz_pair(i2,i1,l,1) =dsfuncdxyz_pair(i2,i1,l,1) + temp1
            dsfuncdxyz_pair(i2,jcount,n,1) =dsfuncdxyz_pair(i2,jcount,n,1) + temp2
            dsfuncdxyz_pair(i2,jcount,m,1) =dsfuncdxyz_pair(i2,jcount,m,1) + temp3
          endif
!!
          if(ldostress)then
            strs_pair(1,1,i2,i1)=strs_pair(1,1,i2,i1)+deltaxij*temp2 +deltaxik*temp3
            strs_pair(2,1,i2,i1)=strs_pair(2,1,i2,i1)+deltayij*temp2 +deltayik*temp3
            strs_pair(3,1,i2,i1)=strs_pair(3,1,i2,i1)+deltazij*temp2 +deltazik*temp3
          endif ! ldostress
!!
!! dsfunc/dy
          if(ldoforces) then
            temp1 = dfcutjkdyi * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdyi * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadyi * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdyi * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdyi * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdyi * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdyi
!!
            temp2 = dfcutjkdyj * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdyj * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadyj * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdyj * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdyj * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdyj * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdyj
!! 
            temp3 = dfcutjkdyk * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdyk * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadyk * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdyk * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdyk * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdyk * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdyk
!!
            dsfuncdxyz_pair(i2,i1,l,2) =dsfuncdxyz_pair(i2,i1,l,2) + temp1
            dsfuncdxyz_pair(i2,jcount,n,2) =dsfuncdxyz_pair(i2,jcount,n,2) + temp2
            dsfuncdxyz_pair(i2,jcount,m,2) =dsfuncdxyz_pair(i2,jcount,m,2) + temp3
          endif
!!
          if(ldostress)then
            strs_pair(1,2,i2,i1)=strs_pair(1,2,i2,i1)+deltaxij*temp2 +deltaxik*temp3
            strs_pair(2,2,i2,i1)=strs_pair(2,2,i2,i1)+deltayij*temp2 +deltayik*temp3
            strs_pair(3,2,i2,i1)=strs_pair(3,2,i2,i1)+deltazij*temp2 +deltazik*temp3
          endif ! ldostress
!!
!! dsfunc/dz
          if(ldoforces) then
            temp1 = dfcutjkdzi * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdzi * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadzi * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdzi * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdzi * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdzi * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdzi
!!
            temp2 = dfcutjkdzj * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdzj * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadzj * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdzj * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdzj * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdzj * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdzj
!! 
            temp3 = dfcutjkdzk * expjk     * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * dexpjkdzk * costheta     * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * dcosthetadzk * fcutik     * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * dfcutikdzk * fcutij     * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * dfcutijdzk * expij     * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * dexpijdzk * expik + &
                     fcutjk    * expjk     * costheta     * fcutik     * fcutij     * expij     * dexpikdzk
!!
            dsfuncdxyz_pair(i2,i1,l,3) =dsfuncdxyz_pair(i2,i1,l,3) + temp1
            dsfuncdxyz_pair(i2,jcount,n,3) =dsfuncdxyz_pair(i2,jcount,n,3) + temp2
            dsfuncdxyz_pair(i2,jcount,m,3) =dsfuncdxyz_pair(i2,jcount,m,3) + temp3
          endif
!!
          if(ldostress)then
            strs_pair(1,3,i2,i1)=strs_pair(1,3,i2,i1)+deltaxij*temp2 +deltaxik*temp3
            strs_pair(2,3,i2,i1)=strs_pair(2,3,i2,i1)+deltayij*temp2 +deltayik*temp3
            strs_pair(3,3,i2,i1)=strs_pair(3,3,i2,i1)+deltazij*temp2 +deltazik*temp3
          endif ! ldostress
!!
        endif ! sympelement(i2,1,iindex) ==  pair_lste(i3
      enddo
!!
      return
      end
