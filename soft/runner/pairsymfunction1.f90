!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconepairfunction.f90
!!
      subroutine pairsymfunction1(i1,i2,jcount,num_atoms,listdim,iindex,&
        pair_lstc,maxnum_funcvaluesp,npairs,max_num_pairs,&
        max_num_atoms,pairs_atom,pair_lste,sympelement,&
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
      integer iindex                                                           ! in, pairindex
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
      integer pair_lste(listdim)                                               ! in
      integer sympelement(maxnum_funcvaluesp,2,npairs)                         ! in
!!
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
      real*8 temp1,temp2,temp3                                           ! internal
      real*8 strs_pair(3,3,maxnum_funcvaluesp,max_num_pairs)                   ! out
      real*8 dsfuncdxyz_pair(maxnum_funcvaluesp,max_num_pairs,max_num_atoms,3) ! out
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs)                    ! out
!!
      logical ldoforces
      logical ldostress

!! loop over all neighbors of pair i1
      do i3 = pair_lsta(i1,1),pair_lsta(i1,2)
!! check if the neighbor has the right nuclear charge for this symmetry function 
        if(sympelement(i2,1,iindex).eq.pair_lste(i3))then
          n = pairs_atom(i1,1)    ! atom A in pair 
          m = pairs_atom(i1,2)    ! atom B in pair
          l = pair_lstc(i3)       ! neighbor
!!
          atomsinpairs(m)=atomsinpairs(m)+1
          atomsinpairs(n)=atomsinpairs(n)+1
!!
!! calculate distance between A=j and neighbor i
          deltaxij=(pair_lstb(i3,1)-pairs_xyz(jcount,1,1))
          deltayij=(pair_lstb(i3,2)-pairs_xyz(jcount,1,2))
          deltazij=(pair_lstb(i3,3)-pairs_xyz(jcount,1,3))
          rij = deltaxij**2 + deltayij**2 + deltazij**2
          rij = sqrt(rij)
!!
          fcutij    =0.0d0 
          fcutjk    =0.0d0
          fcutik    =0.0d0
!!
          drijdxi   =0.0d0 
          drijdyi   =0.0d0 
          drijdzi   =0.0d0
          drijdxj   =0.0d0 
          drijdyj   =0.0d0 
          drijdzj   =0.0d0
          drijdxk   =0.0d0 
          drijdyk   =0.0d0 
          drijdzk   =0.0d0
          dfcutijdxi=0.0d0 
          dfcutijdyi=0.0d0 
          dfcutijdzi=0.0d0
          dfcutijdxj=0.0d0 
          dfcutijdyj=0.0d0 
          dfcutijdzj=0.0d0
          dfcutijdxk=0.0d0 
          dfcutijdyk=0.0d0 
          dfcutijdzk=0.0d0
!!
!! check if rij is small enough for this pair of type iindex
          if(rij.le.funccutoffp(i2,iindex))then
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
            fcutij=0.5d0*(dcos(pi*rij/funccutoffp(i2,iindex))+1.d0)
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
          endif ! rij
!!
!! calculate distance between B=k and neighbor i
          deltaxik=(pair_lstb(i3,1) - pairs_xyz(jcount,2,1))
          deltayik=(pair_lstb(i3,2) - pairs_xyz(jcount,2,2))
          deltazik=(pair_lstb(i3,3) - pairs_xyz(jcount,2,3))
          rik= deltaxik**2 + deltayik**2 + deltazik**2
          rik= sqrt(rik)
!!
          drikdxi   =0.0d0 
          drikdyi   =0.0d0 
          drikdzi   =0.0d0
          drikdxk   =0.0d0 
          drikdyk   =0.0d0 
          drikdzk   =0.0d0
          drikdxj   =0.0d0 
          drikdyj   =0.0d0 
          drikdzj   =0.0d0
          dfcutikdxi=0.0d0 
          dfcutikdyi=0.0d0 
          dfcutikdzi=0.0d0
          dfcutikdxj=0.0d0 
          dfcutikdyj=0.0d0 
          dfcutikdzj=0.0d0
          dfcutikdxk=0.0d0 
          dfcutikdyk=0.0d0 
          dfcutikdzk=0.0d0
!!
!! check if rik is small enough for this pair of type iindex
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
            fcutik=0.5d0*(dcos(pi*rik/funccutoffp(i2,iindex))+1.d0)
            temp1=0.5d0*(-dsin(pi*rik/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutikdxi=temp1*drikdxi
            dfcutikdyi=temp1*drikdyi
            dfcutikdzi=temp1*drikdzi
            dfcutikdxj=0.0d0
            dfcutikdyj=0.0d0
            dfcutikdzj=0.0d0
            dfcutikdxk=-1.d0*dfcutikdxi
            dfcutikdyk=-1.d0*dfcutikdyi
            dfcutikdzk=-1.d0*dfcutikdzi
          endif ! rik
!!
!Third (A & B)
!! calculate distance in pair between A=j and B=k 
          rjk=(pairs_xyz(i1,1,1)-pairs_xyz(jcount,2,1))**2+ &
              (pairs_xyz(i1,1,2)-pairs_xyz(jcount,2,2))**2+ &
              (pairs_xyz(i1,1,3)-pairs_xyz(jcount,2,3))**2
          rjk = sqrt(rjk)

          drjkdxj   =0.0d0 
          drjkdyj   =0.0d0 
          drjkdzj   =0.0d0
          drjkdxk   =0.0d0 
          drjkdyk   =0.0d0 
          drjkdzk   =0.0d0
          drjkdxi   =0.0d0 
          drjkdyi   =0.0d0 
          drjkdzi   =0.0d0
          dfcutjkdxj=0.0d0 
          dfcutjkdyj=0.0d0 
          dfcutjkdzj=0.0d0
          dfcutjkdxk=0.0d0 
          dfcutjkdyk=0.0d0
          dfcutjkdzk=0.0d0
          dfcutjkdxi=0.0d0 
          dfcutjkdyi=0.0d0 
          dfcutjkdzi=0.0d0

!! check if rjk is small enough for this pair of type iindex
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
            fcutjk=0.5d0*(dcos(pi*rjk/funccutoffp(i2,iindex))+1.d0)
            temp1=0.5d0*(-dsin(pi*rjk/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
            dfcutjkdxj=temp1*drjkdxj
            dfcutjkdyj=temp1*drjkdyj
            dfcutjkdzj=temp1*drjkdzj
            dfcutjkdxk=-1.d0*dfcutjkdxj
            dfcutjkdyk=-1.d0*dfcutjkdyj
            dfcutjkdzk=-1.d0*dfcutjkdzj
            dfcutjkdxi=0.0d0
            dfcutjkdyi=0.0d0
            dfcutjkdzi=0.0d0
          endif ! rjk
!!
          symfunctionp(i2,i1)=symfunctionp(i2,i1)+((fcutij*fcutjk)+(fcutik*fcutjk))      
!!
!! dsfunc/dx
          if(ldoforces)then      
            temp1 = (dfcutijdxi*fcutjk) + (fcutij*dfcutjkdxi)+  (fcutik*dfcutjkdxi) + (dfcutikdxi*fcutjk)
            temp2 = (dfcutijdxj*fcutjk) + (fcutij*dfcutjkdxj)+  (fcutik*dfcutjkdxj) + (dfcutikdxj*fcutjk)
            temp3 = (dfcutijdxk*fcutjk) + (fcutij*dfcutjkdxk)+  (fcutik*dfcutjkdxk) + (dfcutikdxk*fcutjk)
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

!! dsfunc/dy
          if(ldoforces)then !! Calculation of derivatives for forces     
            temp1 = (dfcutijdyi*fcutjk) + (fcutij*dfcutjkdyi)+  (fcutik*dfcutjkdyi) + (dfcutikdyi*fcutjk)
            temp2 = (dfcutijdyj*fcutjk) + (fcutij*dfcutjkdyj)+  (fcutik*dfcutjkdyj) + (dfcutikdyj*fcutjk)
            temp3 = (dfcutijdyk*fcutjk) + (fcutij*dfcutjkdyk)+  (fcutik*dfcutjkdyk) + (dfcutikdyk*fcutjk)
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

!! dsfunc/dz
          if(ldoforces)then !! Calculation of derivatives for forces     
            temp1 = (dfcutijdzi*fcutjk) + (fcutij*dfcutjkdzi)+  (fcutik*dfcutjkdzi) + (dfcutikdzi*fcutjk)
            temp2 = (dfcutijdzj*fcutjk) + (fcutij*dfcutjkdzj)+  (fcutik*dfcutjkdzj) + (dfcutikdzj*fcutjk)
            temp3 = (dfcutijdzk*fcutjk) + (fcutij*dfcutjkdzk)+  (fcutik*dfcutjkdzk) + (dfcutikdzk*fcutjk)
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
        endif ! sympelement(i2,1,iindex) == pair_lste(i3)
      enddo ! i3
!!
      return
      end
