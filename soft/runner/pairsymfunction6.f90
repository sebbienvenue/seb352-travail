!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconepairfunction.f90
!!
      subroutine pairsymfunction6(i1,i2,jcount,num_atoms,listdim,iindex,&
        maxnum_funcvaluesp,npairs,max_num_pairs,&
        max_num_atoms,pairs_atom,etap,rshiftp,pairs_rr,&
        pairs_xyz,funccutoffp,strs_pair,&
        dsfuncdxyz_pair,symfunctionp,ldoforces,ldostress)
!!
      use fileunits
      use nnconstants
!!
      implicit none
!!
      integer m,n
      integer i1                                                               ! in, pair
      integer i2                                                               ! in, symfunction
      integer jcount                                                           ! in
      integer listdim                                                          ! in
      integer num_atoms                                                        ! in
      integer pairs_atom(listdim,2)                                            ! in
      integer atomsinpairs(num_atoms)                                          ! internal
      integer maxnum_funcvaluesp                                               ! in
      integer max_num_pairs                                                    ! in
      integer max_num_atoms                                                    ! in
      integer npairs                                                           ! in
      integer iindex                                                           ! in
!!
      real*8 pairs_rr(listdim)                                                 ! in
      real*8 etap(maxnum_funcvaluesp,npairs)                                   ! in
      real*8 rshiftp(maxnum_funcvaluesp,npairs)                                ! in
      real*8 dfcutjk
      real*8 dexpjkdxi,dexpjkdyi,dexpjkdzi
      real*8 dexpjkdxj,dexpjkdyj,dexpjkdzj
      real*8 dexpjkdxk,dexpjkdyk,dexpjkdzk
      real*8 expjk
      real*8 fcutjk
      real*8 pairs_xyz(listdim,2,3)                                            ! in
      real*8 rjk                                                               ! internal
      real*8 drjkdxj, drjkdyj, drjkdzj
      real*8 drjkdxk, drjkdyk, drjkdzk
      real*8 dfcutjkdxi,dfcutjkdyi,dfcutjkdzi
      real*8 dfcutjkdxj,dfcutjkdyj,dfcutjkdzj
      real*8 dfcutjkdxk,dfcutjkdyk,dfcutjkdzk
      real*8 funccutoffp(maxnum_funcvaluesp,npairs)          ! in
      real*8 temp1,temp2                                                       ! internal
      real*8 strs_pair(3,3,maxnum_funcvaluesp,max_num_pairs)                   ! out
      real*8 dsfuncdxyz_pair(maxnum_funcvaluesp,max_num_pairs,max_num_atoms,3) ! out
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs)                    ! out
      real*8 deltaxjk,deltayjk,deltazjk
!!
      logical ldoforces
      logical ldostress
!!
      m = pairs_atom(i1,1)  ! atom A in pair                                                                                     
      n = pairs_atom(i1,2)  ! atom B in pair                                                                                     
!!
      atomsinpairs(m)=atomsinpairs(m)+1
      atomsinpairs(n)=atomsinpairs(n)+1
!!                                                                                                                   
      rjk  = pairs_rr(i1)
!!                                                                                                                    
      if(rjk.le.funccutoffp(i2,iindex))then
!!                                                                                                                    
        fcutjk =(0.5d0*(dcos(pi*rjk/funccutoffp(i2,iindex))+1.d0))
        dfcutjk=0.5d0*(-dsin(pi*rjk/funccutoffp(i2,iindex)))*(pi/funccutoffp(i2,iindex))
!!                                                                                                                    
        deltaxjk=(pairs_xyz(jcount,1,1) - pairs_xyz(jcount,2,1))
        deltayjk=(pairs_xyz(jcount,1,2) - pairs_xyz(jcount,2,2))
        deltazjk=(pairs_xyz(jcount,1,3) - pairs_xyz(jcount,2,3))
        drjkdxj=deltaxjk/rjk
        drjkdyj=deltayjk/rjk
        drjkdzj=deltazjk/rjk
        drjkdxk=-1.d0*drjkdxj
        drjkdyk=-1.d0*drjkdyj
        drjkdzk=-1.d0*drjkdzj
!!                                                                                                                    
        dfcutjkdxj=dfcutjk*drjkdxj
        dfcutjkdyj=dfcutjk*drjkdyj
        dfcutjkdzj=dfcutjk*drjkdzj
        dfcutjkdxk=-1.d0*dfcutjkdxj
        dfcutjkdyk=-1.d0*dfcutjkdyj
        dfcutjkdzk=-1.d0*dfcutjkdzj
        dfcutjkdxi= 0.0d0
        dfcutjkdyi= 0.0d0
        dfcutjkdzi= 0.0d0
!!                                                                                                                    
        expjk=(dexp(-1.d0*etap(i2,iindex)* (rjk-rshiftp(i2,iindex))**2))
!!                                                                                                                    
        dexpjkdxj= (-2.0d0*expjk*deltaxjk* (rjk-rshiftp(i2,iindex)) *etap(i2,iindex))/rjk
        dexpjkdyj= (-2.0d0*expjk*deltayjk* (rjk-rshiftp(i2,iindex)) *etap(i2,iindex))/rjk
        dexpjkdzj= (-2.0d0*expjk*deltazjk* (rjk-rshiftp(i2,iindex)) *etap(i2,iindex))/rjk
        dexpjkdxk=-1.0d0*dexpjkdxj
        dexpjkdyk=-1.0d0*dexpjkdyj
        dexpjkdzk=-1.0d0*dexpjkdzj
        dexpjkdxi= 0.0d0
        dexpjkdyi= 0.0d0
        dexpjkdzi= 0.0d0
!!                                                                                                                    
        symfunctionp(i2,i1)=fcutjk*expjk
!!                                                                                                                    
!! dsfunc/dx                                                                                                          
        if(ldoforces)then !! Calculation of derivatives for forces                                             
          temp1 = (dfcutjkdxj*expjk) + (fcutjk*dexpjkdxj)
          temp2 = (dfcutjkdxk*expjk) + (fcutjk*dexpjkdxk)
!!                                                                                                                    
          dsfuncdxyz_pair(i2,jcount,m,1) =dsfuncdxyz_pair(i2,jcount,m,1) + temp1 
          dsfuncdxyz_pair(i2,jcount,n,1) =dsfuncdxyz_pair(i2,jcount,n,1) + temp2 
!!                                                                                                                    
          if(ldostress)then
            strs_pair(1,1,i2,i1)=strs_pair(1,1,i2,i1)+deltaxjk*temp2
            strs_pair(2,1,i2,i1)=strs_pair(2,1,i2,i1)+deltayjk*temp2
            strs_pair(3,1,i2,i1)=strs_pair(3,1,i2,i1)+deltazjk*temp2
          endif ! ldostress
!! dsfunc/dy
          temp1 = (dfcutjkdyj*expjk) + (fcutjk*dexpjkdyj)
          temp2 = (dfcutjkdyk*expjk) + (fcutjk*dexpjkdyk)
!!
          dsfuncdxyz_pair(i2,jcount,m,2) =dsfuncdxyz_pair(i2,jcount,m,2) + temp1
          dsfuncdxyz_pair(i2,jcount,n,2) =dsfuncdxyz_pair(i2,jcount,n,2) + temp2
!!
          if(ldostress)then
            strs_pair(1,2,i2,i1)=strs_pair(1,2,i2,i1)+deltaxjk*temp2
            strs_pair(2,2,i2,i1)=strs_pair(2,2,i2,i1)+deltayjk*temp2
            strs_pair(3,2,i2,i1)=strs_pair(3,2,i2,i1)+deltazjk*temp2
          endif ! ldostress
!! dsfunc/dz
          temp1 = (dfcutjkdzj*expjk) + (fcutjk*dexpjkdzj)
          temp2 = (dfcutjkdzk*expjk) + (fcutjk*dexpjkdzk)
!!
          dsfuncdxyz_pair(i2,jcount,m,3) =dsfuncdxyz_pair(i2,jcount,m,3) + temp1
          dsfuncdxyz_pair(i2,jcount,n,3) =dsfuncdxyz_pair(i2,jcount,n,3) + temp2
!!
          if(ldostress)then
            strs_pair(1,3,i2,i1)=strs_pair(1,3,i2,i1)+deltaxjk*temp2
            strs_pair(2,3,i2,i1)=strs_pair(2,3,i2,i1)+deltayjk*temp2
            strs_pair(3,3,i2,i1)=strs_pair(3,3,i2,i1)+deltazjk*temp2
          endif ! ldostress
!!
        endif ! ldoforces
      endif  ! rij.le.funccutoff
!!
      return
      end
