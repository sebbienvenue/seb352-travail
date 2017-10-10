!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calconepairfunction.f90
!!
      subroutine pairsymfunction5(i1,i2,jcount,num_atoms,listdim,&
        maxnum_funcvaluesp,max_num_pairs,&
        max_num_atoms,pairs_atom,pairs_rr,&
        pairs_xyz,strs_pair,&
        dsfuncdxyz_pair,symfunctionp,ldoforces,ldostress)
!!
      use fileunits
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
!!
      real*8 pairs_rr(listdim)                                                 ! in
      real*8 pairs_xyz(listdim,2,3)                                            ! in
      real*8 deltaxij,deltayij,deltazij                                        ! internal
      real*8 rij                                                               ! internal
      real*8 drijdxi, drijdyi, drijdzi                                         ! internal
      real*8 drijdxj, drijdyj, drijdzj                                         ! internal
      real*8 strs_pair(3,3,maxnum_funcvaluesp,max_num_pairs)                   ! out
      real*8 dsfuncdxyz_pair(maxnum_funcvaluesp,max_num_pairs,max_num_atoms,3) ! out
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs)                    ! out
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
      rij  = pairs_rr(i1)
!!                                                                                                                    
      symfunctionp(i2,i1)= rij
!!
      deltaxij=(pairs_xyz(jcount,1,1) - pairs_xyz(jcount,2,1))
      deltayij=(pairs_xyz(jcount,1,2) - pairs_xyz(jcount,2,2))
      deltazij=(pairs_xyz(jcount,1,3) - pairs_xyz(jcount,2,3))
      drijdxi=deltaxij/rij
      drijdyi=deltayij/rij
      drijdzi=deltazij/rij
      drijdxj=-1.d0*drijdxi
      drijdyj=-1.d0*drijdyi
      drijdzj=-1.d0*drijdzi
!!                                                                                                                    
      if(ldoforces)then !! Calculation of derivatives for forces                                             
!!                                                                                                                    
!! dsfunc/dx                                                                                                          
        dsfuncdxyz_pair(i2,jcount,m,1) = drijdxi
        dsfuncdxyz_pair(i2,jcount,n,1) = drijdxj
!!                                                                                                                    
!! dsfunc/dy                                                                                                          
        dsfuncdxyz_pair(i2,jcount,m,2) = drijdyi
        dsfuncdxyz_pair(i2,jcount,n,2) = drijdyj
!!                                                                                                                    
!! dsfunc/dz                                                                                                          
        dsfuncdxyz_pair(i2,jcount,m,3) = drijdzi
        dsfuncdxyz_pair(i2,jcount,n,3) = drijdzj
!!                                                                            
!STRESS                                                                                                               
        if(ldostress)then !! Calculation of derivatives for stress      
!dsfunc/dx                                                                                                          
          strs_pair(1,1,i2,i1) = deltaxij*(drijdxj)
          strs_pair(2,1,i2,i1) = deltayij*(drijdxj)
          strs_pair(3,1,i2,i1) = deltazij*(drijdxj)
!!                                                                                                                    
!! dsfunc/dy                                                                                                          
          strs_pair(1,2,i2,i1) = deltaxij*(drijdyj)
          strs_pair(2,2,i2,i1) = deltayij*(drijdyj)
          strs_pair(3,2,i2,i1) = deltazij*(drijdyj)
!!                                                                                                                    
!! dsfunc/dz                                                                                                          
          strs_pair(1,3,i2,i1) = deltaxij*(drijdzj)
          strs_pair(2,3,i2,i1) = deltayij*(drijdzj)
          strs_pair(3,3,i2,i1) = deltazij*(drijdzj)
        endif ! ldostress                                                                                   
      endif ! ldoforces
!!
      return
      end
