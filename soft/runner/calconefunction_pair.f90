!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate the pair symmetry functions for one structure

!! FIXME: parallel case probably not yet working

!! called by:
!! - optimize_short_combinedpair.f90
!! - predictionpair.f90 
!! - calcpairfunctions.f90
!! - getallshortforcespair.f90
!!
      subroutine calconefunction_pair(n_start,&
         istruct,num_atoms,zelem,&
         num_pairs,pairs_charge,&
         lattice,xyzstruct,symfunctionp,&
         dsfuncdxyz_pair,strs_pair,&
         lperiodic,ldoforces_local,lrmin)
!!
      use fileunits
      use globaloptions
      use symfunctions
      use nnshort_pair
!!
      implicit none
!!
      integer n_start                                                          ! in
      integer i1,i2                                                            ! internal
      integer iindex                                                           ! internal
      integer num_atoms                                                        ! in
      integer jcount  
      integer istruct                                                          ! in, number of structure  

      integer lsta(2,max_num_atoms)                                            ! out
      integer num_pairs                                                        ! out, needed in mode 1
      integer zelem(max_num_atoms)                                             ! in

      integer pairs_charge(2,listdim)                                          ! in
      integer pair_lsta(listdim,2)                                             ! in
      integer pair_lstc(listdim)                                               ! in
      integer pairs_atom(listdim,2)                                            ! in   
      integer pair_lste(listdim)                                               ! internal

      real*8 lattice(3,3)                                                      ! in
      real*8 dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3) ! out
      real*8 lstb(listdim,4) 
      real*8 strs_pair(3,3,maxnum_funcvalues_short_pair,max_num_pairs)                   ! out
      real*8 symfunctionp(maxnum_funcvalues_short_pair,max_num_pairs)                       ! in   check new dimension
      real*8 xyzstruct(3,max_num_atoms)                                        ! in

      real*8  pairs_rr(listdim)                                                 ! in
      real*8  pairs_xyz(listdim,2,3)                                            ! in
      real*8  pair_lstb(listdim,5)                                              ! in

      logical lperiodic
      logical ldoforces_local
      logical lrmin                          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! initialization
      symfunctionp(:,:)       = 0.0d0
      dsfuncdxyz_pair(:,:,:,:)= 0.0d0
      strs_pair(:,:,:,:)      = 0.0d0
      lsta(:,:)               = 0
      lstb(:,:)               = 0.0d0
      lrmin                   =.true.
!!
!! get pairs and neighbors of the pairs 
      call neighborpair_para(num_atoms,num_pairs,zelem,&
        lsta,lstb,maxcutoff_short_pair,lattice,xyzstruct,&
        pairs_rr,pairs_atom,pairs_xyz,pairs_charge,&
        pair_lsta,pair_lstb,pair_lstc,pair_lste,lperiodic)
!!
!! check for obviously wrong structures with too short bonds
!! do this only for the atoms of this process here
      do i1=1,num_atoms
        do i2=lsta(1,i1),lsta(2,i1)
          if(lstb(i2,4).le.rmin)then
            lrmin=.false.
          endif
        enddo
      enddo
!!
!! check for isolated atoms without neighbors in the structure
      do i1=1,num_atoms
        if(lsta(1,i1).eq.0)then
          write(ounit,*)'ERROR: atom has no neighbors ',i1,istruct
          stop
        endif
      enddo
!!
!! Calculate the Symmetry Functions
      jcount=n_start
      do i1=1,num_pairs ! loop over all pairs of this process
!!
!! determine which element pair we have
        iindex=pairindex(pairs_charge(1,jcount),pairs_charge(2,jcount))
!!
        do i2=1,num_funcvalues_short_pair(iindex) ! loop over all symmetry functions of this element pair 

          if(function_type_short_pair(i2,iindex).eq.1)then                 ! 2nd Type of SF
!!
            call pairsymfunction1(i1,i2,jcount,num_atoms,listdim,iindex,&
              pair_lstc,maxnum_funcvalues_short_pair,npairs,max_num_pairs,&
              max_num_atoms,pairs_atom,pair_lste,symelement_short_pair,&
              pair_lsta,pair_lstb,pairs_xyz,funccutoff_short_pair,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!
          elseif(function_type_short_pair(i2,iindex).eq.2)then ! Radial Function (1st Type of SF)
!!
            call pairsymfunction2(i1,i2,jcount,num_atoms,listdim,iindex,&
              maxnum_funcvalues_short_pair,npairs,max_num_pairs,&
              max_num_atoms,pairs_atom,pairs_rr,&
              pairs_xyz,funccutoff_short_pair,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!
          elseif(function_type_short_pair(i2,iindex).eq.3)then
!!
            call pairsymfunction3(i1,i2,jcount,num_atoms,listdim,iindex,&
              pair_lstc,maxnum_funcvalues_short_pair,npairs,max_num_pairs,&
              pair_lste,symelement_short_pair,&
              max_num_atoms,pairs_atom,rshift_short_pair,eta_short_pair,&
              pair_lsta,pair_lstb,pairs_xyz,funccutoff_short_pair,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!
          elseif(function_type_short_pair(i2,iindex).eq.4)then
!!
            call pairsymfunction4(i1,i2,jcount,num_atoms,listdim,iindex,&
              pair_lstc,maxnum_funcvalues_short_pair,npairs,max_num_pairs,&
              pair_lste,symelement_short_pair,&
              max_num_atoms,pairs_atom,pairs_charge,eta_short_pair,zeta_short_pair,lambda_short_pair,&
              pair_lsta,pair_lstb,pairs_xyz,funccutoff_short_pair,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!
          elseif(function_type_short_pair(i2,iindex).eq.5)then ! Radial Function (Only Distance)                               
!!                                                                                                                    
            call pairsymfunction5(i1,i2,jcount,num_atoms,listdim,&
              maxnum_funcvalues_short_pair,max_num_pairs,&
              max_num_atoms,pairs_atom,pairs_rr,&
              pairs_xyz,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!                                                                                                                    
          elseif(function_type_short_pair(i2,iindex).eq.6)then ! Exp * Fc                                                       
!!                                                                                                                    
            call pairsymfunction6(i1,i2,jcount,num_atoms,listdim,iindex,&
              maxnum_funcvalues_short_pair,npairs,max_num_pairs,&
              max_num_atoms,pairs_atom,eta_short_pair,rshift_short_pair,pairs_rr,&
              pairs_xyz,funccutoff_short_pair,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!
          else
           write(ounit,*)'ERROR: Unknown pair symmetry function type ',function_type_short_pair(i2,iindex)
           stop
          endif 
!!
        enddo  ! i2
        jcount=jcount+1
      enddo   ! i1
!!
      return
      end
