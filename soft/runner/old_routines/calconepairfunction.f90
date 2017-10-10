!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calcpairfunctions.f90
!! - getallshortforcespair.f90
!!
      subroutine calconepairfunction(&
         istruct,num_atoms,zelem,&
         num_pairs,pairs_charge,&
         lattice,xyzstruct,symfunctionp,&
         dsfuncdxyz_pair,strs_pair,&
         lperiodic,ldoforces_local,lrmin)
!!
      use fileunits
      use globaloptions
      use symfunctions
      use nnpair
!!
      implicit none
!!
      integer i1,i2                                                            ! internal
      integer iindex                                                           ! internal
      integer num_atoms                                                        ! in
      integer istruct                                                          ! in
      integer atomsinpairs(num_atoms)                                          ! internal

      integer lsta(2,max_num_atoms)                                            ! out
      integer lstc(listdim)                                                    ! in
      integer lste(listdim)                                                    ! in
      integer num_pairs                                                        ! out
      integer zelem(max_num_atoms)                                             ! in

      integer pairs_charge(2,listdim)                                          ! out 
      integer pair_lstc(listdim)                                               ! in
      integer pairs_atom(listdim,2)                                            ! internal
      integer pair_lste(listdim)                                               ! internal

      real*8 lattice(3,3)                                                      ! in
      real*8 dsfuncdxyz_pair(maxnum_funcvaluesp,max_num_pairs,max_num_atoms,3) ! out
      real*8 lstb(listdim,4) 
      real*8 pi                                                                ! internal
      real*8 strs_pair(3,3,maxnum_funcvaluesp,max_num_pairs)                   ! out
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs)                    ! in
      real*8 xyzstruct(3,max_num_atoms)                                        ! in

      real*8  pairs_rr(listdim)                                                ! in
      real*8  pairs_xyz(listdim,2,3)                                           ! in
      real*8  pair_lsta(listdim,2)                                             ! in
      real*8  pair_lstb(listdim,5)                                             ! in
      real*8  rad2deg

      logical lperiodic
      logical ldoforces_local
      logical lrmin                          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!    initialization
      pi                      = 3.141592654d0
      symfunctionp(:,:)       = 0.0d0
      dsfuncdxyz_pair(:,:,:,:)= 0.0d0
      strs_pair(:,:,:,:)      = 0.0d0
      lsta(:,:)               = 0
      lstb(:,:)               = 0.0d0
      lstc(:)                 = 0
      lste(:)                 = 0
      atomsinpairs(:)         = 0
      rad2deg                 = 180.d0/pi
      lrmin                   =.true.
      num_pairs               = 0
!!
!! Get Independent Pairs & Its Neighbors
     call neighborpair(num_atoms,num_pairs,zelem,&
        lsta,lstb,lstc,lste,maxcutoffp,lattice,xyzstruct,&
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
      do i1=1,num_pairs  ! loop over all pairs
        iindex=pairindex(pairs_charge(1,i1),pairs_charge(2,i1))
        do i2=1,num_funcvaluesp(iindex) ! loop over all symmetry functions of this pair
!! 
          if(function_typep(i2,iindex).eq.1)then  
!!
            call pairsymfunction1(i1,i2,i1,num_atoms,listdim,iindex,&
              pair_lstc,maxnum_funcvaluesp,npairs,max_num_pairs,&
              max_num_atoms,pairs_atom,pair_lste,sympelement,&
              pair_lsta,pair_lstb,pairs_xyz,funccutoffp,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,pi,ldoforces_local,ldostress)
!!
          elseif(function_typep(i2,iindex).eq.2)then ! just cutoff function of distance, no neighbors 
!!
            call pairsymfunction2(i1,i2,i1,num_atoms,listdim,iindex,&
              maxnum_funcvaluesp,npairs,max_num_pairs,&
              max_num_atoms,pairs_atom,pairs_rr,&
              pairs_xyz,funccutoffp,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,pi,ldoforces_local,ldostress)
!!
          elseif(function_typep(i2,iindex).eq.3)then
!!
            call pairsymfunction3(i1,i2,i1,num_atoms,listdim,iindex,&
              pair_lstc,maxnum_funcvaluesp,npairs,max_num_pairs,&
              pair_lste,sympelement,&
              max_num_atoms,pairs_atom,rshiftp,etap,&
              pair_lsta,pair_lstb,pairs_xyz,funccutoffp,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,pi,ldoforces_local,ldostress)
!!
          elseif(function_typep(i2,iindex).eq.4)then
!!
            call pairsymfunction4(i1,i2,i1,num_atoms,listdim,iindex,&
              pair_lstc,maxnum_funcvaluesp,npairs,max_num_pairs,&
              pair_lste,sympelement,&
              max_num_atoms,pairs_atom,pairs_charge,etap,zetap,lambdap,&
              pair_lsta,pair_lstb,pairs_xyz,funccutoffp,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,pi,ldoforces_local,ldostress)
!!
          elseif(function_typep(i2,iindex).eq.5)then ! just distance in pair, no neighbors                               
!!                                                                                                                    
            call pairsymfunction5(i1,i2,i1,num_atoms,listdim,&
              maxnum_funcvaluesp,max_num_pairs,&
              max_num_atoms,pairs_atom,pairs_rr,&
              pairs_xyz,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!                                                                                                                    
          elseif(function_typep(i2,iindex).eq.6)then ! Exp * Fc                                                       
!!                                                                                                                    
            call pairsymfunction6(i1,i2,i1,num_atoms,listdim,iindex,&
              maxnum_funcvaluesp,npairs,max_num_pairs,&
              max_num_atoms,pairs_atom,etap,rshiftp,pairs_rr,&
              pairs_xyz,funccutoffp,strs_pair,&
              dsfuncdxyz_pair,symfunctionp,pi,ldoforces_local,ldostress)
!!
          else
            write(ounit,*)'ERROR: unknown pair symmetry function type in calconepairfunction ',function_typep(i2,iindex)
            stop
          endif ! End type of symmetry fucntions
!!
        enddo ! over symmetry functions
      enddo ! over pairs
!!
      return
      end
