!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate the pair symmetry functions for one structure

!! FIXME: parallel case probably not yet working

!! called by:
!! - optimize_short_combinedpair.f90
!! - prediction_hextoff.f90 
!! - calcpairfunctions.f90
!! - getallshortforcespair.f90
!!
!! CMH
!! Unlike the other versions, this section requires the preprocessing of the structures for there
!! to be a prediction performed/analysis of the structure for training.
!! For training we have preprocessed already by having all training examples in the correct orientation
!! but for mode 3 this will be different. 
!!
!!
      subroutine calconefunction_hextoff(n_start,&
         istruct,num_atoms,zelem,&
         num_pairs,pairs_charge,&
         lattice,xyzstruct,symfunction_hextoff,&
         dsfuncdxyz_hextoff,strs_hextoff,&
         lperiodic,ldoforces_local,lrmin)
!!
      use fileunits
      use globaloptions
      use symfunctions
      use nnham
      use nnflags
!!
      implicit none
!!
      integer n_start                                                          ! in
      integer i1,i2,i3                                                         ! internal
      integer iindex,jindex                                                    ! internal
      integer num_atoms                                                        ! in
      integer jcount  
      integer istruct                                                          ! in, number of structure  

      integer lsta(2,max_num_atoms)                                            ! out
      integer num_pairs                                                        ! out, needed in mode 1
      integer zelem(max_num_atoms)                                             ! in

      integer pairs_charge(2,listdim)                                          ! in
!!      integer triplet_charge(2,listdim)                                        ! in
      integer pair_lsta(listdim,2)                                             ! in
      integer pair_lstc(listdim)                                               ! in
      integer pairs_atom(listdim,2)                                            ! in   
      integer pair_lste(listdim)                                               ! internal
      integer tripletstore((max_num_atoms*(max_num_atoms-1)),3)


      real*8 lattice(3,3)                                                      ! in
      real*8 dsfuncdxyz_hextoff(maxnum_funcvalues_hextoff,max_num_atoms,max_num_triplets,3) ! out
      real*8 lstb(listdim,4) 
      real*8 pi                                                                ! internal
      real*8 strs_hextoff(3,3,maxnum_funcvalues_hextoff,max_num_pairs)                   ! out
      real*8 symfunction_hextoff(maxnum_funcvalues_hextoff,max_num_pairs)                       ! in   check new dimension
      real*8 xyzstruct(3,max_num_atoms)                                        ! in

      real*8  pairs_rr(listdim)                                                 ! in
      real*8  pairs_xyz(listdim,2,3)                                            ! in
      real*8  pair_lstb(listdim,5)                                              ! in
      real*8 rad2deg
      integer calc_flag                                                         ! in
      logical lperiodic
      logical ldoforces_local
      logical lrmin                          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! initialization
      pi                      = 3.141592654d0
      symfunction_hextoff(:,:)       = 0.0d0
      dsfuncdxyz_hextoff(:,:,:,:)= 0.0d0
      strs_hextoff(:,:,:,:)      = 0.0d0
      lsta(:,:)               = 0
      lstb(:,:)               = 0.0d0
      rad2deg                 = 180.d0/pi
      lrmin                   =.true.
!!
!! get pairs and neighbors of the pairs 
      call neighborpair_para(num_atoms,num_pairs,zelem,&
        lsta,lstb,maxcutoff_hextoff,lattice,xyzstruct,&
        pairs_rr,pairs_atom,pairs_xyz,pairs_charge,&
        pair_lsta,pair_lstb,pair_lstc,pair_lste,lperiodic)
!! get triplets and indentify the unique triplets
!!
      !!write(*,*) 'pair index', pairs_atom(1,:), pairs_atom(2,:), pairs_atom(3,:)
      !!write(*,*) 'pair_lstc', pair_lstc(1:12)
      !!write(*,*)
      !! CMH
      call neighbortriplet_para(num_atoms,num_pairs,zelem,&
          xyzstruct,maxcutoff_hextoff,lperiodic,tripletstore)
!!        lsta,lstb,maxcutoff_hextoff,lattice,xyzstruct,&
!!        pairs_rr,pairs_atom,pairs_xyz,pairs_charge,&
!!        pair_lsta,pair_lstb,pair_lstc,pair_lste,lperiodic)
      !!write(*,*) tripletstore(1,:)
      !!write(*,*) storetriplets(1,:)
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
!!  CMH we need to check what pairs are found. We only need to loop over certain pairs and (so the origin, and z atom)
!! We need to identify those atoms, for use and then find the atoms they make symmetry functions with other atoms in the triplet.
!! so looping is very different and no existent in the determination of the symmetry function
!!
!!      write(*,*) 'Checking calconefunctions_hextoff' !! CMH


!!
!! Calculate the Symmetry Functions
      jcount=n_start
      calc_flag = 0
      do i1=1,num_pairs ! loop over all pairs of this process
!!
!! determine which element pair we have !!! CMH This is now going to be the triplet pair. The important thing to note that the triplet
!!                                            is only recognized in the specific order of origin, zaxis, neighbour
        iindex=pairindex(pairs_charge(1,jcount),pairs_charge(2,jcount))
!!        write(*,*) pair_lste(jcount)
!!
!!      CMH given we now have determined the pair we are looking at we need to loop over all the triplets that they form
!!      of course within the given cut off
!!        write(*,*) size(pair_lste)
!!
!!      Loop over the number of triplets (thus number of triplets each pair forms with their neighbours)
!!
        do i3=1,size(pair_lste) 
!!          write(*,*) calc_flag, 'calc_flag'
          if(pair_lste(i3).eq.0) then
          else
          
            jindex=tripletindex(pairs_charge(1,jcount),pairs_charge(2,jcount),pair_lste(i3))
        
!!            write(*,*) jindex
            
!!        STOP
!!
!!
!!        If mode 1 then check for centrality of the pair - we have no preprocessing to do to rotation of the triplet
!!
            if(mode.eq.1) then
!!              write(*,*) i1, 'i1', iindex, 'iindex', i3, 'i3', jindex, 'jindex'
!!              write(*,*) pairs_charge(1,jcount), 'atom 1 charge', pairs_charge(2,jcount), 'atom 2 charge'
!!              write(*,*) pair_lste(i3)
!!              write(*,*) pairs_xyz(jcount,1,:), 'atom 1', pairs_xyz(jcount,2,:), 'atom 2'
              if(abs((pairs_xyz(jcount,1,1)).le.(0.0d-10))&
                  .and.(abs(pairs_xyz(jcount,1,2)).le.(0.0d-10))&
                  .and.(abs(pairs_xyz(jcount,1,3)).le.(0.0d-10))&
                  .and.(abs(pairs_xyz(jcount,2,1)).le.(0.0d-10))&
                  .and.(abs(pairs_xyz(jcount,2,2)).le.(0.0d-10))) then
!!                write(*,*) 'This pair is central'
!!                write(*,*) i1, 'i1', iindex, 'iindex', i3, 'i3', jindex, 'jindex'
!!                write(*,*) pairs_atom(jcount,1), 'atom 1', pairs_atom(jcount,2), 'atom 2'
!!                write(*,*) pairs_charge(1,jcount), 'atom 1 charge', pairs_charge(2,jcount), 'atom 2 charge'
!!                write(*,*) pair_lste(i3)
!!                write(*,*) pairs_xyz(jcount,1,:), 'atom 1', pairs_xyz(jcount,2,:), 'atom 2'
!!                write(*,*) pair_lsta(1:3,1)
!!                write(*,*) pair_lsta(1:3,2)
!!                write(*,*) pair_lstb(1,1:3)
!!                write(*,*) pair_lstc(1:3)
!!                write(*,*) tripletmode12
                 if((pairs_charge(1,jcount).eq.hextoff_training_triplet(1))&
                     .and.(pairs_charge(2,jcount).eq.hextoff_training_triplet(2))&
                     .and.(pair_lste(i3).eq.hextoff_training_triplet(3))) then
!!                   write(*,*) pairs_charge(1,n_start), 'atom 1 charge', pairs_charge(2,n_start), 'atom 2 charge', pair_lste(i3), 'atom 3 charge'
!!                   PAUSE
                   calc_flag = 1
                 endif
              endif
            elseif(mode.eq.2) then
!!
!!       If mode 2 we need to preprocess the xyz coords and bring the pair and neighbour into alignment for calculation
!!          calc_flag = 1
            endif
!!
!! We now know which of our pairs in the system is the central one, and of course we we know which is the neighbor to the pair - because we only have
!! three atoms in the system.
!!
!! So lets now find the function values
!!
            if(calc_flag.eq.1) then
              do i2=1,num_funcvalues_hextoff(jindex) ! loop over all symmetry functions of this element pair 
!!
                if(function_type_hextoff(i2,jindex).eq.1)then                 ! 2nd Type of SF
                  write(*,*) 'function 1 starts'
!!
                  call hextoffsymfunction1(i1,i2,jcount,num_atoms,listdim,iindex,&
                    jindex,&
                    pair_lstc,maxnum_funcvalues_hextoff,npairs,ntriplets,max_num_pairs,&
                    max_num_atoms,pairs_atom,pair_lste,symelement_hextoff,&
                    pair_lsta,pair_lstb,pairs_xyz,eta_hextoff,funccutoff_hextoff,strs_hextoff,&
                    dsfuncdxyz_hextoff,symfunction_hextoff,pi,ldoforces_local,ldostress)
!!
                elseif(function_type_hextoff(i2,jindex).eq.2)then ! Radial Function (1st Type of SF)
                  write(*,*) 'function 2 starts'

!!
!!              call hextoffsymfunction2(i1,i2,jcount,num_atoms,listdim,iindex,&
!!                jindex,&
!!                pair_lstc,maxnum_funcvalues_hextoff,npairs,ntriplets,max_num_pairs,&
!!                max_num_atoms,pairs_atom,pair_lste,symelement_hextoff,&
!!                pair_lsta,pair_lstb,pairs_xyz,eta_hextoff,funccutoff_hextoff,strs_hextoff,&
!!                dsfuncdxyz_hextoff,symfunction_hextoff,pi,ldoforces_local,ldostress)
!!
                elseif(function_type_hextoff(i2,jindex).eq.3)then
                write(*,*) 'function 3 starts'
!!
!!
!!             call hextoffsymfunction3(i1,i2,jcount,num_atoms,listdim,iindex,&
!!                jindex,&
!!                pair_lstc,maxnum_funcvalues_hextoff,npairs,ntriplets,max_num_pairs,&
!!                max_num_atoms,pairs_atom,pair_lste,symelement_hextoff,&
!!                pair_lsta,pair_lstb,pairs_xyz,eta_hextoff,funccutoff_hextoff,strs_hextoff,&
!!                dsfuncdxyz_hextoff,symfunction_hextoff,pi,ldoforces_local,ldostress)
!!
!!          elseif(function_type_short_pair(i2,iindex).eq.4)then
!!
!!            call pairsymfunction4(i1,i2,jcount,num_atoms,listdim,iindex,&
!!              pair_lstc,maxnum_funcvalues_short_pair,npairs,max_num_pairs,&
!!              pair_lste,symelement_short_pair,&
!!              max_num_atoms,pairs_atom,pairs_charge,eta_short_pair,zeta_short_pair,lambda_short_pair,&
!!              pair_lsta,pair_lstb,pairs_xyz,funccutoff_short_pair,strs_pair,&
!!              dsfuncdxyz_pair,symfunctionp,pi,ldoforces_local,ldostress)
!!
!!          elseif(function_type_short_pair(i2,iindex).eq.5)then ! Radial Function (Only Distance)                               
!!                                                                                                                    
!!           call pairsymfunction5(i1,i2,jcount,num_atoms,listdim,&
!!              maxnum_funcvalues_short_pair,max_num_pairs,&
!!              max_num_atoms,pairs_atom,pairs_rr,&
!!              pairs_xyz,strs_pair,&
!!              dsfuncdxyz_pair,symfunctionp,ldoforces_local,ldostress)
!!                                                                                                                    
!!          elseif(function_type_short_pair(i2,iindex).eq.6)then ! Exp * Fc                                                       
!!                                                                                                                    
!!            call pairsymfunction6(i1,i2,jcount,num_atoms,listdim,iindex,&
!!              maxnum_funcvalues_short_pair,npairs,max_num_pairs,&
!!              max_num_atoms,pairs_atom,eta_short_pair,rshift_short_pair,pairs_rr,&
!!              pairs_xyz,funccutoff_short_pair,strs_pair,&
!!              dsfuncdxyz_pair,symfunctionp,pi,ldoforces_local,ldostress)
!!
                else
                  write(ounit,*)'ERROR: Unknown pair symmetry function type ',function_type_short_pair(i2,iindex)
                  stop
                endif              
              enddo !! i2
!!
            endif
          endif
          calc_flag = 0
        enddo ! i3
        
        jcount = jcount + 1
      enddo   ! i1
!!
      write(*,*) 'End calconefunction_hextoff'
      STOP
      return
      end
