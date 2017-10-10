!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by: -getsymmetryfunctions.f90
!!
      subroutine calcpairfunctions(npoints,ndone,&
        iseed,numtrain,numtest,numrej,pointnumber,&
        minvaluee,maxvaluee,avvaluee)
!!
      use fileunits
      use globaloptions 
      use mode1options
      use symfunctions
      use nnshort
      use nnpair
      use nnewald
      use structures
!!
      implicit none
!!
      integer npoints                                                     ! in
      integer num_atoms                                                   ! internal
      integer zelem(max_num_atoms)                                        ! internal
      integer iseed                                                       ! in/out
      integer numtrain                                                    ! in/out
      integer numtest                                                     ! in/out
      integer numrej                                                      ! in/out
      integer ndone                                                       ! in
      integer pointnumber                                                 ! in/out
      integer i1,i2,i3                                                    ! internal
      integer istruct                                                     ! internal
      integer num_pairs                                                   ! internal
      integer pair_charge_list(2,listdim,nblock)                          ! internal 
      integer pairs_charge(2,listdim)                                     ! internal
      integer, allocatable :: lsta(:,:)                          ! numbers of neighbors
      integer, allocatable :: lstc(:)                            ! identification of atom
      integer, allocatable :: lste(:)                            ! nuclear charge of atom
      integer, allocatable :: num_neighbors(:)                            
      integer max_num_neighbors
      integer, allocatable :: neighboridx(:,:)          
      integer, allocatable :: invneighboridx(:,:)  
!!
      real*8 lattice(3,3)                                                 ! internal
      real*8 totalcharge                                                  ! internal
      real*8 totalenergy                                                  ! internal
      real*8 ewaldenergy                                                  ! internal
      real*8 atomcharge(max_num_atoms)                                    ! internal
      real*8 atomenergy(max_num_atoms)                                    ! internal
      real*8 xyzforce(3,max_num_atoms)                                    ! internal
      real*8 ewaldforce(3,max_num_atoms)                                  ! internal
      real*8 ran0                                                         ! internal
      real*8 random                                                       ! internal
      real*8 symfunction(maxnum_funcvalues,max_num_atoms)                 ! internal
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs)               ! internal
      real*8 symfunctione(maxnum_funcvaluese,max_num_atoms)               ! internal
      real*8 symfunctionework(maxnum_funcvaluese,max_num_atoms)           ! internal
      real*8 strs(3,3,maxnum_funcvalues,max_num_atoms)                    ! internal
      real*8 strs_pair(3,3,maxnum_funcvaluesp,max_num_pairs)              ! internal
      real*8 strse(3,3,maxnum_funcvaluese,max_num_atoms)                  ! internal
      real*8, allocatable :: dsfuncdxyz(:,:,:,:)                          ! internal
      real*8, allocatable :: dsfuncdxyz_pair(:,:,:,:)                     ! internal
      real*8, allocatable :: dsfuncdxyze(:,:,:,:)                         ! internal
      real*8, allocatable :: dchargedxyz(:,:,:)                           ! internal
      real*8 minvaluee(nelem,maxnum_funcvaluese)                          ! in 
      real*8 maxvaluee(nelem,maxnum_funcvaluese)                          ! in
      real*8 avvaluee(nelem,maxnum_funcvaluese)                           ! in
      real*8 nnatomcharge(max_num_atoms)                                  ! internal 
      real*8, allocatable :: lstb(:,:)                                  ! xyz and r_ij 

      character*2 elementsymbol(max_num_atoms)                            ! internal

      logical lperiodic                                                   ! internal
      logical ldoforces                                                   ! internal here!!!
      logical lrmin(npoints)                                              ! internal
!!
!! initializations
      dchargedxyz(:,:,:)    = 0.0d0 
      ewaldforce(:,:)       = 0.0d0
      ewaldforce_list(:,:,:)= 0.0d0
      lrmin(:)              = .true.
      ldoforces             = luseforces
      symfunction(:,:)      = 0.0d0
      symfunctionp(:,:)     = 0.0d0
      symfunctione(:,:)     = 0.0d0
!!
!! prepare input for calculating the symmetry functions of one structure
!!
      do i1=1,npoints
        num_atoms        = num_atoms_list(i1)
        zelem(:)         = zelem_list(i1,:) 
        lattice(:,:)     = lattice_list(:,:,i1)
        totalcharge      = totalcharge_list(i1)
        totalenergy      = totalenergy_list(i1)
        atomcharge(:)    = atomcharge_list(i1,:)
        atomenergy(:)    = atomenergy_list(i1,:)
        xyzforce(:,:)    = xyzforce_list(:,:,i1)
        elementsymbol(:) = elementsymbol_list(i1,:)
        lperiodic        = lperiodic_list(i1)
        istruct          = ndone+i1
!!
        allocate(lsta(2,max_num_atoms))
        allocate(lstc(listdim))
        allocate(lste(listdim))
        allocate(lstb(listdim,4))
        allocate(num_neighbors(num_atoms))
!!
        call getneighbors(1,num_atoms,num_atoms,listdim,&
          num_atoms,num_neighbors,zelem,max_num_atoms,max_num_neighbors,&
          lsta,lstc,lste,&
          maxcutoffe,lattice,xyzstruct_list(1,1,i1),lstb,lperiodic)
!!
        if(nn_type.eq.1)then
          allocate(dsfuncdxyz(maxnum_funcvalues,max_num_atoms,0:max_num_neighbors,3))
        elseif(nn_type.eq.2)then
          allocate(dsfuncdxyz_pair(maxnum_funcvaluesp,max_num_pairs,max_num_atoms,3)) 
        endif
        allocate(dsfuncdxyze(maxnum_funcvaluese,max_num_atoms,0:max_num_neighbors,3))
        allocate(neighboridx(num_atoms,0:max_num_neighbors))  
        allocate(invneighboridx(num_atoms,max_num_atoms))  
        call getneighboridx(1,num_atoms,num_atoms,listdim,&
          max_num_atoms,max_num_neighbors,&
          lsta,lstc,neighboridx,invneighboridx)
!!
!! calculate the short range symmetry functions for one structure here 
        if(nn_type.eq.1)then
          call calconefunction_para(cutoff_type,max_num_neighbors,&
            max_num_atoms,1,num_atoms,num_atoms,max_num_atoms,elementindex,&
            maxnum_funcvalues,num_funcvalues, &
            nelem,num_atoms,zelem,listdim,&
            lsta,lstc,lste,invneighboridx,&
            function_type,symelement,lattice,&
            xyzstruct_list(1,1,i1),symfunction,maxcutoff,rmin,&
            funccutoff,eta,rshift,lambda,zeta,dsfuncdxyz,strs,lstb,&
            lperiodic,ldoforces,ldostress,lrmin(i1))
        elseif(nn_type.eq.2)then
          call calconepairfunction_para(1,max_num_pairs,max_num_pairs,&
            istruct,num_atoms,zelem,&
            num_pairs,pairs_charge,&
            lattice,xyzstruct_list(1,1,i1),symfunctionp,&
            dsfuncdxyz_pair,strs_pair,&
            lperiodic,ldoforces,lrmin(i1))
        endif
!!
!! calculate the electrostatic symmetry functions for one structure here 
        call calconefunction_para(cutoff_type,max_num_neighbors,&
          max_num_atoms,1,num_atoms,num_atoms,max_num_atoms,elementindex,&
          maxnum_funcvaluese,num_funcvaluese,&
          nelem,num_atoms,zelem,listdim,&
          lsta,lstc,lste,invneighboridx,&
          function_typee,symeelement,lattice,&
          xyzstruct_list(1,1,i1),symfunctione,maxcutoffe,rmin,&
          funccutoffe,etae,rshifte,lambdae,zetae,dsfuncdxyze,strse,lstb,&
          lperiodic,ldoforces,ldostress,lrmin(i1))
!!
        deallocate(lsta)
        deallocate(lstc)
        deallocate(lste)
        deallocate(lstb)
        deallocate(invneighboridx)  
        deallocate(dsfuncdxyz_pair)  
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the electrostatic contribution to the total energy and forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        if(lelec)then
          allocate(dchargedxyz(max_num_atoms,0:max_num_neighbors,3)) 
          dchargedxyz(:,:,:)=0.0d0
!! fixed charges (not environment-dependent)
          if(lfixedcharges)then
            if(lperiodic) then
              call getewaldenergy(max_num_neighbors,&
                neighboridx,num_neighbors,num_atoms,zelem,&
                lattice,xyzstruct_list(1,1,i1),atomcharge,ewaldenergy,&
                dchargedxyz,ewaldforce,.false.)
            else ! not periodic
!! calculate electrostatic energy for non-periodic system
              call electrostatic(num_atoms,&
                atomcharge,xyzstruct_list(1,1,i1),ewaldenergy)
!! calculate electrostatic force for non-periodic system and fixed charges
              if(ldoforces)then
                call splitcoulombforces(&
                  num_atoms,atomcharge,xyzstruct_list(1,1,i1),ewaldforce)
              endif ! ldoforces
            endif ! lperiodic
!! environment-dependent charges:
          else
!!
!! We need to get the charges and forces from a prepared NN fit here!
!!
!! Step 1: prepare symmetry functions
!! for the scaling we have to make a working copy of symfuncione -> symfunctionework
!! (symfunctione itself must be written to file unscaled!)
!! for dsfuncdxyze this is not necessary, because it is not needed later and can be modified
            symfunctionework(:,:)=symfunctione(:,:)
!!
!! scale the symmetry functions for the charge prediction
            call scalesymone(nelem,&
              maxnum_funcvaluese,num_funcvaluese,num_atoms,&
              zelem,symfunctionework,&
              minvaluee,maxvaluee,avvaluee,&
              scmin_ewald,scmax_ewald)
!!
!! we also need to scale the derivative terms dsfuncdxyze and strse 
            call scaledsfunc(max_num_neighbors,&
              maxnum_funcvaluese,num_funcvaluese,&
              nelem,num_atoms,minvaluee,maxvaluee,&
              scmin_ewald,scmax_ewald,&
              zelem,dsfuncdxyze,strse)
!!
!! predict the atomic charges 'nnatomcharge' for this structure
            nnatomcharge(:)=0.0d0
            call calconecharge(num_atoms,&
              zelem,symfunctionework,nnatomcharge)
!!
!! Step 2: calculate dchargedxyz array for forces
            if(lperiodic)then
              call getdchargedxyz(max_num_neighbors,&
                num_neighbors,neighboridx,num_atoms,zelem,&
                dsfuncdxyze,dchargedxyz,symfunctionework)
            else ! not periodic
              call getcoulombdchargedxyz(max_num_neighbors,&
                num_neighbors,neighboridx,num_atoms,zelem,&
                dsfuncdxyze,symfunctionework,dchargedxyz)
            endif ! lperiodic
!!
!! FIXME: Dirty workaround in case we are intending to do the first charge fit without short range part
!! (if we want to fit short range part then NN electrostatic E and F should be removed, but otherwise not)
            if(.not.lshort)then
              nnatomcharge(:)=atomcharge(:)
            endif
!!
!! Step 3: calculate the NN electrostatic energy and forces
            if(lperiodic) then
              call getewaldenergy(max_num_neighbors,&
                neighboridx,num_neighbors,num_atoms,zelem,&
                lattice,xyzstruct_list(1,1,i1),nnatomcharge,ewaldenergy,&
                dchargedxyz,ewaldforce,.false.)
            else ! not periodic
!! calculate ewald energy for non-periodic system
              call electrostatic(num_atoms,&
                nnatomcharge,xyzstruct_list(1,1,i1),ewaldenergy)
!! calculate ewald force for non-periodic system and fixed charges
              if(ldoforces)then
                call getcoulombforcesone(max_num_neighbors,&
                  num_atoms,dchargedxyz,&
                  nnatomcharge,xyzstruct_list(1,1,i1),ewaldforce)
              endif ! ldoforces
            endif ! lperiodic
!!
          endif ! lfixedcharges
!!
          ewaldforce_list(:,:,i1) =ewaldforce(:,:)
          shortforce_list(:,:,i1) =xyzforce_list(:,:,i1)-ewaldforce_list(:,:,i1)
          deallocate(dchargedxyz) 
!!
        else ! no electrostatics
          ewaldenergy=0.0d0
          shortforce_list(:,:,i1) =xyzforce_list(:,:,i1)
        endif ! lelec
        deallocate(num_neighbors)
        deallocate(neighboridx)
        deallocate(dsfuncdxyze)
!!
        symfunctionp_list(:,:,i1)= symfunctionp(:,:)
        symfunctione_list(:,:,i1)= symfunctione(:,:)
        ewaldenergy_list(i1)     = ewaldenergy
        shortenergy_list(i1)     = totalenergy_list(i1)-ewaldenergy_list(i1)
        num_pairs_list(i1)       = num_pairs
        pair_charge_list(:,:,i1) = pairs_charge(:,:)
      enddo ! i1 over npoints
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write points to files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      do i1=1,npoints
        pointnumber=pointnumber+1
!!
!! decide if the point should be used (no high energies and no too close atoms)
        if(((lfitethres.and.(totalenergy_list(i1).lt.(num_atoms_list(i1)*fitethres)))&
          .or.(.not.lfitethres)).and.lrmin(i1))then
!!
!! get random number for splitting in training and test set
          random=ran0(iseed)
!!
!! check if write format statements are sufficient for the number of short range symmetry functions
          if(lshort)then
            do i2=1,npairs
              if(num_funcvaluesp(i2).gt.500)then
                write(ounit,*)'Error: only 500 funcvalues possible'
                stop
              endif
            enddo
          endif ! lshort
          if(lelec.and.(etype.eq.1))then
            do i2=1,nelem
              if(num_funcvaluese(i2).gt.500)then
                write(ounit,*)'Error: only 500 funcvaluese possible'
                stop
              endif
            enddo
          endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(random.gt.splitthres)then ! this is a training point
!!
!! write function.data
            if(lshort)then
              call writepairsymfunctions(symunit,i1,pair_charge_list,&
                maxnum_funcvaluesp,num_funcvaluesp,symfunctionp_list)
            endif ! lshort
!!
!! write functione.data
            if(lelec.and.(etype.eq.1).and.(.not.lfixedcharges))then
              call writeatomicsymfunctions(symeunit,i1,&
                maxnum_funcvaluese,num_funcvaluese,symfunctione_list)
            endif ! lelec
!!
!! write trainstruct.data
!! CHANGE ANDI: GFORTRAN: gfortran has different default width for logicals (1 instead of 2 for ifort),
!!                        this results in a missing space between i8 and l -> problem when reading file,
!!                        so I inserted an additional space. 
           !write(trainstructunit,'(i8,l)')numtrain+1,lperiodic_list(i1)
            write(trainstructunit,'(i8,tr1,l)')numtrain+1,lperiodic_list(i1)
!! END CHANGE

            if(lperiodic_list(i1))then
              do i2=1,3
                write(trainstructunit,'(3f20.14)')(lattice_list(i2,i3,i1),i3=1,3)
              enddo
            endif
            do i2=1,num_atoms_list(i1)
!! Here we write the total forces 
              write(trainstructunit,'(i3,8f15.10)')zelem_list(i1,i2),&
                (xyzstruct_list(i3,i2,i1),i3=1,3),atomcharge_list(i1,i2),&
                atomenergy_list(i1,i2),(xyzforce_list(i3,i2,i1),i3=1,3)
            enddo ! i2
!!
            if(luseforces)then
!!
!! write trainforces.data
              if(lshort)then
                write(trainfunit,'(i8)')numtrain+1
                do i2=1,num_atoms_list(i1)
                  write(trainfunit,'(3f15.10)')(shortforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lshort
!!
!! write trainforcese.data
!! FOR FITTING WE DO NOT NEED THIS FILE
              if(lelec.and.(.not.lfixedcharges))then
                write(trainfeunit,'(i8)')numtrain+1
                do i2=1,num_atoms_list(i1)
                  write(trainfeunit,'(3f15.10)')(ewaldforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lelec
!!
            endif ! luseforces
!!
            numtrain=numtrain+1
            write(ounit,*)pointnumber,' Point is used for training ',numtrain
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else ! point is in test set
!!
!! write testing.data
            if(lshort)then
              call writepairsymfunctions(tymunit,i1,pair_charge_list,&
                maxnum_funcvaluesp,num_funcvaluesp,symfunctionp_list)
            endif ! lshort
!!
!! write testinge.data
            if(lelec.and.(etype.eq.1).and.(.not.lfixedcharges))then
              call writeatomicsymfunctions(tymeunit,i1,&
                maxnum_funcvaluese,num_funcvaluese,symfunctione_list)
            endif ! lelec
!!
!! write teststruct.data
!! CHANGE ANDI: GFORTRAN: gfortran has different default width for logicals (1 instead of 2 for ifort),
!!                        this results in a missing space between i8 and l -> problem when reading file,
!!                        so I inserted an additional space. 
           !write(teststructunit,'(i8,l)')numtest+1,lperiodic_list(i1)
            write(teststructunit,'(i8,tr1,l)')numtest+1,lperiodic_list(i1)
!! END CHANGE

            if(lperiodic_list(i1))then
              do i2=1,3
                write(teststructunit,'(3f20.14)')(lattice_list(i2,i3,i1),i3=1,3)
              enddo
            endif ! lperiodic
            do i2=1,num_atoms_list(i1)
              write(teststructunit,'(i3,8f15.10)')zelem_list(i1,i2),&
                (xyzstruct_list(i3,i2,i1),i3=1,3),atomcharge_list(i1,i2),&
                atomenergy_list(i1,i2),(xyzforce_list(i3,i2,i1),i3=1,3)
            enddo ! i2
!!
            if(luseforces)then
!!
!! write testforces.data
              if(lshort)then
                write(testfunit,'(i8)')numtest+1
                do i2=1,num_atoms_list(i1)
                  write(testfunit,'(3f15.10)')(shortforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lshort
!!
!! write testforcese.data
!! FOR FITTING WE DO NOT NEED THIS FILE
              if(lelec.and.(.not.lfixedcharges))then
                write(testfeunit,'(i8)')numtest+1
                do i2=1,num_atoms_list(i1)
                  write(testfeunit,'(3f15.10)')(ewaldforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lelec
            endif ! luseforces
!!
            numtest=numtest+1
            write(ounit,*)pointnumber,' Point is used for testing ',numtest
          endif ! random
!!
        else ! point is rejected because of lfitethres or lrmin
!!
          numrej=numrej+1
          if(.not.lrmin(i1))then
            write(ounit,*)pointnumber,' Point is rejected (too short bond) ',numrej
          else !'
            write(ounit,*)pointnumber,' Point is rejected (high E) ',numrej
          endif
!!
        endif ! lfitethres
!!
      enddo ! i1 loop over all structures
!!
      return
!!
      end
