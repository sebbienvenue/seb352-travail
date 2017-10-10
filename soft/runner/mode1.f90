!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!! - main.f90
!!
      subroutine mode1(iseed)
!!
      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions
      use mode1options
      use nnshort_atomic
      use nnshort_pair
      use nnham
      use nnewald
      use symfunctions
      use structures
      use timings
!!
      implicit none
!!
      integer iseed                                                        ! in/out
      integer numtrain                                                     ! internal 
      integer numtest                                                      ! internal
      integer numrej                                                       ! internal 
      integer ndim                                                         ! internal 
      integer pointnumber                                                  ! internal
      integer i1,i2                                                        ! internal
      character*200 :: filenametemp                                        ! internal


!!
!!===================================================================================
!! measure total time of mode 1 
!!===================================================================================
      call zerotime(daymode1,timemode1start,timemode1end) 
      call abstime(timemode1start,daymode1)
!!
!!===================================================================================
!! do everything in this subroutine only for non-parallel case
!!===================================================================================
      if(mpisize.gt.1)then
        write(ounit,*)'ERROR: mode 1 is implemented only for non-parallel case'
        stop !'
      endif
!!
!!===================================================================================
!! write output 
!!===================================================================================
      write(ounit,*)'Maximum number of atoms: ',max_num_atoms
      write(ounit,*)'Maximum number of pairs: ',max_num_pairs
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Calculating Symmetry Functions'
      write(ounit,*)'for ',totnum_structures,' structures'
      write(ounit,*)'-------------------------------------------------------------'
!!===================================================================================
!!'
!!===================================================================================
!!  allocate arrays of structures module for a set of nblock structures
!!===================================================================================
      call allocatestructures()
!!
!!===================================================================================
!! initializations
!!===================================================================================
      numtrain                = 0
      numtest                 = 0
      numrej                  = 0
      pointnumber             = 0
!!
!!===================================================================================
!! consistency check if the number of function values  
!! the same as input nodes in input.nn
!!===================================================================================
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_short_atomic(i1).ne.nodes_short_atomic(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues inconsistent with nodes_short_atomic ',element(i1)
            stop !'
          endif
        enddo ! i1
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        do i1=1,npairs
          if(num_funcvalues_short_pair(i1).ne.nodes_short_pair(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues inconsistent with nodes_short_pair ',i1
            stop !'
          endif
        enddo ! i1
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_elec(i1).ne.nodes_elec(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues inconsistent with nodes_elec ',element(i1)
            stop !'
          endif
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(0))then
        do i1=1,npairs
          if(num_funcvalues_ham(i1).ne.nodes_ham(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues_ham inconsistent with nodes_ham ',i1
            stop !'
          endif
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(1))then
        do i1=1,npairs
          if(num_funcvalues_s(i1).ne.nodes_s(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues_s inconsistent with nodes_s ',i1
            stop !'
          endif
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(2))then
        do i1=1,npairs
          if(num_funcvalues_hexton(i1).ne.nodes_hexton(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues_hexton inconsistent with nodes_hexton ',i1
            stop !'
          endif
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(3).and.(mode.gt.1))then
        if(mode.eq.2)then
          ndim=1
        elseif(mode.eq.3)then
          ndim=ntriplets
        endif
        do i1=1,ndim
          if(num_funcvalues_hextoff(i1).ne.nodes_hextoff(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues_hextoff inconsistent with nodes_hextoff ',i1
            stop !'
          endif
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(4))then
        do i1=1,npairs
          if(num_funcvalues_dens(i1).ne.nodes_dens(0,i1))then
            write(ounit,*)'ERROR: num_funcvalues_dens inconsistent with nodes_dens ',i1
            stop !'
          endif
        enddo ! i1
      endif

!!
!!===================================================================================
!! determine maxcutoffs 
!!===================================================================================
      maxcutoff_short_atomic  = 0.0d0
      maxcutoff_short_pair    = 0.0d0
      maxcutoff_elec          = 0.0d0
      maxcutoff_ham           = 0.0d0
      maxcutoff_s             = 0.0d0
      maxcutoff_hexton        = 0.0d0
      maxcutoff_hextoff       = 0.0d0
      maxcutoff_dens          = 0.0d0

      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          do i2=1,num_funcvalues_short_atomic(i1)
            maxcutoff_short_atomic=max(maxcutoff_short_atomic,funccutoff_short_atomic(i2,i1))
          enddo ! i2
        enddo ! i1
      elseif(lshort.and.(nn_type_short.eq.2))then
        do i1=1,npairs
          do i2=1,num_funcvalues_short_pair(i1)
            maxcutoff_short_pair=max(maxcutoff_short_pair,funccutoff_short_pair(i2,i1))
          enddo ! i2
        enddo ! i1
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        do i2=1,nelem
          do i1=1,num_funcvalues_elec(i2)
            maxcutoff_elec=max(maxcutoff_elec,funccutoff_elec(i1,i2))
          enddo ! i1
        enddo ! i2
      endif
      if(lnntb.and.nntb_flag(0))then
        do i1=1,npairs
          do i2=1,num_funcvalues_ham(i1)
            maxcutoff_ham=max(maxcutoff_ham,funccutoff_ham(i2,i1))
          enddo ! i2
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(1))then
        do i1=1,npairs
          do i2=1,num_funcvalues_s(i1)
            maxcutoff_s=max(maxcutoff_s,funccutoff_s(i2,i1))
          enddo ! i2
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(2))then
        do i1=1,npairs
          do i2=1,num_funcvalues_hexton(i1)
            maxcutoff_hexton=max(maxcutoff_hexton,funccutoff_hexton(i2,i1))
          enddo ! i2
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(3))then
        if(mode.lt.3)then
          ndim=1
        else
          ndim=ntriplets
        endif
        do i1=1,ndim
          do i2=1,num_funcvalues_hextoff(i1)
            maxcutoff_hextoff=max(maxcutoff_hextoff,funccutoff_hextoff(i2,i1))
          enddo ! i2
        enddo ! i1
      endif
      if(lnntb.and.nntb_flag(4))then
        do i1=1,npairs
          do i2=1,num_funcvalues_dens(i1)
            maxcutoff_dens=max(maxcutoff_dens,funccutoff_dens(i2,i1))
          enddo ! i2
        enddo ! i1
      endif

!!
!!===================================================================================
!! open files
!!===================================================================================

      open(dataunit,file='input.data',form='formatted')
      rewind(dataunit)
      open(trainstructunit,file='trainstruct.data',form='formatted',status='replace')
      rewind(trainstructunit)
      open(teststructunit,file='teststruct.data',form='formatted',status='replace')
      rewind(teststructunit)
      if(lshort)then
        open(symunit,file='function.data',form='formatted',status='replace')
        rewind(symunit)
        open(tymunit,file='testing.data',form='formatted',status='replace')
        rewind(tymunit)
        open(trainfunit,file='trainforces.data',form='formatted',status='replace')
        rewind(trainfunit)
        open(testfunit,file='testforces.data',form='formatted',status='replace')
        rewind(testfunit)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then ! if we train a separate charge NN:
        open(symeunit,file='functione.data',form='formatted',status='replace')
        rewind(symeunit)
        open(tymeunit,file='testinge.data',form='formatted',status='replace')
        rewind(tymeunit)
        open(trainfeunit,file='trainforcese.data',form='formatted',status='replace')
        rewind(trainfeunit)
        open(testfeunit,file='testforcese.data',form='formatted',status='replace')
        rewind(testfeunit) !'
      endif
      if(lnntb)then
        if(nntb_flag(0))then
!!          write(filenametemp,'(A,I3.3,A,I3.3,A)') &
!!          'function_ham.',tripletmode12(1),'.',tripletmode12(1),'.',tripletmode12(3),'.data'
!!          open(symhamunit,file=filenametemp,form='formatted',status='replace')
!!          rewind(symhamunit)
!!          write(filenametemp,'(A,I4.4,A,I4.4,A)') &
!!          'testing_ham.',tripletmode12(1),'.',tripletmode12(1),'.',tripletmode12(3),'.data'
!!          open(tymhamunit,file=filenametemp,form='formatted',status='replace')
!!          rewind(tymhamunit)
        elseif(nntb_flag(1))then
          open(symoverunit,file='function_s.data',form='formatted',status='replace')
          rewind(symoverunit)
          open(tymoverunit,file='testing_s.data',form='formatted',status='replace')
          rewind(tymoverunit)
        elseif(nntb_flag(2))then
          open(symhextonunit,file='function_hexton.data',form='formatted',status='replace')
          rewind(symhextonunit)
          open(tymhextonunit,file='testing_hexton.data',form='formatted',status='replace')
          rewind(tymhextonunit)
        elseif(nntb_flag(3))then
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
          'function_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          open(symhextoffunit,file=filenametemp,form='formatted',status='replace')
          rewind(symhextoffunit)
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
          'testing_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          open(tymhextoffunit,file=filenametemp,form='formatted',status='replace')
          rewind(tymhextoffunit)
        elseif(nntb_flag(4))then
          open(symdensunit,file='function_dens.data',form='formatted',status='replace')
          rewind(symdensunit)
          open(tymdensunit,file='testing_dens.data',form='formatted',status='replace')
          rewind(tymdensunit)
        endif

      endif
!!
!!===================================================================================
!! calculate and write the symmetry functions 
!!===================================================================================
      call getsymmetryfunctions(iseed,numtrain,numtest,numrej)
!!
!!===================================================================================
!! close files
!!===================================================================================
      close(dataunit)
      close(trainstructunit)
      close(teststructunit)
      if(lshort)then
        close(symunit)
        close(tymunit)
        close(trainfunit)
        close(testfunit)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        close(symeunit)
        close(tymeunit)
        close(trainfeunit)
        close(testfeunit)
      endif
      if(lnntb)then
        if(nntb_flag(0))then
          close(symhamunit)
          close(tymhamunit)
        elseif(nntb_flag(1))then
          close(symoverunit)
          close(tymoverunit)
        elseif(nntb_flag(2))then
          close(symhextonunit)
          close(tymhextonunit)
        elseif(nntb_flag(3))then
          close(symhextoffunit)
          close(tymhextoffunit)
        elseif(nntb_flag(4))then
          close(symdensunit)
          close(tymdensunit)
        endif
      endif
!!
!!===================================================================================
!! write summary to output file
!!===================================================================================
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Number of fitting points: ',numtrain
      write(ounit,*)'Number of testing points: ',numtest
      write(ounit,*)'Number of rejected points:',numrej
!! '
!!===================================================================================
!!    deallocate arrays of structures module
!!===================================================================================
      call deallocatestructures()
!!
      call abstime(timemode1end,daymode1)
!!
      return
      end
