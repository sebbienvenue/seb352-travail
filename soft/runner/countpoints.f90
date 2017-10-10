!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose?

!! called by:
!! - fitting.f90
!! - fitting_batch.f90
!! - fittingpair.f90
!!
      subroutine countpoints(iswitch,ndim,&
        ntrain,ntest,ntrainatoms,ntestatoms,&
        train_local,test_local)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer iswitch               ! in
      integer ndim                  ! in
      integer ntrain
      integer ntraine
      integer ntrainatoms
      integer ntrainatomse
      integer train_local(ndim)
      integer test_local(ndim)
      integer train_locale(ndim)
      integer test_locale(ndim)
      integer ntest
      integer nteste
      integer ntestatoms
      integer ntestatomse
      integer num_atoms
      integer ndummy
      integer ndummy1
      integer ndummy2
      integer i1
      integer num_pairs             ! internal
!!
      real*8 edummy
!!
      logical lexist
!!
      character*200 filenametemp         ! internal
!!
      ntrain      =0
      ntraine     =0
      ntrainatoms =0
      ntrainatomse=0
      ntest       =0
      nteste      =0
      ntestatoms  =0
      ntestatomse =0
      train_local(:)=0
      test_local(:) =0
      train_locale(:)=0
      test_locale(:) =0
!!
!! check if files are present
      if(lshort)then
        inquire(file='function.data',exist=lexist)
        if(.not.lexist)then
          write(ounit,*)'Error: function.data not found'
          stop
        endif
        inquire(file='testing.data',exist=lexist)
        if(.not.lexist)then
          write(ounit,*)'Error: testing.data not found'
          stop
        endif
      endif ! lshort
      if(lelec.and.(nn_type_elec.eq.1))then
        inquire(file='functione.data',exist=lexist)
        if(.not.lexist)then
          write(ounit,*)'Error: functione.data not found'
          stop
        endif
        inquire(file='testinge.data',exist=lexist)
        if(.not.lexist)then
          write(ounit,*)'Error: testinge.data not found'
          stop
        endif
      endif ! lelec
      if(lnntb)then
        if(nntb_flag(3))then
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
            'function_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          inquire(file=filenametemp,exist=lexist)
          if(.not.lexist)then
            write(ounit,*)'Error: ',filenametemp,'not found'
            stop
          endif
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
            'testing_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          inquire(file=filenametemp,exist=lexist)
          if(.not.lexist)then
            write(ounit,*)'Error: ',filenametemp,'not found'
            stop
          endif
        endif

      endif
!!
!! open function.data & read
      if(lshort)then
        open(symunit,file='function.data',form='formatted',status='old')
        rewind(symunit)
 10     continue
        if(iswitch.eq.1)then
          read(symunit,*,end=30) num_atoms
          do i1=1,num_atoms
            read(symunit,*)ndummy
            train_local(elementindex(ndummy))=train_local(elementindex(ndummy))+1
          enddo
        elseif(iswitch.eq.2)then
          read(symunit,*,end=30) num_atoms,num_pairs
          do i1=1,num_pairs
            read(symunit,*)ndummy1,ndummy2
            train_local(pairindex(ndummy1,ndummy2))=train_local(pairindex(ndummy1,ndummy2))+1
          enddo
        else
          write(ounit,*)'ERROR unknown iswitch in countpoints ',iswitch
          stop
        endif
        read(symunit,*)edummy
        ntrain=ntrain+1
        ntrainatoms=ntrainatoms+num_atoms   ! out 
        goto 10
 30     continue
        close(symunit)
!!
        open(tymunit,file='testing.data',form='formatted',status='old')
        rewind(tymunit)
 11     continue
        if(iswitch.eq.1)then
          read(tymunit,*,end=31) num_atoms
          do i1=1,num_atoms
            read(tymunit,*)ndummy
            test_local(elementindex(ndummy))=test_local(elementindex(ndummy))+1
          enddo
        elseif(iswitch.eq.2)then
          read(tymunit,*,end=31) num_atoms,num_pairs
          do i1=1,num_pairs
            read(tymunit,*)ndummy1,ndummy2
            test_local(pairindex(ndummy1,ndummy2))=test_local(pairindex(ndummy1,ndummy2))+1
          enddo
        else
          write(ounit,*)'ERROR unknown iswitch in countpoints ',iswitch
          stop
        endif
        read(tymunit,*)edummy
        ntest=ntest+1
        ntestatoms=ntestatoms+num_atoms    
        goto 11
 31     continue
        close(tymunit)
      endif ! lshort
!!
      
      if(lelec.and.(nn_type_elec.eq.1))then
        open(symeunit,file='functione.data',form='formatted',status='old')
        rewind(symeunit) !'
 20     continue
        read(symeunit,*,end=40) num_atoms
        do i1=1,num_atoms
          read(symeunit,*)ndummy
          train_locale(elementindex(ndummy))=train_locale(elementindex(ndummy))+1
        enddo
        read(symeunit,*)edummy
        ntraine=ntraine+1
        ntrainatomse=ntrainatomse+num_atoms   
        goto 20
 40     continue
        close(symeunit)
!!
        open(tymeunit,file='testinge.data',form='formatted',status='old')
        rewind(tymeunit)
 21     continue
        read(tymeunit,*,end=41) num_atoms
        do i1=1,num_atoms
          read(tymeunit,*)ndummy
          test_locale(elementindex(ndummy))=test_locale(elementindex(ndummy))+1
        enddo
        read(tymeunit,*)edummy
        nteste=nteste+1
        ntestatomse=ntestatomse+num_atoms    ! out
        goto 21
 41     continue
        close(tymeunit)
      endif ! lelec
!!
      if(lnntb)then
        if(nntb_flag(3))then
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
            'function_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          open(symhextoffunit,file=filenametemp,form='formatted',status='old')
          rewind(symhextoffunit) !'
 50       continue
          read(symhextoffunit,*,end=60) edummy
          read(symhextoffunit,*,end=70) edummy
          ntrain=ntrain+1
          goto 50
 70       continue
          write(ounit,*) 'Line missing when expected in file,', filenametemp
 60       continue
          close(symhextoffunit)
!!
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
            'testing_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          open(tymhextoffunit,file=filenametemp,form='formatted',status='old')
          rewind(tymhextoffunit)
 51       continue
          read(tymhextoffunit,*,end=61)edummy
          read(tymhextoffunit,*,end=71)edummy
          ntest=ntest+1
          goto 51
 71       continue
          write(ounit,*) 'Line missing when expected in file,', filenametemp
 61       continue
          close(tymhextoffunit)
        endif
      endif ! lnntb
      
!!
      if(lshort.and.(ntest.eq.0))then
        write(ounit,*)'Error: No testing points found in testing.data'
        stop
      endif
      if(lshort.and.(ntrain.eq.0))then
        write(ounit,*)'Error: No training points found in function.data'
        stop
      endif
      if(lelec.and.(nn_type_elec.eq.1).and.(nteste.eq.0))then
        write(ounit,*)'Error: No testing points found in testinge.data'
        stop
      endif
      if(lelec.and.(nn_type_elec.eq.1).and.(ntraine.eq.0))then
        write(ounit,*)'Error: No training points found in functione.data'
        stop
      endif
      if(lnntb)then
        if(nntb_flag(3))then
          if(ntest.eq.0)then
            write(ounit,*)'Error: No training points found in function_hextoff.',&
                           hextoff_training_triplet(1),'.',hextoff_training_triplet(2),&
                           '.',hextoff_training_triplet(3),'.data'
            stop
          endif
          if(ntrain.eq.0)then
            write(ounit,*)'Error: No training points found in training_hextoff.',&
                           hextoff_training_triplet(1),'.',hextoff_training_triplet(2),&
                           '.',hextoff_training_triplet(3),'.data'
            stop
          endif
        endif
      endif
!!
!! make sure that we also know the number of training and testing points if lshort is off
      if((.not.lshort).and.lelec.and.(nn_type_elec.eq.1))then
        ntrain=ntraine
        ntest =nteste
        ntrainatoms    =ntrainatomse
        ntestatoms     =ntestatomse
        train_local(:) =train_locale(:)
        test_local(:)  =test_locale(:)
      endif
!!
      return
      end
