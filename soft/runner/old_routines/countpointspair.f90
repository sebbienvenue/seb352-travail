!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine countpointspair(nelem,npairs,pairindex,&
        ntrain,ntest,ntrainatoms,ntestatoms,&
        trainelemp,testelemp,ldebug)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2                 ! internal
      integer nelem                 ! in
      integer npairs                ! in  
      integer ntrain                ! out
      integer ntrainatoms           ! out
      integer trainelemp(npairs)    ! out
      integer testelemp(npairs)     ! out
      integer ntest                 ! out
      integer ntestatoms            ! out
      integer num_atoms             ! internal
      integer nsizepair             ! internal
      integer ndummy1               ! internal
      integer ndummy2               ! internal 
      integer pairindex(102,102)    ! in
!!
      real*8 rdummy                 ! internal
!!
      logical lexist              
      logical ldebug                ! in
!!
      ntrain=0
      ntrainatoms=0
      ntest=0
      ntestatoms=0
      trainelemp(:)=0
      testelemp(:)=0
!!
!! check if files are present
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

!! open function.data & read
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)
 10   continue
      read(symunit,*,end=30) num_atoms,nsizepair

      do i1=1,nsizepair 
        read(symunit,*)ndummy1,ndummy2
        trainelemp(pairindex(ndummy1,ndummy2))=trainelemp(pairindex(ndummy1,ndummy2))+1
      enddo
      read(symunit,*)rdummy

      ntrain=ntrain+1
      ntrainatoms=ntrainatoms+num_atoms   ! out 
      goto 10
 30   continue
      close(symunit)

!! open testing.data & read
      open(symunit,file='testing.data',form='formatted',status='old')
      rewind(symunit)
 11   continue
      read(symunit,*,end=31) num_atoms,nsizepair
      do i1=1,nsizepair  
        read(symunit,*)ndummy1,ndummy2
        testelemp(pairindex(ndummy1,ndummy2))=testelemp(pairindex(ndummy1,ndummy2))+1
      enddo
      read(symunit,*)rdummy
      ntest=ntest+1
      ntestatoms=ntestatoms+num_atoms    ! out
      goto 11
 31   continue
      close(symunit)

!! If the files are absent !!
      if(ntest.eq.0)then
        write(ounit,*)'Error: No testing points found in testing.data'
        stop
      endif
      if(ntrain.eq.0)then
        write(ounit,*)'Error: No training points found in function.data'
        stop
      endif
!!
      return
      end
