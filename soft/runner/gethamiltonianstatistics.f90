!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fitting_batch.f90
!! - fittingpair.f90
!!
      subroutine gethamiltonianstatistics(iswitch,&
        hamav,hamstddev,hammin,hammax,matrixsize)
!!
      use fileunits
      use nnflags
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer i1
      integer iswitch               ! in
      integer matrixsize            ! in
      integer natoms                ! internal
      integer ndim                  ! internal
      integer nsizepair             ! internal
      integer idummy                ! internal
      integer nstruct               ! internal
!!
      real*8 hammin(matrixsize)              ! out
      real*8 hammax(matrixsize)              ! out
      real*8 hamav(matrixsize)               ! out
      real*8 hamstddev(matrixsize)           ! out 
      real*8 ham(matrixsize)                    ! internal
      real*8 tounit                 ! internal
      real*8 zdummy                 ! internal
      real*8 tempsumham(matrixsize)          ! internal
      character*100 filenametemp      ! internal
!!
!!
!! initializations
      hammin            = 1000000000.d0
      hammax            = -1000000000.0d0
      hamav             = 0.0d0
      hamstddev         = 0.0d0
      tempsumham        = 0.0d0
      tounit               = 27.211d0   ! energy conversion Ha to eV
      nstruct              = 0
!!
!! FIXME: This does not yet for for pair NN with only electrostatics (no short range NN)
      if(nntb_flag(3))then
        write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
          'function_hextoff.',hextoff_training_triplet(1),'.'&
           ,hextoff_training_triplet(2),'.'&
           ,hextoff_training_triplet(3),'.data'
        open(symhextoffunit,file=filenametemp,form='formatted',status='old')
        rewind(symhextoffunit)
      endif
!!
 11   continue
      read(symhextoffunit,*,END=10)zdummy
      read(symhextoffunit,*,END=10)(ham(i1),i1=1,matrixsize)
      do i1 = 1,matrixsize
        hammin(i1)=min(hammin(i1),ham(i1))
        hammax(i1)=max(hammax(i1),ham(i1))
        hamav(i1)   =hamav(i1)+ham(i1)
      enddo
      
      nstruct=nstruct+1
      goto 11
 10   continue
!!
      close(symhextoffunit)
      if(nstruct.eq.0) then
        write(*,*) 'Error: ', filenametemp, ' is empty'
        stop
      endif
!!
!! get average energies 
       do i1=1,matrixsize
         hamav(i1)  =hamav(i1)/dble(nstruct)
       enddo
!!
      open(symhextoffunit,file=filenametemp,form='formatted',status='old')
      rewind(symhextoffunit)
!!
 21   continue
      read(symhextoffunit,*,END=20)zdummy
      read(symhextoffunit,*,END=20)(ham(i1),i1=1,matrixsize)
      do i1=1,matrixsize
        tempsumham(i1)=tempsumham(i1)+(ham(i1)-hamav(i1))**2.d0
      enddo
      goto 21
 20   continue
!!
      close(symhextoffunit)
!!
!! write the results
      do i1=1,matrixsize
        tempsumham(i1)=tempsumham(i1)/dble(nstruct)
        hamstddev(i1)=dsqrt(tempsumham(i1))
      enddo
!!
!! write results
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Hamiltonian Element in training set (Ha) (Need to check the units):'
      write(ounit,'(a)')'                   Emin          Emax          average        stddev          range'  
      do i1=1,matrixsize
        write(ounit,'(x,a7,2x,5f15.6)')'Hextoff',hammin(i1),hammax(i1),hamav(i1),hamstddev(i1),abs(hammax(i1)-hammin(i1))
      enddo
      
      return
      end
