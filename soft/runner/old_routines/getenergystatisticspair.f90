!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine getenergystatisticspair(nelem,nucelem,&
        eshortav,eshortstddev,eshortmin,eshortmax,&
        element,lshort,ldebug)
!!
      use fileunits
!!
      implicit none
!!
      integer i1
      integer nelem                 ! in
      integer nucelem(nelem)        ! in
      integer natoms                ! internal
      integer nsizepair             ! internal
      integer idummy                ! internal
      integer nstruct               ! internal
!!
      real*8 eshortmin              ! out
      real*8 eshortmax              ! out
      real*8 eshortav               ! out
      real*8 eshortstddev           ! out 
      real*8 etotmin                ! internal
      real*8 etotmax                ! internal
      real*8 etotav                 ! internal
      real*8 etotstddev             ! internal
      real*8 eewaldmin              ! internal
      real*8 eewaldmax              ! internal
      real*8 eewaldav               ! internal
      real*8 eewaldstddev           ! internal
      real*8 etot                   ! internal
      real*8 eshort                 ! internal
      real*8 eewald                 ! internal
      real*8 tounit                 ! internal
      real*8 zdummy                 ! internal
      real*8 tempsumeshort          ! internal
      real*8 tempsumeewald          ! internal
      real*8 tempsumetot            ! internal
      real*8 etotdev                ! internal
!!
      character*2 element(nelem)    ! in
!!
      logical lshort                ! in
      logical ldebug                ! in
!!
!! initializations
      etotmin              = 1000000000.d0
      etotmax              = -1000000000.0d0
      etotav               = 0.0d0
      etotdev              = 0.0d0
      eshortmin            = 1000000000.d0
      eshortmax            = -1000000000.0d0 
      eshortav             = 0.0d0
      eshortstddev         = 0.0d0
      eewaldmin            = 1000000000.d0
      eewaldmax            = 0.0d0
      eewaldav             = 0.0d0
      eewaldstddev         = 0.0d0
      tempsumeshort        = 0.0d0
      tempsumeewald        = 0.0d0
      tempsumetot          = 0.0d0
      tounit               = 27.211d0   ! energy conversion Ha to eV
      nstruct              = 0
!!
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)
!!
 11   continue
      read(symunit,*,END=10)natoms,nsizepair
      nstruct=nstruct+1
      do i1=1,nsizepair !natoms
        read(symunit,*)idummy
      enddo
      read(symunit,*)zdummy,etot,eshort,eewald

      etotmin  =min(etotmin,etot)
      eshortmin=min(eshortmin,eshort)
      eewaldmin=min(eewaldmin,eewald)

      etotmax  =max(etotmax,etot)
      eshortmax=max(eshortmax,eshort)
      eewaldmax=max(eewaldmax,eewald)

      etotav   =etotav+etot
      eshortav =eshortav+eshort
      eewaldav =eewaldav+eewald
      goto 11
 10   continue
!!
      close(symunit)
!!
!! get average energies 
      etotav  =etotav/dble(nstruct)
      eshortav=eshortav/dble(nstruct)
      eewaldav=eewaldav/dble(nstruct)
!!
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)
!!
 21   continue
      read(symunit,*,END=20)natoms,nsizepair
      do i1=1,nsizepair 
        read(symunit,*)idummy
      enddo
      read(symunit,*)zdummy,etot,eshort,eewald

      tempsumeshort=tempsumeshort+(eshort-eshortav)**2.d0
      tempsumeewald=tempsumeewald+(eewald-eewaldav)**2.d0
      tempsumetot= tempsumetot+(etot-etotav)**2.d0
      goto 21
 20   continue
!!
      close(symunit)
!!
!! write the results
      tempsumeshort=tempsumeshort/dble(nstruct)
      tempsumeewald=tempsumeewald/dble(nstruct)
      tempsumetot  =tempsumetot/dble(nstruct)

      eshortstddev=dsqrt(tempsumeshort)
      eewaldstddev=dsqrt(tempsumeewald)
      etotstddev  =dsqrt(tempsumetot)
!!
!! write results
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Energies in training set (Ha/atom):'
      write(ounit,'(a)')'                   Emin          Emax          average        stddev'  
        write(ounit,'(x,a7,2x,4f15.6)')'Eshort ',eshortmin,eshortmax,eshortav,eshortstddev
        write(ounit,'(x,a7,2x,4f15.6)')'Eewald ',eewaldmin,eewaldmax,eewaldav,eewaldstddev
        write(ounit,'(x,a7,2x,4f15.6)')'Etot   ',etotmin,etotmax,etotav,etotstddev
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Energies in training set (eV/atom):'
      write(ounit,'(a)')'                   Emin          Emax          average        stddev'  
        write(ounit,'(x,a7,2x,4f15.6)')'Eshort ',eshortmin*tounit,eshortmax*tounit,eshortav*tounit,eshortstddev*tounit
        write(ounit,'(x,a7,2x,4f15.6)')'Eewald ',eewaldmin*tounit,eewaldmax*tounit,eewaldav*tounit,eewaldstddev*tounit
        write(ounit,'(x,a7,2x,4f15.6)')'Etot   ',etotmin*tounit,etotmax*tounit,etotav*tounit,etotstddev*tounit
!!
      return
      end
