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
      subroutine getenergystatistics(iswitch,belowmaxenergy,&
        eshortav,eshortstddev,eshortmin,eshortmax)
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
      integer natoms                ! internal
      integer ndim                  ! internal
      integer nsizepair             ! internal
      integer idummy                ! internal
      integer nstruct               ! internal
      integer belowmaxenergy        ! out
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
      eewaldmax            = -100000000.0d0
      eewaldav             = 0.0d0
      eewaldstddev         = 0.0d0
      tempsumeshort        = 0.0d0
      tempsumeewald        = 0.0d0
      tempsumetot          = 0.0d0
      tounit               = 27.211d0   ! energy conversion Ha to eV
      nstruct              = 0
      belowmaxenergy       = 0
!!
!! FIXME: This does not yet for for pair NN with only electrostatics (no short range NN)
      if(lshort)then
        open(symunit,file='function.data',form='formatted',status='old')
      elseif(lelec.and.(nn_type_elec.eq.1))then ! make sure it works also for lshort=F
        open(symunit,file='functione.data',form='formatted',status='old')
      endif
      rewind(symunit)
!!
 11   continue
      if(iswitch.eq.1)then
        read(symunit,*,END=10)natoms
        ndim=natoms
      elseif((iswitch.eq.2))then
        read(symunit,*,END=10)natoms,nsizepair
        ndim=nsizepair
      else
        write(ounit,*)'ERROR: unknown iswitch in getenergystatistics ',iswitch
        stop
      endif
      nstruct=nstruct+1
      do i1=1,ndim
        read(symunit,*)idummy
      enddo
      read(symunit,*)zdummy,etot,eshort,eewald
      if(eshort.le.maxenergy)then
        belowmaxenergy=belowmaxenergy+1
      endif

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
      if(iswitch.eq.1)then
        read(symunit,*,END=20)natoms
        ndim=natoms
      elseif(iswitch.eq.2)then
        read(symunit,*,END=20)natoms,nsizepair
        ndim=nsizepair
      else
        write(ounit,*)'ERROR: unknown iswitch in getenergystatistics ',iswitch
        stop
      endif
      do i1=1,ndim
        read(symunit,*)idummy
      enddo
      read(symunit,*)zdummy,etot,eshort,eewald

      tempsumeshort=tempsumeshort+(eshort-eshortav)**2.d0
      tempsumeewald=tempsumeewald+(eewald-eewaldav)**2.d0
      tempsumetot  =tempsumetot  +(etot-etotav)**2.d0
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
      write(ounit,'(a)')'                   Emin          Emax          average        stddev          range'  
        write(ounit,'(x,a7,2x,5f15.6)')'Eshort ',eshortmin,eshortmax,eshortav,eshortstddev,abs(eshortmax-eshortmin)
!! FIXME: is the following Ewald line really meaningful?
        write(ounit,'(x,a7,2x,5f15.6)')'Eelec  ',eewaldmin,eewaldmax,eewaldav,eewaldstddev,abs(eewaldmax-eewaldmin)
        write(ounit,'(x,a7,2x,5f15.6)')'Etot   ',etotmin,etotmax,etotav,etotstddev,abs(etotmax-etotmin)
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Energies in training set (eV/atom):'
      write(ounit,'(a)')'                   Emin          Emax          average        stddev          range'  
        write(ounit,'(x,a7,2x,5f15.6)')'Eshort ',eshortmin*tounit,eshortmax*tounit,eshortav*tounit,&
          eshortstddev*tounit,abs(eshortmax-eshortmin)*tounit
!! FIXME: is the following Ewald line really meaningful?
        write(ounit,'(x,a7,2x,5f15.6)')'Eelec  ',eewaldmin*tounit,eewaldmax*tounit,eewaldav*tounit,&
          eewaldstddev*tounit,abs(eewaldmax-eewaldmin)*tounit
        write(ounit,'(x,a7,2x,5f15.6)')'Etot   ',etotmin*tounit,etotmax*tounit,etotav*tounit,&
          etotstddev*tounit,abs(etotmax-etotmin)*tounit
!!
      return
      end
