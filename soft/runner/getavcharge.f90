!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fitting_batch.f90
!! - fittingpair.f90

      subroutine getavcharge(ntrain,&
          avcharge,stddevcharge,&
          chargemin,chargemax)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none

      integer ntrain
      integer idummy
      integer i0,i1,i2
      integer iatom
      integer zelem
      integer num_atoms_element(nelem)
      integer unit_local

      real*8 avcharge(nelem)              ! out 
      real*8 zdummy
      real*8 charge
      real*8 stddevcharge(nelem)          ! out
      real*8 tempsum(nelem)
      real*8 chargemin(nelem)                                        ! in/out
      real*8 chargemax(nelem)                                        ! in/out

      logical lperiodic

      num_atoms_element(:)=0
      avcharge(:)=0.0d0
      tempsum(:)=0.0d0
      stddevcharge(:)=0.0d0

      if(nn_type_elec.eq.1)then
        unit_local=symeunit ! functione.data
        open(unit_local,file='functione.data',form='formatted',status='old')
      elseif(nn_type_elec.eq.2)then
        unit_local=symunit ! function.data
        open(unit_local,file='function.data',form='formatted',status='old')
      endif
      rewind(unit_local)

      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      rewind(trainstructunit)
!!
      write(ounit,*)'-------------------------------------------------------------'
!!'
      do i0=1,ntrain
!! get the number of atoms 
        read(unit_local,*)iatom 
        do i1=1,iatom
          read(unit_local,*)idummy
        enddo ! i1
        read(unit_local,*)zdummy
!! get the charges
        read(trainstructunit,*)idummy,lperiodic
        if(lperiodic)then
          do i2=1,3
            read(trainstructunit,*)zdummy
          enddo
        endif
        do i1=1,iatom
          read(trainstructunit,*)zelem,zdummy,zdummy,zdummy,charge
          avcharge(elementindex(zelem))&
            =avcharge(elementindex(zelem))+charge 
          num_atoms_element(elementindex(zelem))&
            =num_atoms_element(elementindex(zelem))+1
          chargemax(elementindex(zelem))=max(chargemax(elementindex(zelem)),charge)
          chargemin(elementindex(zelem))=min(chargemin(elementindex(zelem)),charge)
        enddo
      enddo ! i0

      close(unit_local)
      close(trainstructunit)
!!
      do i1=1,nelem
!!        write(ounit,*)'-------------------------------------------------------------'
        avcharge(i1)=avcharge(i1)/dble(num_atoms_element(i1))
      enddo


!! get the standard deviation
      if(nn_type_elec.eq.1)then
        unit_local=symeunit ! functione.data
        open(unit_local,file='functione.data',form='formatted',status='old')
      elseif(nn_type_elec.eq.2)then
        unit_local=symunit ! function.data
        open(unit_local,file='function.data',form='formatted',status='old')
      endif
      rewind(unit_local)

      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      rewind(trainstructunit)

      do i0=1,ntrain
!! get the number of atoms 
        read(unit_local,*)iatom 
        do i1=1,iatom
          read(unit_local,*)idummy
        enddo ! i1
        read(unit_local,*)zdummy
!! get the charges
        read(trainstructunit,*)idummy,lperiodic
        if(lperiodic)then
          do i2=1,3
            read(trainstructunit,*)zdummy
          enddo
        endif
        do i1=1,iatom
          read(trainstructunit,*)zelem,zdummy,zdummy,zdummy,charge
          tempsum(elementindex(zelem))=&
            tempsum(elementindex(zelem))+(charge-avcharge(elementindex(zelem)))**2.d0
        enddo
      enddo ! i0

      close(unit_local)
      close(trainstructunit)

      do i1=1,nelem
        tempsum(i1)=tempsum(i1)/dble(num_atoms_element(i1))
        stddevcharge(i1)=dsqrt(tempsum(i1))
      enddo
!!
!! better output:
      write(ounit,*)'Charges in training set (e):'
      write(ounit,*)'              Qmin          Qmax        average        stddev     range'
      do i1=1,nelem
        write(ounit,'(x,a2,4x,5f14.6)')element(i1),chargemin(i1),chargemax(i1),avcharge(i1),&
          stddevcharge(i1),abs(chargemax(i1)-chargemin(i1)) 
      enddo
!!
      return
      end
