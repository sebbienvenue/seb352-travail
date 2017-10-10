!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine getforcestatisticspair()
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1
      integer natoms               ! internal
      integer num_pairs               ! internal
      integer ielem                 ! internal
      integer idummy                ! internal
      integer num_atoms_element(nelem) ! internal
!!
      real*8 fmin(nelem)            ! internal
      real*8 fmax(nelem)            ! internal
      real*8 tounit                 ! internal
      real*8 zdummy                 ! internal
      real*8 fx                     ! internal
      real*8 fy                     ! internal
      real*8 fz                     ! internal
      real*8 f                      ! internal
      real*8 avforce(nelem)         ! internal
      real*8 stddevforce(nelem)     ! internal 
      real*8 tempsum(nelem)         ! internal
!!
      logical lperiodic             ! internal
!!
!! initializations
      fmin(:)              =1000.d0
      fmax(:)              =0.0d0
      avforce(:)           =0.0d0
      num_atoms_element(:) =0
      stddevforce(:)       =0.0d0
      tempsum(:)           =0.0d0
      tounit               = 27.211d0   ! energy conversion Ha to eV
!!
      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(trainstructunit)
      rewind(symunit)
!!
 11   continue
      read(symunit,*,END=10)natoms,num_pairs
      do i1=1,num_pairs+1
        read(symunit,*)idummy
      enddo
      read(trainstructunit,*)idummy,lperiodic
      if(lperiodic)then
        do i1=1,3
          read(trainstructunit,*)zdummy
        enddo
      endif
      do i1=1,natoms
        read(trainstructunit,*)ielem,zdummy,zdummy,zdummy,zdummy,zdummy,fx,fy,fz
        f=fx**2+fy**2+fz**2
        f=dsqrt(f)
        fmin(elementindex(ielem))=min(fmin(elementindex(ielem)),abs(f))
        fmax(elementindex(ielem))=max(fmax(elementindex(ielem)),abs(f))
        num_atoms_element(elementindex(ielem))=num_atoms_element(elementindex(ielem))+1
        avforce(elementindex(ielem))=avforce(elementindex(ielem))+f
      enddo
      goto 11
 10   continue
!!
      close(trainstructunit)
      close(symunit)

!!
!! get average forces
      do i1=1,nelem
        avforce(i1)=avforce(i1)/dble(num_atoms_element(i1))
      enddo

!!
      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(trainstructunit)
      rewind(symunit)
!!
 21   continue
      read(symunit,*,END=20)natoms,num_pairs
      do i1=1,num_pairs+1
        read(symunit,*)idummy
      enddo
      read(trainstructunit,*)idummy,lperiodic
      if(lperiodic)then
        do i1=1,3
          read(trainstructunit,*)zdummy
        enddo
      endif
      do i1=1,natoms
        read(trainstructunit,*)ielem,zdummy,zdummy,zdummy,zdummy,zdummy,fx,fy,fz
        f=fx**2+fy**2+fz**2
        f=dsqrt(f)
        tempsum(elementindex(ielem))=tempsum(elementindex(ielem))+(f-avforce(elementindex(ielem)))**2.d0
      enddo
      goto 21
 20   continue
!!
      close(trainstructunit)
      close(symunit)

!!
!! write the results
      do i1=1,nelem
        tempsum(i1)=tempsum(i1)/dble(num_atoms_element(i1))
        stddevforce(i1)=dsqrt(tempsum(i1))
      enddo
!!
!! write results
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Forces in training set (Ha/Bohr):'
      write(ounit,'(a)')'               Fmin          Fmax          average        stddev'  
      do i1=1,nelem
        write(ounit,'(x,a2,2x,4f15.6)')element(i1), fmin(i1),  fmax(i1),  avforce(i1),  stddevforce(i1)
      enddo
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Forces in training set (eV/Bohr):'
      write(ounit,'(a)')'               Fmin          Fmax          average        stddev'  
      do i1=1,nelem
        write(ounit,'(x,a2,2x,4f15.6)')element(i1),  fmin(i1)*tounit,  fmax(i1)*tounit,  avforce(i1)*tounit,  stddevforce(i1)*tounit
      enddo

!!
      return
      end
