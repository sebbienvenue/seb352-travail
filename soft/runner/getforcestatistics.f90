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
      subroutine getforcestatistics(belowfmax,fmin,fmax,fvecmin,fvecmax)
!!
      use fileunits
      use nnflags
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer i1
      integer num_atoms             ! internal
      integer num_pairs             ! internal
      integer ndim                  ! internal
      integer ielem                 ! internal
      integer idummy                ! internal
      integer num_atoms_element(nelem) ! internal
      integer belowfmax(nelem)      ! out
!!
      real*8 fmin(nelem)            ! out 
      real*8 fmax(nelem)            ! out 
      real*8 fvecmin(nelem)         ! out 
      real*8 fvecmax(nelem)         ! out 
      real*8 tounit                 ! internal
      real*8 zdummy                 ! internal
      real*8 fx                     ! internal
      real*8 fy                     ! internal
      real*8 fz                     ! internal
      real*8 f                      ! internal
      real*8 avforce(nelem)         ! internal
      real*8 stddevforce(nelem)     ! internal 
      real*8 tempsum(nelem)         ! internal
      real*8 edummy                 ! internal
!!
      logical lperiodic             ! internal
!!
!! initializations
      fmin(:)              = 1000.d0
      fmax(:)              = 0.0d0
      fvecmin(:)           = 1000.d0
      fvecmax(:)           = 0.0d0
      avforce(:)           = 0.0d0
      num_atoms_element(:) = 0
      stddevforce(:)       = 0.0d0
      tempsum(:)           = 0.0d0
      tounit               = 27.211d0   ! energy conversion Ha to eV
      belowfmax(:)         = 0
!!
      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(trainstructunit)
      rewind(symunit)
!!
 11   continue
      if(nn_type_short.eq.1)then
        read(symunit,*,END=10)num_atoms
        ndim=num_atoms
      elseif(nn_type_short.eq.2)then
        read(symunit,*,END=10)num_atoms,num_pairs
        ndim=num_pairs
      else
        write(ounit,*)'ERROR in getforcestatistics, unknown nn_type_short '&
          ,nn_type_short
        stop !'
      endif
      do i1=1,ndim
        read(symunit,*)idummy
      enddo
      read(symunit,*)edummy
      read(trainstructunit,*)idummy,lperiodic
      if(lperiodic)then
        do i1=1,3
          read(trainstructunit,*)zdummy
        enddo
      endif
      do i1=1,num_atoms
        read(trainstructunit,*)ielem,zdummy,zdummy,zdummy,zdummy,zdummy,fx,fy,fz
        f=fx**2+fy**2+fz**2
        f=dsqrt(f)
        fvecmin(elementindex(ielem))=min(fmin(elementindex(ielem)),abs(f))
        fvecmax(elementindex(ielem))=max(fmax(elementindex(ielem)),abs(f))
        fmin(elementindex(ielem))=min(fmin(elementindex(ielem)),abs(fx))
        fmax(elementindex(ielem))=max(fmax(elementindex(ielem)),abs(fx))
        fmin(elementindex(ielem))=min(fmin(elementindex(ielem)),abs(fy))
        fmax(elementindex(ielem))=max(fmax(elementindex(ielem)),abs(fy))
        fmin(elementindex(ielem))=min(fmin(elementindex(ielem)),abs(fz))
        fmax(elementindex(ielem))=max(fmax(elementindex(ielem)),abs(fz))
        num_atoms_element(elementindex(ielem))=num_atoms_element(elementindex(ielem))+1
        avforce(elementindex(ielem))=avforce(elementindex(ielem))+f
        if(abs(fx).le.maxforce)belowfmax(elementindex(ielem))=belowfmax(elementindex(ielem))+1
        if(abs(fy).le.maxforce)belowfmax(elementindex(ielem))=belowfmax(elementindex(ielem))+1
        if(abs(fz).le.maxforce)belowfmax(elementindex(ielem))=belowfmax(elementindex(ielem))+1
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
      if(nn_type_short.eq.1)then
        read(symunit,*,END=20)num_atoms
        ndim=num_atoms
      elseif(nn_type_short.eq.2)then
        read(symunit,*,END=20)num_atoms,num_pairs
        ndim=num_pairs
      else
        write(ounit,*)'ERROR: unknown nn_type_short in getforcestatistics ',nn_type_short
        stop
      endif
      do i1=1,ndim+1

!!
!! CHANGE ANDI: GFORTRAN: gfortran crashes here because in line number ndim+1 there is a not an integer but a real,
!!                        to skip lines one can also completely omit the argument of the read statement.
!!
       !read(symunit,*)idummy
        read(symunit,*)
!! END CHANGE

      enddo
      read(trainstructunit,*)idummy,lperiodic
      if(lperiodic)then
        do i1=1,3
          read(trainstructunit,*)zdummy
        enddo
      endif
      do i1=1,num_atoms
        read(trainstructunit,*)ielem,zdummy,zdummy,zdummy,zdummy,zdummy,fx,fy,fz
        f=fx**2+fy**2+fz**2
        f=dsqrt(f)
        tempsum(elementindex(ielem))=&
          tempsum(elementindex(ielem))+(f-avforce(elementindex(ielem)))**2.d0
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
      write(ounit,*)'Force vectors in training set (Ha/Bohr):'
      write(ounit,'(a)')'               Fmin          Fmax          average        stddev          range'  
      do i1=1,nelem
        write(ounit,'(x,a2,2x,5f15.6)')&
          element(i1),fvecmin(i1),fvecmax(i1),avforce(i1),&
          stddevforce(i1),abs(fvecmax(i1)-fvecmin(i1))
      enddo
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Force vectors in training set (eV/Bohr):'
      write(ounit,'(a)')'               Fmin          Fmax          average        stddev          range'  
      do i1=1,nelem
        write(ounit,'(x,a2,2x,5f15.6)')&
          element(i1),fvecmin(i1)*tounit,fvecmax(i1)*tounit,&
          avforce(i1)*tounit,stddevforce(i1)*tounit,&
          abs(fvecmax(i1)-fvecmin(i1))*tounit
      enddo
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Force components in training set (Ha/Bohr):'
      write(ounit,'(a)')'               Fmin          Fmax          range'
      do i1=1,nelem
        write(ounit,'(x,a2,2x,3f15.6)')&
          element(i1),fmin(i1),fmax(i1),abs(fmax(i1)-fmin(i1))
      enddo
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Force components in training set (eV/Bohr):'
      write(ounit,'(a)')'               Fmin          Fmax          range'
      do i1=1,nelem
        write(ounit,'(x,a2,2x,3f15.6)')&
          element(i1),fmin(i1)*tounit,fmax(i1)*tounit,&
          abs(fmax(i1)-fmin(i1))*tounit
      enddo
!!
      return
      end
