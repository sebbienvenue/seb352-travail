!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - geterror.f90
!! - precondition.f90
!! - fitting_batch.f90
!! - fittingpair.f90 
!! - geterrorpair.f90
!! - preconditionpair.f90
!!
      subroutine readfunctions_hextoff(iswitch,unit,npoints,ndim,&
         max_num,maxnum_funcvalues_local,num_funcvalues_local,&
         symfunction_list_local)
!!
      use fileunits
      use globaloptions
      use structures
      use basismod
!!
      implicit none
!!     
      integer unit                                                  ! internal
      integer ndim                                                  ! in
      integer iswitch                                               ! in
      integer npoints                                               ! in
      integer max_num                                               ! in
      integer maxnum_funcvalues_local                                     ! in
      integer num_funcvalues_local(ndim)                                  ! in
      integer i1,i2,i3,i4                                           ! internal
      real*8  dummy,dummy1,dummy2,dummy3,dummy4                     ! internal
!!
      integer matrixsize
      real*8, dimension(:), allocatable :: hextoffin                ! internal
      real*8 symfunction_list_local(maxnum_funcvalues_local,nblock)     ! out
!!
      matrixsize  = num_basis(elementindex(hextoff_training_triplet(1)))*&
                         num_basis(elementindex(hextoff_training_triplet(2)))
      allocate(hextoffin(matrixsize))

      do i1=1,npoints
          num_atoms_list(i1) = 3
          read(unit,*)(symfunction_list_local(i3,i1),i3=1,num_funcvalues_local(1))
          read(unit,*)(hextoffin(i2),i2=1,matrixsize)
          i4 = 0
          do i2=1,num_basis(elementindex(hextoff_training_triplet(1)))
            do i3=1,num_basis(elementindex(hextoff_training_triplet(2)))
              i4 = i4 + 1
              hextoff_list(i1,i2,i3) = hextoffin(i4)
            enddo
          enddo
      enddo ! i1
!!
      deallocate(hextoffin)
      return
      end
