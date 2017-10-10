!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! called by: 
!! - fitting.f90
!! - fitting_batch.f90
!! - fittingpair.f90
!! - geterror.f90
!! - geterrorpair.f90
!!
      subroutine readforces(unit,npoints,force_list_local)
!!
      use globaloptions
      use structures
!! don't use module fileunits here, because unit can have several values
!!
      implicit none
!!
      integer unit                                               ! in
      integer npoints                                            ! in
      integer dummy                                              ! internal
      integer i1,i2,i3                                           ! internal
!!
      real*8 force_list_local(3,max_num_atoms,nblock)
!!       
      force_list_local(:,:,:)=0.0d0
!!
      do i1=1,npoints
        read(unit,*)dummy
        do i2=1,num_atoms_list(i1)
          read(unit,*)(force_list_local(i3,i2,i1),i3=1,3)
!          write(*,*) force_list_local(:,i2,i1)
        enddo
      enddo ! i1 
!!
      return
      end
