!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - ewaldrecip.f90
!!
      subroutine getreclat(lattice,reclattice,determinant)
!!
      implicit none
!!
      integer i,j
!!
      real*8 lattice(3,3)
      real*8 reclattice(3,3)
      real*8 determinant
      real*8 twopi
!!
      parameter (twopi=6.283185307d0)
!!
!! a x b
      reclattice(3,1)=lattice(1,2)*lattice(2,3)&
                     -lattice(1,3)*lattice(2,2)
      reclattice(3,2)=lattice(1,3)*lattice(2,1)&
                     -lattice(1,1)*lattice(2,3)
      reclattice(3,3)=lattice(1,1)*lattice(2,2)&
                     -lattice(1,2)*lattice(2,1)
!!
!! b x c
      reclattice(1,1)=lattice(2,2)*lattice(3,3)&
                     -lattice(2,3)*lattice(3,2)
      reclattice(1,2)=lattice(2,3)*lattice(3,1)&
                     -lattice(2,1)*lattice(3,3)
      reclattice(1,3)=lattice(2,1)*lattice(3,2)&
                     -lattice(2,2)*lattice(3,1)
!!
!! c x a
      reclattice(2,1)=lattice(3,2)*lattice(1,3)&
                     -lattice(3,3)*lattice(1,2)
      reclattice(2,2)=lattice(3,3)*lattice(1,1)&
                     -lattice(3,1)*lattice(1,3)
      reclattice(2,3)=lattice(3,1)*lattice(1,2)&
                     -lattice(3,2)*lattice(1,1)
!!
!! calculation of the determinant
      determinant=0.0d0
      determinant=lattice(1,1)*reclattice(1,1)&
                 +lattice(1,2)*reclattice(1,2)&
                 +lattice(1,3)*reclattice(1,3)
!!
!!      write(*,*)'determinant ',determinant
!!
      do i=1,3
        do j=1,3
!!          reclattice(i,j)=reclattice(i,j)/determinant
          reclattice(i,j)=reclattice(i,j)*twopi/determinant
        enddo
      enddo
!!
!! debug
!!      write(*,*)'reciprocal lattice'
!!      do i=1,3
!!        write(*,'(3f14.6)')(reclattice(i,j),j=1,3)
!!      enddo
!!
      return
      end
