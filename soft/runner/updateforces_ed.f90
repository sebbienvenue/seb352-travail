!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - optimize_short_combined.f90
!!
      subroutine updateforces_ed(corrdim,num_weights,num_atoms,&
                global_scaling,weights,dfdw,corrmatrix,errorf,scalef)
!!
      use fileunits
!! Don't use fittingoptions, because updatekalman is used for different quantities
!!
      implicit none
!!
      integer num_weights                                         ! in
      integer num_atoms                                           ! in
      integer i1                                                  ! internal
      integer corrdim                                             ! in
!! 
      real*8 weights(num_weights)                                 ! in/out
      real*8 dfdw(num_weights)                                    ! in
      real*8 corrmatrix(corrdim)                                  ! in
      real*8 errorf                                               ! in
      real*8 scalef                                               ! in
      real*8 global_scaling                                       ! in
      real*8 coh(num_weights)                                     ! internal
      real*8 ddot                                                 ! internal
!!
!! initialization
      coh(:)=0.0d0
!! to get 1/3N scaling, divide errorf by 3, and since both dF/dw and errorf are
!! multiplied by scaling factor, remove it once.
!! In addition, since dF/dw is also scaled by 1/N and we do not want 1/N^2
!! behaviour, multiply by N once.
!!D      WRITE(*,*) errorf, "I"
!!      errorf = errorf / 3d0 / scalef * dble( num_atoms )
!! calculate coh=corrmatrix*dfdw
      call dspmv('l',num_weights,1.d0,corrmatrix,dfdw,1,0.d0,coh,1)
!!
!! result: coh = P(n-1)*J(n)
!!
!!
      coh(:) = coh(:) / 3d0 / scalef * dble( num_atoms )
!!
!! Update the weights
      do i1=1,num_weights
        weights(i1)=weights(i1)+errorf*coh(i1)*global_scaling
      enddo
!!D      WRITE(*,*) DOT_PRODUCT((errorf*coh(:)*global_scaling),(errorf*coh(:)*global_scaling)), "V1"
!!D      WRITE(*,*) errorf, "II"
!!
      return
      end
