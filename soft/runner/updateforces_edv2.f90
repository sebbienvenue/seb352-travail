!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Update of the weights using 2nd variant of ED force update.
!!
!!  w = w + P*dF/dw*B*errorf,
!!
!! where P is the covariance matrix (wxw), dF/dw derivative of forces with
!! respect to weights (wx3), B is a local scaling matrix (3x3) and errorf is the
!! errorvector of the forces (3x1)
!!
!! called by:
!! - optimize_short_combined.f90
!!
      subroutine updateforces_edv2(corrdim,num_weights,num_atoms,&
                local_scaling,weights,dfdw_xyz,corrmatrix,errorf_xyz)
!!
      use fileunits
!!
      implicit none
!!
      integer num_weights                                         ! in
      integer num_atoms                                           ! in
      integer i1                                                  ! internal
      integer corrdim                                             ! in
!! 
      real*8 weights(num_weights)                                 ! in/out
      real*8 dfdw_xyz(num_weights,3)                              ! in
      real*8 corrmatrix(corrdim)                                  ! in
      real*8 errorf_xyz(3)                                        ! in
      real*8 local_scaling(3,3)                                   ! in

      real*8 scaled_errorf(3)                                     ! internal
      real*8 scaled_dfdw(num_weights)                             ! internal
      real*8 wupdate(num_weights)                                 ! internal
!! Since errorf_xyz contains 1/Natoms and dF/dw 1/Natoms*scalef, we have to
!! rescale errorvector of forces here to obtain desired scaling of
!! scalef/(3*Natoms)
      errorf_xyz(:) = errorf_xyz(:) * dble(num_atoms) / 3d0
!! 
!! 1) First compute B*errorf, which gives vector of dimension 3
      scaled_errorf(:) = 0d0
!!
      do i1=1,3
        scaled_errorf(i1) = local_scaling(i1,1)*errorf_xyz(1) +&
                            local_scaling(i1,2)*errorf_xyz(2) +&
                            local_scaling(i1,3)*errorf_xyz(3)
      enddo !i1
!! 
!! 2) Then compute matrix vector product dF/dw*scaled_errorf to obtain wx1 matrix 
      scaled_dfdw(:) = 0d0
!!
      do i1=1,num_weights
        scaled_dfdw(i1) = dfdw_xyz(i1,1)*scaled_errorf(1) +& 
                          dfdw_xyz(i1,2)*scaled_errorf(2) +& 
                          dfdw_xyz(i1,3)*scaled_errorf(3) 
      enddo !i1
!!
!! 3) Then compute product symmetric packed matrix corrmatrix times scaled_dfdw
!!    to finally obtain update vector wupdate (wx1) using DSPMV (see e.g.
!!    updatekalman.f90 for detailed description of subroutine
      wupdate(:)=0.0d0
      call dspmv('l',num_weights,1.d0,corrmatrix,scaled_dfdw,1,0.d0,wupdate,1)
!!      WRITE(*,*) DOT_PRODUCT(wupdate,wupdate), "V2"
!!
!! 4) Finally, update the weights
      do i1=1,num_weights
        weights(i1) = weights(i1) + wupdate(i1)
      enddo
!!
      return
      end
