!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! MG: Subroutine for inversion of 3x3 matrix used in computation of
!! MG: local_scaling matrix.
!! MG: Inverse is computed as M^-1 = 1 / det(M) * (cofactor matrix of M)^T
!! MG: Transpose is obtained by swapping indices

!! called by:
!! - getlocalscaling.f90
!!
    subroutine invert3x3(M)
!!
      use fileunits
!!
      implicit none
!!
      real*8 M(3,3)                                               ! in/out
      real*8 M_inv(3,3)                                           ! internal
      real*8 M_det                                                ! internal
      real*8 inv_limit                                            ! internal
!! Set threshold for detection of zero determinants
      inv_limit = 1.0d-9
!!
      M_inv(:,:) = 0d0
!! Compute cofactor matrices
      M_inv(1,1) = M(2,2)*M(3,3) - M(3,2)*M(2,3)
      M_inv(1,2) = M(3,2)*M(1,3) - M(1,2)*M(3,3)
      M_inv(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
!! Compute determinant of matrix
      M_det = M_inv(1,1)*M(1,1) + M_inv(1,2)*M(2,1) + M_inv(1,3)*M(3,1)
!!
!! Check if matrix is invertible:
      if (abs(M_det).le.inv_limit) then
        write(ounit,*) M
        write(ounit,*) M_det
        write(ounit,*) "Error: 3x3 Matrix not invertible"
        stop  
      endif
!!
!! Compute first row of inverse matrix
      M_inv(1,1) = M_inv(1,1) / M_det
      M_inv(1,2) = M_inv(1,2) / M_det
      M_inv(1,3) = M_inv(1,3) / M_det
!! Compute second row
      M_inv(2,1) = ( M(2,3)*M(3,1) - M(2,1)*M(3,3) ) / M_det
      M_inv(2,2) = ( M(1,1)*M(3,3) - M(3,1)*M(1,3) ) / M_det
      M_inv(2,3) = ( M(2,1)*M(1,3) - M(1,1)*M(2,3) ) / M_det
!! Compute third row
      M_inv(3,1) = ( M(2,1)*M(3,2) - M(2,2)*M(3,1) ) / M_det
      M_inv(3,2) = ( M(3,1)*M(1,2) - M(1,1)*M(3,2) ) / M_det
      M_inv(3,3) = ( M(1,1)*M(2,2) - M(1,2)*M(2,1) ) / M_det
!!
!! Overwrite old matrix
      M(:,:) = M_inv(:,:)
!!
      return
!!
    end subroutine
