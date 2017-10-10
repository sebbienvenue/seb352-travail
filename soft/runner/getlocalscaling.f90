!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Subroutine for computation of local scaling matrix B used in variant 2
!! of the ED-GEKF force update.
!! The matrix B is a 3x3 matrix computed according to:
!! B = [ lambda*I + \sum^Nelem_z (dF/dw_z)^T P_z dF/dw_z ]^-1

!! called by:
!! - optimize_short_combined.f90
!!
      subroutine getlocalscaling(num_weights_short_atomic_free,&
                dfshortdw_xyz,corrmatrixf_list,maxcorrdim,corrdim,&
                num_atoms_element,num_atoms,local_scaling,scalef)
!!
      use fileunits
      use fittingoptions !! lupdatebyelement, nucelem, elemupdate
      use globaloptions  !! maxnum_weights_short_atomic
!!
      implicit none
!!
      integer i2                                                  ! internal
      integer i3                                                  ! internal
      integer i4                                                  ! internal
      integer corrdim_counter                                     ! internal
      integer maxcorrdim                                          ! in
      integer num_atoms                                           ! in
      integer num_atoms_element(nelem)                            ! in
      integer num_weights_short_atomic_free(nelem)                ! in
      integer corrdim(nelem)                                      ! in
!!
      real*8 local_scaling(3,3)                                   ! in/out
      real*8 scalef                                               ! in
      real*8 corrmatrixf_list(maxcorrdim,nelem)                    ! in
      real*8 dfshortdw_xyz(maxnum_weights_short_atomic,3,nelem)   ! in
      real*8 ddot                                                 ! internal
      real*8, dimension(:,:), allocatable :: dfshortdw_xyz_temp   ! internal
      real*8, dimension(:,:), allocatable :: corrmatrixf_dfshortdw ! internal
      real*8, dimension(:,:), allocatable :: corrmatrixf_unpacked  ! internal
!!
      local_scaling(:,:) = 0d0
!!
      do i2=1,nelem
!!
!        if(lupdatebyelement.and.(nucelem(i2).ne.elemupdate))then
!          continue
!        else
          if(num_atoms_element(i2).gt.0)then

!!
!! Allocate temporary arrays for unpacked covariance matrix and reduced
!! dfshortdw_xyz:
            allocate( corrmatrixf_unpacked(num_weights_short_atomic_free(i2),num_weights_short_atomic_free(i2)) )
            allocate( dfshortdw_xyz_temp(num_weights_short_atomic_free(i2),3) )
            allocate( corrmatrixf_dfshortdw(num_weights_short_atomic_free(i2),3) )
!! 
!! Unpack the covariance matrix, which is stored as a columnwise lower triangular matrix
            corrmatrixf_unpacked(:,:) = 0d0
            corrdim_counter = 1
            do i3=1,num_weights_short_atomic_free(i2)
              do i4=i3,num_weights_short_atomic_free(i2)
                corrmatrixf_unpacked(i4,i3) = corrmatrixf_list(corrdim_counter,i2)
                corrdim_counter = corrdim_counter + 1
              enddo ! i4
            enddo ! i3
!!
!! Reduce dfshortdw_xyz to contain only the free weights
            dfshortdw_xyz_temp(:,:) = 0d0
            do i3=1,num_weights_short_atomic_free(i2)
              do i4=1,3
                dfshortdw_xyz_temp(i3,i4) = dfshortdw_xyz(i3,i4,i2)
              enddo ! i4
            enddo ! i3
!! 
!! Since force derivative is scaled by 1/Natoms and scalef, descale:
            dfshortdw_xyz_temp(:,:) = dfshortdw_xyz_temp(:,:) * dble(num_atoms) / scalef
!!
!! Use DSYMM subroutine to compute matrix product P*dF/dw, where P is symmetric.
!! DSYMM performs the matrix matrix operation
!!
!!   C = alpha*A*B + beta*C,
!!
!! where A is a symmetric m by m matrix and B and C are m by n matrices and
!! alpha and beta are scalars
!!
!! call DSYMM(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
!!
!!  SIDE:  'L', symmetric matrix is multiplied from the left
!!  UPLO:  'L', lower triangular matrix of symmetric matrix is used
!!  M:     rows in matrix c, in our case n_weights
!!  N:     columns in matrix c, here 3 (xyz)
!!  ALPHA: 1d0, scalar factor 
!!  A:     m by m symmetric matrix, here unpacked covariance matrix
!!  LDA:   rows of A
!!  B:     m by n matrix, here dF_xyz/dw (n_weights x 3)
!!  LDB:   rows of B
!!  BETA:  0d0, scalar factor
!!  C:     final m by n matrix, here corrmatrixf_dfshortdw (n_weights x 3)
!!  LDC:   rows of C
!!
            corrmatrixf_dfshortdw(:,:) = 0d0
            call dsymm('L','L',num_weights_short_atomic_free(i2),3,1d0,&
                   corrmatrixf_unpacked,num_weights_short_atomic_free(i2),&
                   dfshortdw_xyz_temp,num_weights_short_atomic_free(i2),0d0,&
                   corrmatrixf_dfshortdw,num_weights_short_atomic_free(i2))
!!
!! Use DGEMM subroutine to compute matrix product (dF/dw)^T * corrmatrixf_dfshortdw
!! DGEMM performs the general matrix matrix operation
!!
!!   C = alpha*op(A)*op(B) + beta*C,
!!
!! CALL DGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
!!
!!  TRANSA: Operation performed on A, here transpose 'T' is computed
!!  TRANSB: Operation on B, here no operation 'N' is performed
!!  M:      rows of matrix op(A) and C, here A^T with 3 rows (xyz)
!!  N:      columns of matrix op(B) and C, here 3
!!  K:      columns of op(A) and rows of op(B), here n_weights
!!  ALPHA:  scalar factor, 1d0
!!  A:      first matrix, in this case dFxyz/dw
!!  LDA:    leading dimension of untransformed A, here n_weights
!!  B:      second matrix, in this case corrmatrixf_dfshortdw
!!  LBD:    n_weights
!!  BETA:   scalar factor for addition, since we want to compute a sum of
!!          elemental contributions, set to 1d0
!!  C:      final matrix, local_scaling
!!  LDC:    3
!!
            call dgemm('T','N',3,3,num_weights_short_atomic_free(i2),1d0,&
                   dfshortdw_xyz_temp,num_weights_short_atomic_free(i2),corrmatrixf_dfshortdw,&
                   num_weights_short_atomic_free(i2),1d0,local_scaling,3)
!!
            deallocate(corrmatrixf_unpacked,dfshortdw_xyz_temp,corrmatrixf_dfshortdw)
!! 
          endif ! num_atoms_element(i2).gt.0
!        endif ! lupdatebyelement.and.(nucelem(i2).ne.elemupdate)
      enddo ! i2
!!
!! Compute lambda*I + local_scaling
      local_scaling(1,1) = local_scaling(1,1) + kalmanlambda(1)
      local_scaling(2,2) = local_scaling(2,2) + kalmanlambda(1)
      local_scaling(3,3) = local_scaling(3,3) + kalmanlambda(1)
!!
!! Invert local_scaling
      call invert3x3( local_scaling )
!!
      return
      end
