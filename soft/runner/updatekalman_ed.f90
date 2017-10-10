!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! MG: Element decoupled update routine for Kalman filter
!! MG: Using the element decoupled filter, computation of the
!! MG: Kalman gain becomes:
!! MG:    K = P*J*A,
!! MG: with
!! MG     A = { lambda + SUM^nelem[ J^T(elem)*P(elem)*J(elem) ] }^-1
!! MG: Global scaling matrix A is computed externally
!! MG: (This routine is just an adapted copy of the original
!! MG:  updatekalman.f90)

!! called by:
!! - optimize_short_combined.f90
!!
      subroutine updatekalman_ed(kaldim,corrdim,&
                num_weights,kalmanlambda,kalmannue,&
                weights,dedw,corrmatrix,error,global_scaling,num_atoms)
!!
      use fileunits
!!
      implicit none
!!
      integer num_weights                                         ! in
      integer i1                                                  ! internal
      integer kaldim                                              ! in
      integer corrdim                                             ! in
      integer num_atoms                                           ! in
!! 
      real*8 kalmanlambda                                         ! in/out
      real*8 kalmannue                                            ! in
      real*8 weights(num_weights)                                 ! in/out
      real*8 dedw(num_weights)                                    ! in
!! MG: Fix2: temporary array for rescaled gradient:
      real*8 dedw_tmp(num_weights)                                ! in
      real*8 kalgainmat(kaldim)                                   ! internal
      real*8 corrmatrix(corrdim)                                  ! in/out
      real*8 error                                                ! in
!! MG: Fix2: temporary double for rescaled error:
      real*8 error_tmp                                            ! in
      real*8 alpha                                                ! internal
      real*8 invlambda                                            ! internal
      real*8 coh(num_weights)                                     ! internal
!! MG: global scaling matrix
      real*8 global_scaling                                       ! in
!!
!! debugging:
!!      write(*,*)'kalmanlambda ',kalmanlambda
!!      write(*,*)'kalmannue ',kalmannue
!!
!! initialization
      coh(:)=0.0d0
!! MG: Fix2: Rescale gradients and errors
      error_tmp = error*dble(num_atoms)
      dedw_tmp = dedw*dble(num_atoms)
!! Kalman forgetting schedule
      invlambda=1.d0/kalmanlambda
!!
!! 1. Calculation of the Kalman gain matrix
!!-----------------------------------------
!!
!! Calculation of lambda^-1*P(n-1)*J(n)   ! where does lambda enter???
!! output of dspmv is "coh"=corrmatrix*dedw
!! DSPMV  performs the matrix-vector operation
!!
!!     y := alpha*A*x + beta*y,
!!
!!  where alpha and beta are scalars, x and y are n element vectors and
!!  A is an n by n symmetric matrix, supplied in packed form.
!!  SUBROUTINE DSPMV (UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY)
!!
!!  UPLO: 'l' lower triangle of A is stored
!!  N:    order of matrix A                       - num_weights
!!  alpha: scalar, unchanged
!!  AP: array of dimensions n*(n+1)/2 , unchanged - corrmatrix
!!  X:    vector, unchanged                       - dedw
!!  INCX: increment of the elements of X
!!  beta: scalar, unchanged
!!  Y: output                                     - coh = corrmatrix*dedw 
!!
!!      call omp_set_num_threads(8)
!!      call mkl_set_num_threads(8)
!!
!! calculate coh=corrmatrix*dedw
!! MG: Fix2: Use scaled temporary array
      call dspmv('l',num_weights,1.d0,corrmatrix,dedw_tmp,1,0.d0,coh,1)

!! result: coh = P(n-1)*J(n)
!!
!! MG: Use global scaling factor A to compute Kalman gain matrix
!! MG:   K = P(n-1)*J(n)*A

      kalgainmat(:)=global_scaling*coh(:) 

!! debug
!!      do i1=1,num_weights
!!        write(ounit,*)' kalgainmat ',i1,kalgainmat(i1)
!!      enddo

!! 2. Calculation of the updated weights 
!!-----------------------------------------
!! update the weights: Eq. (5.35) p. 53 thesis SL
!! weights(n)=weights(n-1) + K(n)*[d(n)-h(i(n))]

      do i1=1,num_weights
!! MG: Fix2: Use scaled temporary error
        weights(i1)=weights(i1)+error_tmp*kalgainmat(i1)
      enddo
!!
!! 3. Calculation of the updated correlation matrix 
!!-------------------------------------------------
!! Calculation of the covariance matrix P(n):
!! P(n) = lambda^(-1)*P(n-1) -lambda^(-1)K(n)*J^T(n)*P(n-1)]
!!      = lambda^(-1)*[P(n-1)-K(n)*J^T(n)*P(n-1)]
!!
!! with (P*J)^T=J^T*P^T=J^T*P, P symmetric
!!
!! SUBROUTINE DSPR (UPLO, N, ALPHA, X, INCX, AP)
      alpha=-global_scaling
!!
!! Perform the symmetric rank 1 operation.
!!     A := alpha*x*x' + A,
!!  where alpha is a real scalar, x is an n element vector and A is an
!!  n by n symmetric matrix, supplied in packed form.
!! output is AP = corrmatrix
!!
!! my version

!!      call omp_set_num_threads(8)
!!      call mkl_set_num_threads(8)
      call dspr('l',num_weights,alpha,&
                  coh,1,corrmatrix)
       corrmatrix(:)=invlambda*corrmatrix(:)

!! debug
!!        do i1=1,corrdim
!!         write(ounit,*)'corrmat ',i1,corrmatrix(i1)
!!        enddo
!!        stop

!! 4. Calculation of the updated lambda 
!!-------------------------------------------------
!! update of the Kalman damping factor kalmanlambda
      kalmanlambda=kalmannue*kalmanlambda + 1.d0 - kalmannue
!!
      return
      end
