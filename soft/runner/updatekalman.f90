!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - optimize_short.f90
!! - optimize_short_combinedpair.f90
!! - optimize_ewald.f90
!! - fitforcesshort.f90
!!
      subroutine updatekalman(kaldim,corrdim,&
                num_weights,kalmanlambda,kalmannue,&
                weights,dedw,corrmatrix,error)
!!
      use fileunits
!! Don't use fittingoptions, because updatekalman is used for different quantities
!!
      implicit none
!!
      integer num_weights                                         ! in
      integer i1                                                  ! internal
      integer kaldim                                              ! in
      integer corrdim                                             ! in
!! 
      real*8 kalmanlambda                                         ! in/out
      real*8 kalmannue                                            ! in
      real*8 weights(num_weights)                                 ! in/out
      real*8 dedw(num_weights)                                    ! in
      real*8 kalgainmat(kaldim)                                   ! internal
      real*8 corrmatrix(corrdim)                                  ! in/out
      real*8 error                                                ! in
      real*8 alpha                                                ! internal
      real*8 inverse                                              ! internal
      real*8 invlambda                                            ! internal
      real*8 coh(num_weights)                                     ! internal
      real*8 ddot                                                 ! internal
!!
!! debugging:
!!      write(*,*)'kalmanlambda ',kalmanlambda
!!      write(*,*)'kalmannue ',kalmannue
!!
!! initialization
      coh(:)=0.0d0
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
      call dspmv('l',num_weights,1.d0,corrmatrix,dedw,1,0.d0,coh,1)

!! result: coh = P(n-1)*J(n)
!!
!! Calculation of K(n)=lambda^(-1)*P(n-1)*J(n)*[I+lambda^(-1)*J^T(n)*P(n-1)*J(n)]^(-1) ! thesis SL page 53
!!                    =lambda^(-1)*coh        *[I+lambda^(-1)*J^T(n)*P(n-1)*J(n)]^(-1) ! thesis SL page 53
!!                    = coh *[lambda + J^T(n)*P(n-1)*J(n)]^(-1)
!! 1. Matrix Inversion including weight matrix
!! ddot=Calculation of J^T(n)*P(n-1)*J(n)
!!       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!! forms the dot product of two vectors: dedw and coh
      inverse = 1.d0/(kalmanlambda+ddot(num_weights,dedw,1,coh,1))

!!
!! calculation of K(n) in thesis SL page 53
!! multiplication by lambda^-1

      kalgainmat(:)=inverse*coh(:) 

!! debug
!!      do i1=1,num_weights
!!        write(ounit,*)' kalgainmat ',i1,kalgainmat(i1)
!!      enddo


!! result: kalgainmat is now complete
!!
!!
!! 2. Calculation of the updated weights 
!!-----------------------------------------
!! update the weights: Eq. (5.35) p. 53 thesis SL
!! weights(n)=weights(n-1) + K(n)*[d(n)-h(i(n))]


      do i1=1,num_weights
        weights(i1)=weights(i1)+error*kalgainmat(i1)
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
      alpha=-inverse
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
