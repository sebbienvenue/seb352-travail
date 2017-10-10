      program dgemmtest

! compile on breda with
!ifort -free dgemmtest.f90  -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t -lpthread
! fort -free -check all -traceback dgemmtest.f90  -L/opt/intel/mkl721/lib/em64t -lguide -lmkl_em64t -lpthread
      implicit none

      integer m,n,k
      integer lda,ldb,ldc
      
      real*8 one
!      real*8 a(1,3)
      real*8 a(2,3)
!      real*8 b(3,1)
      real*8 b(3,2)
!      real*8 c(1)
      real*8 c(2,2)
      real*8 alpha,beta

      one=1.0d0
      alpha=one
      beta=one

! matrix A
! ( 1 2 3
!   4 5 6 )

      a(1,1)=1.0d0
      a(1,2)=2.0d0
      a(1,3)=3.0d0
      a(2,1)=4.0d0
      a(2,2)=5.0d0
      a(2,3)=6.0d0

! matrix B
! ( 7 8 
!   9 10
!  11 12 )

      b(1,1)=7.0d0
      b(1,2)=8.0d0
      b(2,1)=9.0d0
      b(2,2)=10.0d0
      b(3,1)=11.0d0
      b(3,2)=12.0d0

! matrix C to be added to A*B

      c(1,1)=0.0d0
      c(2,1)=0.0d0
      c(1,2)=0.0d0
      c(2,2)=0.0d0

! result
! c_11 = 1*7 + 2*9 + 3* 11 = 58
! c_12 = 1*8 + 2*10 + 3*12 = 64
! c_21 = 4*7 + 5*9 + 6*11 = 139
! c_22 = 4*8 + 5*10 + 6*12 = 154

! matrix C result
! ( 58 64
!   139 154 )

! note: with write(*,*) a the elements are written column-wise


!       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!      C := alpha*op( A )*op( B ) + beta*C,
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.

      m=2
      n=2
      k=3
! the leading dimension is the number of numbers per column
      lda=2
      ldb=3
      ldc=2

!      m=1
!      n=1
!      k=3
!      lda=3
!      ldb=3
!      ldc=3
      call dgemm('n','n',m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

      write(*,*)'a'
      write(*,*)a
      write(*,*)'b'
      write(*,*)b
      write(*,*)'c'
      write(*,*)c

      end
