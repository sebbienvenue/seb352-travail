!*> \brief \b DSPMV
!*
!* =========== DOCUMENTATION ===========
!*
!* Online html documentation available at
!* http://www.netlib.org/lapack/explore-html/
!*
!* Definition:
!* ===========
!*
!* SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
!*
!* .. Scalar Arguments ..
!* DOUBLE PRECISION ALPHA,BETA
!* INTEGER INCX,INCY,N
!* CHARACTER UPLO
!* ..
!* .. Array Arguments ..
!* DOUBLE PRECISION AP(*),X(*),Y(*)
!* ..
!*
!*
!*> \par Purpose:
!* =============
!*>
!*> \verbatim
!*>
!*> DSPMV performs the matrix-vector operation
!*>
!*> y := alpha*A*x + beta*y,
!*>
!*> where alpha and beta are scalars, x and y are n element vectors and
!*> A is an n by n symmetric matrix, supplied in packed form.
!*> \endverbatim
!*
!* Arguments:
!* ==========
!*
!*> \param[in] UPLO
!*> \verbatim
!*> UPLO is CHARACTER*1
!*> On entry, UPLO specifies whether the upper or lower
!*> triangular part of the matrix A is supplied in the packed
!*> array AP as follows:
!*>
!*> UPLO = 'U' or 'u' The upper triangular part of A is
!*> supplied in AP.
!*>
!*> UPLO = 'L' or 'l' The lower triangular part of A is
!*> supplied in AP.
!*> \endverbatim
!*>
!*> \param[in] N
!*> \verbatim
!*> N is INTEGER
!*> On entry, N specifies the order of the matrix A.
!*> N must be at least zero.
!*> \endverbatim
!*>
!*> \param[in] ALPHA
!*> \verbatim
!*> ALPHA is DOUBLE PRECISION.
!*> On entry, ALPHA specifies the scalar alpha.
!*> \endverbatim
!*>
!*> \param[in] AP
!*> \verbatim
!*> AP is DOUBLE PRECISION array of DIMENSION at least
!*> ( ( n*( n + 1 ) )/2 ).
!*> Before entry with UPLO = 'U' or 'u', the array AP must
!*> contain the upper triangular part of the symmetric matrix
!*> packed sequentially, column by column, so that AP( 1 )
!*> contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!*> and a( 2, 2 ) respectively, and so on.
!*> Before entry with UPLO = 'L' or 'l', the array AP must
!*> contain the lower triangular part of the symmetric matrix
!*> packed sequentially, column by column, so that AP( 1 )
!*> contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!*> and a( 3, 1 ) respectively, and so on.
!*> \endverbatim
!*>
!*> \param[in] X
!*> \verbatim
!*> X is DOUBLE PRECISION array of dimension at least
!*> ( 1 + ( n - 1 )*abs( INCX ) ).
!*> Before entry, the incremented array X must contain the n
!*> element vector x.
!*> \endverbatim
!*>
!*> \param[in] INCX
!*> \verbatim
!*> INCX is INTEGER
!*> On entry, INCX specifies the increment for the elements of
!*> X. INCX must not be zero.
!*> \endverbatim
!*>
!*> \param[in] BETA
!*> \verbatim
!*> BETA is DOUBLE PRECISION.
!*> On entry, BETA specifies the scalar beta. When BETA is
!*> supplied as zero then Y need not be set on input.
!*> \endverbatim
!*>
!*> \param[in,out] Y
!*> \verbatim
!*> Y is DOUBLE PRECISION array of dimension at least
!*> ( 1 + ( n - 1 )*abs( INCY ) ).
!*> Before entry, the incremented array Y must contain the n
!*> element vector y. On exit, Y is overwritten by the updated
!*> vector y.
!*> \endverbatim
!*>
!*> \param[in] INCY
!*> \verbatim
!*> INCY is INTEGER
!*> On entry, INCY specifies the increment for the elements of
!*> Y. INCY must not be zero.
!*> \endverbatim
!*
!* Authors:
!* ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2011
!*
!*> \ingroup double_blas_level2
!*
!*> \par Further Details:
!* =====================
!*>
!*> \verbatim
!*>
!*> Level 2 Blas routine.
!*> The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*> -- Written on 22-October-1986.
!*> Jack Dongarra, Argonne National Lab.
!*> Jeremy Du Croz, Nag Central Office.
!*> Sven Hammarling, Nag Central Office.
!*> Richard Hanson, Sandia National Labs.
!*> \endverbatim
!*>
!* =====================================================================
      SUBROUTINE dspmv(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
!*
!* -- Reference BLAS level2 routine (version 3.4.0) --
!* -- Reference BLAS is a software package provided by Univ. of Tennessee, --
!* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!* November 2011
!*
!* .. Scalar Arguments ..
      DOUBLE PRECISION alpha,beta
      INTEGER incx,incy,n
      CHARACTER uplo
!* ..
!* .. Array Arguments ..
      DOUBLE PRECISION ap(*),x(*),y(*)
!* ..
!*
!* =====================================================================
!*
!* .. Parameters ..
      DOUBLE PRECISION one,zero
      parameter(one=1.0d+0,zero=0.0d+0)
!* ..
!* .. Local Scalars ..
      DOUBLE PRECISION temp1,temp2
      INTEGER i,info,ix,iy,j,jx,jy,k,kk,kx,ky
!* ..
!* .. External Functions ..
      LOGICAL lsame
      EXTERNAL lsame
!* ..
!* .. External Subroutines ..
     EXTERNAL xerbla
!* ..
!*
!* Test the input parameters.
!*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
      info = 1
      ELSE IF (n.LT.0) THEN
      info = 2
      ELSE IF (incx.EQ.0) THEN
      info = 6
      ELSE IF (incy.EQ.0) THEN
      info = 9
      END IF
      IF (info.NE.0) THEN
      CALL xerbla('DSPMV ',info)
      RETURN
      END IF
!*
!* Quick return if possible.
!*
      IF ((n.EQ.0) .OR. ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
!*
!* Set up the start points in X and Y.
!*
      IF (incx.GT.0) THEN
      kx = 1
      ELSE
      kx = 1 - (n-1)*incx
      END IF
      IF (incy.GT.0) THEN
      ky = 1
      ELSE
      ky = 1 - (n-1)*incy
      END IF
!*
!* Start the operations. In this version the elements of the array AP
!* are accessed sequentially with one pass through AP.
!*
!* First form y := beta*y.
!*
      IF (beta.NE.one) THEN
      IF (incy.EQ.1) THEN
      IF (beta.EQ.zero) THEN
      DO 10 i = 1,n
      y(i) = zero
   10 CONTINUE
      ELSE
      DO 20 i = 1,n
      y(i) = beta*y(i)
   20 CONTINUE
      END IF
      ELSE
      iy = ky
      IF (beta.EQ.zero) THEN
      DO 30 i = 1,n
      y(iy) = zero
      iy = iy + incy
   30 CONTINUE
      ELSE
      DO 40 i = 1,n
      y(iy) = beta*y(iy)
      iy = iy + incy
   40 CONTINUE
      END IF
      END IF
      END IF
      IF (alpha.EQ.zero) RETURN
      kk = 1
      IF (lsame(uplo,'U')) THEN
!*
!* Form y when AP contains the upper triangle.
!*
      IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
      DO 60 j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      k = kk
      DO 50 i = 1,j - 1
      y(i) = y(i) + temp1*ap(k)
      temp2 = temp2 + ap(k)*x(i)
      k = k + 1
   50 CONTINUE
      y(j) = y(j) + temp1*ap(kk+j-1) + alpha*temp2
      kk = kk + j
   60 CONTINUE
      ELSE
      jx = kx
      jy = ky
      DO 80 j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      ix = kx
      iy = ky
      DO 70 k = kk,kk + j - 2
      y(iy) = y(iy) + temp1*ap(k)
      temp2 = temp2 + ap(k)*x(ix)
      ix = ix + incx
      iy = iy + incy
   70 CONTINUE
      y(jy) = y(jy) + temp1*ap(kk+j-1) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
      kk = kk + j
   80 CONTINUE
      END IF
      ELSE
!*
!* Form y when AP contains the lower triangle.
!*
      IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
      DO 100 j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      y(j) = y(j) + temp1*ap(kk)
      k = kk + 1
      DO 90 i = j + 1,n
      y(i) = y(i) + temp1*ap(k)
      temp2 = temp2 + ap(k)*x(i)
      k = k + 1
   90 CONTINUE
      y(j) = y(j) + alpha*temp2
      kk = kk + (n-j+1)
  100 CONTINUE
      ELSE
      jx = kx
      jy = ky
      DO 120 j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      y(jy) = y(jy) + temp1*ap(kk)
      ix = jx
      iy = jy
      DO 110 k = kk + 1,kk + n - j
      ix = ix + incx
      iy = iy + incy
      y(iy) = y(iy) + temp1*ap(k)
      temp2 = temp2 + ap(k)*x(ix)
  110 CONTINUE
      y(jy) = y(jy) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
      kk = kk + (n-j+1)
  120 CONTINUE
      END IF
      END IF
!*
      RETURN
!*
!* End of DSPMV .
!*
      END
