!*> \brief \b DSPR
!*
!* =========== DOCUMENTATION ===========
!*
!* Online html documentation available at
!* http://www.netlib.org/lapack/explore-html/
!*
!* Definition:
!* ===========
!*
!* SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
!*
!* .. Scalar Arguments ..
!* DOUBLE PRECISION ALPHA
!* INTEGER INCX,N
!* CHARACTER UPLO
!* ..
!* .. Array Arguments ..
!* DOUBLE PRECISION AP(*),X(*)
!* ..
!*
!*
!*> \par Purpose:
!* =============
!*>
!*> \verbatim
!*>
!*> DSPR performs the symmetric rank 1 operation
!*>
!*> A := alpha*x*x**T + A,
!*>
!*> where alpha is a real scalar, x is an n element vector and A is an
!*> n by n symmetric matrix, supplied in packed form.
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
!*> \param[in,out] AP
!*> \verbatim
!*> AP is DOUBLE PRECISION array of DIMENSION at least
!*> ( ( n*( n + 1 ) )/2 ).
!*> Before entry with UPLO = 'U' or 'u', the array AP must
!*> contain the upper triangular part of the symmetric matrix
!*> packed sequentially, column by column, so that AP( 1 )
!*> contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!*> and a( 2, 2 ) respectively, and so on. On exit, the array
!*> AP is overwritten by the upper triangular part of the
!*> updated matrix.
!*> Before entry with UPLO = 'L' or 'l', the array AP must
!*> contain the lower triangular part of the symmetric matrix
!*> packed sequentially, column by column, so that AP( 1 )
!*> contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!*> and a( 3, 1 ) respectively, and so on. On exit, the array
!*> AP is overwritten by the lower triangular part of the
!*> updated matrix.
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
!*>
!*> -- Written on 22-October-1986.
!*> Jack Dongarra, Argonne National Lab.
!*> Jeremy Du Croz, Nag Central Office.
!*> Sven Hammarling, Nag Central Office.
!*> Richard Hanson, Sandia National Labs.
!*> \endverbatim
!*>
!* =====================================================================
      SUBROUTINE dspr(UPLO,N,ALPHA,X,INCX,AP)
!*
!* -- Reference BLAS level2 routine (version 3.4.0) --
!* -- Reference BLAS is a software package provided by Univ. of Tennessee, --
!* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!* November 2011
!*
!* .. Scalar Arguments ..
      DOUBLE PRECISION alpha
      INTEGER incx,n
      CHARACTER uplo
!* ..
!* .. Array Arguments ..
      DOUBLE PRECISION ap(*),x(*)
!* ..
!*
!* =====================================================================
!*
!* .. Parameters ..
      DOUBLE PRECISION zero
      parameter(zero=0.0d+0)
!* ..
!* .. Local Scalars ..
      DOUBLE PRECISION temp
      INTEGER i,info,ix,j,jx,k,kk,kx
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
      info = 5
      END IF
      IF (info.NE.0) THEN
      CALL xerbla('DSPR ',info)
      RETURN
      END IF
!*
!* Quick return if possible.
!*
      IF ((n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
!*
!* Set the start point in X if the increment is not unity.
!*
      IF (incx.LE.0) THEN
      kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
      kx = 1
      END IF
!*
!* Start the operations. In this version the elements of the array AP
!* are accessed sequentially with one pass through AP.
!*
      kk = 1
      IF (lsame(uplo,'U')) THEN
!*
!* Form A when upper triangle is stored in AP.
!*
      IF (incx.EQ.1) THEN
      DO 20 j = 1,n
      IF (x(j).NE.zero) THEN
      temp = alpha*x(j)
      k = kk
      DO 10 i = 1,j
      ap(k) = ap(k) + x(i)*temp
      k = k + 1
   10 CONTINUE
      END IF
      kk = kk + j
   20 CONTINUE
      ELSE
      jx = kx
      DO 40 j = 1,n
      IF (x(jx).NE.zero) THEN
      temp = alpha*x(jx)
      ix = kx
      DO 30 k = kk,kk + j - 1
      ap(k) = ap(k) + x(ix)*temp
      ix = ix + incx
   30 CONTINUE
      END IF
      jx = jx + incx
      kk = kk + j
   40 CONTINUE
      END IF
      ELSE
!*
!* Form A when lower triangle is stored in AP.
!*
      IF (incx.EQ.1) THEN
      DO 60 j = 1,n
      IF (x(j).NE.zero) THEN
      temp = alpha*x(j)
      k = kk
      DO 50 i = j,n
      ap(k) = ap(k) + x(i)*temp
      k = k + 1
   50 CONTINUE
      END IF
      kk = kk + n - j + 1
   60 CONTINUE
      ELSE
      jx = kx
      DO 80 j = 1,n
      IF (x(jx).NE.zero) THEN
      temp = alpha*x(jx)
      ix = jx
      DO 70 k = kk,kk + n - j
      ap(k) = ap(k) + x(ix)*temp
      ix = ix + incx
   70 CONTINUE
      END IF
      jx = jx + incx
      kk = kk + n - j + 1
   80 CONTINUE
      END IF
      END IF
!*
      RETURN
!*
!* End of DSPR .
!*
      END

