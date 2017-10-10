      DOUBLE PRECISION FUNCTION PDDOT( N, X, Y )
!
!     -- BLACS example code --
!     Written by Clint Whaley 7/26/94.
!     ..
!     .. Scalar Arguments ..
      INTEGER N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION X(*), Y(*)
!     ..
!
!  Purpose
!  =======
!  PDDOT is a restricted parallel version of the BLAS routine
!  DDOT.  It assumes that the increment on both vectors is one,
!  and that process {0,0} starts out owning the vectors and
!  has N.  It returns the dot product of the two N-length vectors
!  X and Y, i.e. PDDOT = X' Y.
!
!  Arguments
!  =========
!  N            (input/ouput) INTEGER
!               The length of the vectors X and Y.  Input
!               for {0,0}, output for everyone else.
!
!  X            (input/output) DOUBLE PRECISION array of dimension (N)
!               The vector X of PDDOT = X' Y.  Input for {0,0},
!               output for everyone else.
!
!  Y            (input/output) DOUBLE PRECISION array of dimension (N)
!               The vector Y of PDDOT = X' Y.  Input for {0,0},
!               output for everyone else.
!
!  =====================================================================
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL GRIDINFO, DGEBS2D, DGEBR2D, DGSUM2D
!     ..
!     .. Local Scalars ..
      INTEGER IAM, NPROCS, NPROW, NPCOL, MYPROW, MYPCOL, I, LN
      DOUBLE PRECISION LDDOT
!     ..
!     .. Executable Statements ..
!
!     Find out what grid has been set up, and pretend it's 1-D
!
      CALL GRIDINFO( NPROW, NPCOL, MYPROW, MYPCOL )
      IAM = MYPROW*NPCOL + MYPCOL
      NPROCS = NPROW * NPCOL
!
!     Do bone-headed thing, and just send entire X and Y to
!     everyone
!
      IF ( (MYPROW.EQ.0) .AND. (MYPCOL.EQ.0) ) THEN
         CALL IGEBS2D('All', 'i-ring', 1, 1, N, 1 )
         CALL DGEBS2D('All', 'i-ring', N, 1, X, N )
         CALL DGEBS2D('All', 'i-ring', N, 1, Y, N )
      ELSE
         CALL IGEBR2D('All', 'i-ring', 1, 1, N, 1, 0, 0 )
         CALL DGEBR2D('All', 'i-ring', N, 1, X, N, 0, 0 )
         CALL DGEBR2D('All', 'i-ring', N, 1, Y, N, 0, 0 )
      ENDIF
!
!     Find out the number of local rows to multiply (LN), and
!     where in vectors to start (I)
!
      LN = N / NPROCS
      I = 1 + IAM * LN
!
!     Last process does any extra rows
!
      IF (IAM .EQ. NPROCS-1) LN = LN + MOD(N, NPROCS)
!
!     Figure dot product of my piece of X and Y
!
      LDDOT = DDOT( LN, X(I), 1, Y(I), 1 )
!
!     Add local dot products to get global dot product;
!     give all procs the answer
!
      CALL DGSUM2D( 'All', '1-tree', 1, 1, LDDOT, 1, -1, 0 )

      PDDOT = LDDOT

      RETURN
      END
