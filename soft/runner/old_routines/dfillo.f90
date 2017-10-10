	subroutine dfillo ( fillout, uplo, n, nsq, ap )
****************************************************************************
*                                                                          *
*   DATA PARALLEL BLAS based on MPL                                        *
*                                                                          *
*   Version 1.0   1/9-92 ,                                                 *
*   For MasPar MP-1 computers                                              *
*                                                                          *
*   para//ab, University of Bergen, NORWAY                                 *
*                                                                          *
*   These programs must be called using F90 style array syntax.            *
*   Note that the F77 style calling sequence has been retained             *
*   in this version for compatibility reasons, be aware that               *
*   parameters related to the array dimensions and shape therefore may     *
*   be redundant and without any influence.                                *
*   The calling sequence may be changed in a future version.               *
*   Please report any BUGs, ideas for improvement or other                 *
*   comments to                                                            *
*                    adm@parallab.uib.no                                   *
*                                                                          *
*   Future versions may then reflect your suggestions.                     *
*   The most current version of this software is available                 *
*   from netlib@nac.no , send the message `send index from maspar'         *
*                                                                          *
*   REVISIONS:                                                             *
*                                                                          *
****************************************************************************
c
	integer n, nsq
	double precision, array(:) :: ap
	double precision, array(n,n) :: fillout
	character*1 uplo
	logical  lsame
	external lsame
c
c   local variables
c
	integer mold(2), k
	integer, array(nsq) :: nrow, ncol, whereis, iseqn
	double precision, array(nsq) :: fill1d
*
*     Set up the where array
*
cts	do i = 1 , nsq
cts		iseqn(i) = i
cts	enddo
	forall (i=1:nsq) iseqn(i) = i
	ncol = iseqn - 1
	ncol = ncol / n
 	nrow = iseqn - ncol * n
	ncol = ncol + 1

	if ( lsame (uplo, 'U') ) then
	   where (nrow .lt. ncol)
	      whereis = ncol
	      whereis = (whereis * (whereis-1)) / 2
	      whereis = whereis + nrow
	   elsewhere
	      whereis = nrow
	      whereis = (whereis * (whereis-1)) / 2
	      whereis = whereis + ncol
	   endwhere
c
c   so whereis, in case, n = 4,  now stores
c
c	1	2	4	7
c	2	3	5	8
c	4	5	6	9
c	7	8	9      10
c
	else
c
c  uplo = 'l'
c
	   where (nrow .ge. ncol)
	      whereis = n + 1 - ncol
	      whereis = (n*(n-1)) / 2 - ((whereis*(whereis-1)) / 2)
	      whereis = whereis + nrow
	   elsewhere
	      whereis = n + 1 - nrow
	      whereis = (n*(n-1)) / 2 - ((whereis*(whereis-1)) / 2)
	      whereis = whereis + ncol
	   endwhere
c
c   so if whereis is 4 by 4 then it now stores
c
c	1	2	3	4
c	2	5	6	7
c	3	6	8	9
c	4	7	9      10
c
	endif
c
c   set up fulla by getting data from ap(whereis)
c
	fill1d = ap(whereis)
c
	mold(:) = n
clater 	fillout = reshape(mold,fill1d)
        k = 1
        do j = 1, n
            fillout(1:n,j) = fill1d(k:k+n-1)
            k = k + n
        enddo
c
	return
	end
