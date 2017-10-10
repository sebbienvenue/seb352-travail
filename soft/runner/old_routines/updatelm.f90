!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine updatelm(&
        dimlm,nelem,maxnum_weightsewald,&
        alphalm,covarlm,betalm,dalm,chisq,alambdalm,&
        weights_ewald,&
        ldebug)
!!
      use fileunits
!!
      implicit none
!!
      integer dimlm
      integer i1,i2,i3,i4  
      integer maxnum_weightsewald
      integer nelem

      real*8 alphalm(dimlm,dimlm,nelem)
      real*8 covarlm(dimlm,dimlm,nelem)
      real*8 betalm(dimlm,nelem)
      real*8 chisq(nelem)
      real*8 weights_ewald(maxnum_weightsewald,nelem)
      real*8 alambdalm(nelem)
      real*8 dalm(dimlm,nelem)

      logical ldebug
!!
!! debug
!!      alambdalm(:)=0.0d0
!! debug
!!      do i3=1,nelem
!!        do i1=1,dimlm
!!          do i2=1,dimlm
!!            if(alphalm(i1,i2,i3).ne.alphalm(i2,i1,i3))then
!!              write(*,*)'Error '
!!              stop
!!            endif
!!          enddo ! i2
!!        enddo ! i1
!!      enddo ! i3
!!
!!
      do i3=1,nelem 
        write(ounit,*)'LM alambdalm ',alambdalm(i3)
        do i1=1,dimlm
          do i2=1,dimlm
            covarlm(i1,i2,i3)=alphalm(i1,i2,i3)
          enddo ! i2
          covarlm(i1,i1,i3)=alphalm(i1,i1,i3)*(1.d0+alambdalm(i3))
          dalm(i1,i3)=betalm(i1,i3)
        enddo ! i1
      enddo ! i3
!!
!!      write(ounit,*)dalm
!!      write(ounit,*)covarlm
!!
!! Gaussian Jordan matrix inversion and linear equation solution
      do i1=1,nelem
        write(ounit,*)'LM calling gaussj ',i1
        call gaussj(covarlm(1,1,i1),dimlm,dimlm,dalm(1,i1),1,1)
      enddo
!!
      do i1=1,dimlm
        do i2=1,nelem
          weights_ewald(i1,i2)=weights_ewald(i1,i2)+dalm(i1,i2)
        enddo ! i2
      enddo ! i1
!!
      return
      end

!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE gaussj(a,n,np,b,m,mp)
!!
      use fileunits
!!
      implicit none
!!
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=5000)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
! initialize ipiv
      do 11 j=1,n
        ipiv(j)=0
11    continue
! loop over all weight parameters
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(ounit,*) 'Error: singular matrix in gaussj ipiv ',i,j,k,ipiv(k)
                stop
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) then
          write(ounit,*)'Error: singular matrix in gaussj a'
          stop
        endif
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

