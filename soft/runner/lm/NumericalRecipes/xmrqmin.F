      PROGRAM xmrqmin
C     driver for routine mrqmin
      INTEGER NPT,MA
      REAL SPREAD
      PARAMETER(NPT=100,MA=6,SPREAD=0.001)
      INTEGER i,ia(MA),idum,iter,itst,j,k,mfit
      REAL alamda,chisq,fgauss,gasdev,ochisq,x(NPT),y(NPT),sig(NPT),
     *     a(MA),covar(MA,MA),alpha(MA,MA),gues(MA)
      EXTERNAL fgauss
      DATA a/5.0,2.0,3.0,2.0,5.0,3.0/
      DATA gues/4.5,2.2,2.8,2.5,4.9,2.8/
      idum=-911
C     first try a sum of two Gaussians
      do 12 i=1,NPT
        x(i)=0.1*i
        y(i)=0.0
        do 11 j=1,MA,3
          y(i)=y(i)+a(j)*exp(-((x(i)-a(j+1))/a(j+2))**2)
11      continue
        y(i)=y(i)*(1.0+SPREAD*gasdev(idum))
        sig(i)=SPREAD*y(i)
12    continue
      mfit=MA
      do 13 i=1,mfit
        ia(i)=1
13    continue
      do 14 i=1,MA
        a(i)=gues(i)
14    continue
      do 16 iter=1,2
        alamda=-1
        call mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,
     *       MA,chisq,fgauss,alamda)
        k=1
        itst=0
1       write(*,'(/1x,a,i2,t18,a,f10.4,t43,a,e9.2)') 'Iteration #',k,
     *       'Chi-squared:',chisq,'ALAMDA:',alamda
        write(*,'(1x,t5,a,t13,a,t21,a,t29,a,t37,a,t45,a)') 'A(1)',
     *       'A(2)','A(3)','A(4)','A(5)','A(6)'
        write(*,'(1x,6f8.4)') (a(i),i=1,6)
        k=k+1
        ochisq=chisq
        call mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,
     *       MA,chisq,fgauss,alamda)
        if (chisq.gt.ochisq) then
          itst=0
        else if (abs(ochisq-chisq).lt.0.1) then
          itst=itst+1
        endif
        if (itst.lt.4) then
          goto 1
        endif
        alamda=0.0
        call mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,
     *       MA,chisq,fgauss,alamda)
        write(*,*) 'Uncertainties:'
        write(*,'(1x,6f8.4/)') (sqrt(covar(i,i)),i=1,6)
        write(*,'(1x,a)') 'Expected results:'
        write(*,'(1x,f7.2,5f8.2/)') 5.0,2.0,3.0,2.0,5.0,3.0
        if (iter.eq.1) then
          write(*,*) 'press return to continue with constraint'
          read(*,*)
          write(*,*) 'Holding a(2) and a(5) constant'
          do 15 j=1,MA
            a(j)=a(j)+.1
15        continue
          a(2)=2.0
          ia(2)=0
          a(5)=5.0
          ia(5)=0
        endif
16    continue
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *funcs,alamda)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,mfit
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        call covsrt(alpha,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        do 17 j=1,mfit
          do 16 k=1,mfit
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     *funcs)
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      EXTERNAL funcs
      PARAMETER (MMAX=20)
      INTEGER mfit,i,j,k,l,m
      REAL dy,sig2i,wt,ymod,dyda(MMAX)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      do 16 i=1,ndata
        call funcs(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL covar(npc,npc)
      INTEGER i,j,k
      REAL swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
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
                pause 'singular matrix in gaussj'
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
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
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

cccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE fgauss(x,a,y,dyda,na)
      INTEGER na
      REAL x,y,a(na),dyda(na)
      INTEGER i
      REAL arg,ex,fac
      y=0.
      do 11 i=1,na-1,3
        arg=(x-a(i+1))/a(i+2)
        ex=exp(-arg**2)
        fac=a(i)*ex*2.*arg
        y=y+a(i)*ex
        dyda(i)=ex
        dyda(i+1)=fac/a(i+2)
        dyda(i+2)=fac*arg/a(i+2)
11    continue
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

cccccccccccccccccccccccccccccccccccccc

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END



