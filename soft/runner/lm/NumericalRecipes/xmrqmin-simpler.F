      PROGRAM xmrqmin
C     driver for routine mrqmin
      INTEGER ndata,ma
      REAL SPREAD
      PARAMETER(ndata=100,ma=6,SPREAD=0.001)
      INTEGER i,idum,iter,itst,j,k
      REAL alamda,chisq,fgauss,gasdev,ochisq,x(ndata),y(ndata),
     *     sig(ndata),
     *     a(ma),covar(ma,ma),alpha(ma,ma),gues(ma)
      EXTERNAL fgauss
      DATA a/5.0,2.0,3.0,2.0,5.0,3.0/
C initial guess
      DATA gues/4.5,2.2,2.8,2.5,4.9,2.8/

C set seed for random noise
      idum=-911
C calculate sum of two Gaussians as target function
      do 12 i=1,ndata
        x(i)=0.1*i
        y(i)=0.0
        do 11 j=1,ma,3
          y(i)=y(i)+a(j)*exp(-((x(i)-a(j+1))/a(j+2))**2)
11      continue
C add noise to y(i) to mimic real noisy data:
        y(i)=y(i)*(1.0+SPREAD*gasdev(idum))
        sig(i)=SPREAD*y(i)
12    continue

c prepare initial guess
      do 14 i=1,ma
        a(i)=gues(i)
14    continue

c initializations
        alamda=-1
        call mrqmin(x,y,sig,ndata,a,ma,covar,alpha,
     *       chisq,fgauss,alamda)

c initialize number of iterations
        k=1
c initialize number of good steps
        itst=0
1       write(*,'(/1x,a,i2,t18,a,f10.4,t43,a,e9.2)') 'Iteration #',k,
     *       'Chi-squared:',chisq,'ALAMDA:',alamda
        write(*,'(1x,t5,a,t13,a,t21,a,t29,a,t37,a,t45,a)') 'A(1)',
     *       'A(2)','A(3)','A(4)','A(5)','A(6)'
        write(*,'(1x,6f8.4)') (a(i),i=1,6)
c start new iteration
        k=k+1
        ochisq=chisq
        call mrqmin(x,y,sig,ndata,a,ma,covar,alpha,
     *       chisq,fgauss,alamda)
c check if step was in right direction
        if (chisq.gt.ochisq) then
          itst=0
        else if (abs(ochisq-chisq).lt.0.1) then
          itst=itst+1
        endif
c check convergence
        if (itst.lt.4) then
          goto 1
        endif

c get final results if converged
        alamda=0.0
        call mrqmin(x,y,sig,ndata,a,ma,covar,alpha,
     *       chisq,fgauss,alamda)
        write(*,*) 'Uncertainties:'
        write(*,'(1x,6f8.4/)') (sqrt(covar(i,i)),i=1,6)
        write(*,'(1x,a)') 'Expected results:'
        write(*,'(1x,f7.2,5f8.2/)') 5.0,2.0,3.0,2.0,5.0,3.0
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE mrqmin(x,y,sig,ndata,a,ma,covar,alpha,chisq,
     *funcs,alamda)

      INTEGER ma,ndata,MMAX
      REAL alamda,chisq,funcs,a(ma),alpha(ma,ma),covar(ma,ma),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20)
      INTEGER j,k,l
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da

c initialization only in first iteration
      if(alamda.lt.0.)then
        alamda=0.001
c get alpha, beta and chisq:
        call mrqcof(x,y,sig,ndata,a,ma,alpha,beta,ma,chisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif

c get da and covar from alpha and beta
      do 14 j=1,ma
        do 13 k=1,ma
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
14    continue

c Gauss-Jordan matrix inversion and linear equation solution
c get covar and da
      call gaussj(covar,ma,da)

      if(alamda.eq.0.)then
        return
      endif

c update atry
      do 15 l=1,ma
        atry(l)=a(l)+da(l)
15    continue

c get covar, da and chisq for new atry:
      call mrqcof(x,y,sig,ndata,atry,ma,covar,da,ma,chisq,funcs)

c good step: decrease alambda, keep covar, da and a
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        do 17 j=1,ma
          do 16 k=1,ma
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
c bad step: increase alambda, discard covar, da and a
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c get covar, da and chisq for new atry:
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ma,alpha,beta,nalp,chisq,
     *funcs)

      INTEGER ma,nalp,ndata
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      EXTERNAL funcs
      INTEGER i,j,k,l,m
      REAL dy,sig2i,wt,ymod,dyda(ma)

c initialize beta and triangle of alpha 
      do 13 j=1,ma
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue

      chisq=0.

c loop over all training points
      do 16 i=1,ndata
c get function value for current set of a(i)
        call funcs(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
c get error of this point
        dy=y(i)-ymod

c set up beta and triangle of alpha
        do 15 l=1,ma
          wt=dyda(l)*sig2i
          do 14 m=1,l
            alpha(l,m)=alpha(l,m)+wt*dyda(m)
14        continue
          beta(l)=beta(l)+dy*wt
15      continue
c get squared error
        chisq=chisq+dy*dy*sig2i
16    continue

c symmetrize alpha
      do 18 j=2,ma
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE gaussj(covar,ma,b)
      INTEGER ma
      REAL covar(ma,ma),b(ma)
      INTEGER i,icol,irow,j,k,l,ll,indxc(ma),indxr(ma),ipiv(ma)
      REAL big,dum,pivinv

      do 11 j=1,ma
        ipiv(j)=0
11    continue

      do 22 i=1,ma
        big=0.
        do 13 j=1,ma
          if(ipiv(j).ne.1)then
            do 12 k=1,ma
              if (ipiv(k).eq.0) then
                if (abs(covar(j,k)).ge.big)then
                  big=abs(covar(j,k))
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
          do 14 l=1,ma
            dum=covar(irow,l)
            covar(irow,l)=covar(icol,l)
            covar(icol,l)=dum
14        continue

          dum=b(irow)
          b(irow)=b(icol)
          b(icol)=dum
        endif

        indxr(i)=irow
        indxc(i)=icol

        if (covar(icol,icol).eq.0.) pause 'singular matrix in gaussj'

        pivinv=1./covar(icol,icol)
        covar(icol,icol)=1.

        do 16 l=1,ma
          covar(icol,l)=covar(icol,l)*pivinv
16      continue
        b(icol)=b(icol)*pivinv
        do 21 ll=1,ma
          if(ll.ne.icol)then
            dum=covar(ll,icol)
            covar(ll,icol)=0.
            do 18 l=1,ma
              covar(ll,l)=covar(ll,l)-covar(icol,l)*dum
18          continue
            b(ll)=b(ll)-b(icol)*dum
          endif
21      continue
22    continue

      do 24 l=ma,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,ma
            dum=covar(k,indxr(l))
            covar(k,indxr(l))=covar(k,indxc(l))
            covar(k,indxc(l))=dum
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



