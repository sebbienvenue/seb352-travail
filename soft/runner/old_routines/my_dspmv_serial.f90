      subroutine my_dspmv_serial(uplo,num_weightsshort,corrmatrix,dedw,coh)
!!
      implicit none
!!
      integer num_weightsshort 
      integer i1,i2,i3
      integer nrow
      integer ncol
      integer iseqn(num_weightsshort*num_weightsshort)
      integer whereis(num_weightsshort*num_weightsshort)
!!
      real*8 corrmatrix(num_weightsshort*(num_weightsshort-1)/2)
      real*8 dedw(num_weightsshort)
      real*8 coh(num_weightsshort)
      real*8 fullcorrmatrix(num_weightsshort,num_weightsshort)
      real*8 fill1d(num_weightsshort*num_weightsshort)
!!
      character*1 uplo   
!!
!!      call dfillo(fullcorrmatrix,uplo,num_weightsshort,num_weightsshort*num_weightsshort,corrmatrix) 
!!
!! my version of dfillo:
!!      forall(i1=1:num_weightsshort**2) iseqn(i1)=i1
!!        ncol=iseqn - 1
!!        ncol=ncol/num_weightsshort
!!        nrow=iseqn - ncol*num_weightsshort
!!        ncol=ncol+1
!!
!!        where(nrow.ge.ncol)
!!          whereis = num_weightsshort + 1 - ncol
!!          whereis = (num_weightsshort*(num_weightsshort-1))/2 - ((whereis*(whereis-1))/2)
!!          whereis = whereis + nrow
!!        elsewhere
!!          whereis = num_weightsshort +1 - nrow
!!          whereis = (num_weightsshort*(num_weightsshort-1))/2 - ((whereis*(whereis-1))/2)
!!          whereis = whereis + ncol
!!        endwhere
!!
!!        fill1d = corrmatrix(whereis)
!!
!!      end forall
!!
!!
      do i1=1,num_weightsshort**2
        ncol = i1 - 1
        ncol = ncol/num_weightsshort
        nrow = i1 - ncol*num_weightsshort
        ncol = ncol + 1
        if(nrow.ge.ncol)then
          whereis(i1) = num_weightsshort + 1 - ncol
          whereis(i1) = (num_weightsshort*(num_weightsshort-1))/2 - ((whereis(i1)*(whereis(i1)-1))/2)
          whereis(i1) = whereis(i1) + nrow
        else
          whereis(i1) = num_weightsshort +1 - nrow
          whereis(i1) = (num_weightsshort*(num_weightsshort-1))/2 - ((whereis(i1)*(whereis(i1)-1))/2)
          whereis(i1) = whereis(i1) + ncol
        endif
        fill1d(i1)=corrmatrix(whereis(i1))
      enddo ! i1
!!
      i2 = 1
      do i3=1,num_weightsshort
        fullcorrmatrix(1:num_weightsshort,i2)=fill1d(i2:i2+num_weightsshort-1)
        i2=i2+num_weightsshort
      enddo ! i3
!!
      call dgemv('N',num_weightsshort,num_weightsshort,1.0d0,fullcorrmatrix,&
        num_weightsshort,dedw,1,0.0d0,coh,1)
!!
      return 
      end 
