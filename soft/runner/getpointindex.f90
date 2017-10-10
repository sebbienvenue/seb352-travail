!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine getpointindex(ntrain,oseed,pointindex)
!!
      use fileunits
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer ntrain                ! in
      integer pointindex(ntrain)    ! out
      integer oseed                 ! in/out
      integer i1                    ! internal
      integer icount                ! internal
      integer ncount                ! internal
      integer npoints               ! internal
      integer ndone                 ! internal
      integer, dimension(:)  , allocatable :: idxtemp  ! internal
!!
      real*8 ran0                   ! internal
      real*8 z(ntrain)              ! internal
!!
!!
!! initialization
      do i1=1,ntrain
        pointindex(i1)=i1
      enddo
!!
!! get random order of training points
      if(lmixpoints)then
!! get ntrain random numbers
        do i1=1,ntrain
!!          write(ounit,'(a,i6,i14)')'JBgetpointindex ',i1,oseed
          z(i1)=ran0(oseed)
        enddo
!!
!! get new order of points
        call sort2realint(ntrain,z,pointindex)
!! debugging output
!!        write(ounit,*)'fully mixed points:'
!!        do i1=1,ntrain
!!          write(ounit,'(2i7)')i1,pointindex(i1)
!!        enddo
!!
!! OLD and inefficient:
!! get new order of points
!! 10     continue
!!        do i1=1,ntrain-1
!!          if(z(i1).gt.z(i1+1))then
!!            ztemp=z(i1)
!!            jtemp=pointindex(i1)
!!            z(i1)=z(i1+1)
!!            z(i1+1)=ztemp
!!            pointindex(i1)=pointindex(i1+1)
!!            pointindex(i1+1)=jtemp
!!            goto 10
!!          endif
!!        enddo
!! debugging output
!!        write(ounit,*)'fully mixed points:'
!!        do i1=1,ntrain
!!          write(ounit,'(2i7)')i1,pointindex(i1)
!!        enddo
!!
!!
!! now we need to sort within each block the indices for efficient file reading later
        ncount=ntrain
        icount=0
        ndone=0
 11     continue
        if(ncount.gt.nblock)then
          npoints=nblock
          ncount=ncount-nblock
        else
          npoints=ncount
          ncount=ncount-npoints
        endif
!! prepare temp array for subblock
        allocate(idxtemp(npoints))
        idxtemp(:)=0
!! copy subblock to temporary array
        do i1=1,npoints
          icount=icount+1
          idxtemp(i1)=pointindex(icount)
        enddo
!! sort points in subblock 
!!        call sortint(npoints,idxtemp) !! JB: if we use this line, different values of nblock yield different results
!! copy subblock back to full array
        do i1=1,npoints
          ndone=ndone+1
          pointindex(ndone)=idxtemp(i1)
!!          write(ounit,'(a,i4,i8)')'JBpointindex ',i1,pointindex(ndone)
        enddo
!!
        deallocate(idxtemp)
        if(ncount.gt.0) goto 11 ! do next block of points
      endif ! lmixpoints
!!
!! debugging output
!!      write(ounit,*)'sorted subblocks:'
!!      do i1=1,ntrain
!!        write(ounit,'(2i7)')i1,pointindex(i1)
!!      enddo
!!
!!
!! OLD INEFFICIENT SORTING:
!! now we need to sort within each block the indices for efficient file reading later
!!        icount=0
!!        jcount=0
!!        do i1=1,ntrain
!!          icount=icount+1
!! make a local copy of nblock indices
!!          itemp(icount)=pointindex(i1)
!!          if((icount.eq.nblock).or.(i1.eq.ntrain))then
!! sort this block of icount indices now
!! 30         continue
!!            do i2=1,icount-1
!!              if(itemp(i2).gt.itemp(i2+1))then
!!                jtemp           =itemp(i2)
!!                itemp(i2)       =itemp(i2+1)
!!                itemp(i2+1)     =jtemp
!!                goto 30
!!              endif
!!            enddo ! i2
!! copy sorted block back to pointindex array
!!            do i2=1,icount
!!              jcount=jcount+1
!!              pointindex(jcount)=itemp(i2)
!!            enddo
!!            icount=0
!!          endif
!!        enddo ! i1
!!      endif ! lmixpoints
!!
!! debugging output
!!      write(ounit,*)'sorted subblocks:'
!!      do i1=1,ntrain
!!        write(ounit,'(2i7)')i1,pointindex(i1)
!!      enddo
!!
      return
      end
