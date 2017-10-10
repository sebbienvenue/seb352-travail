!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - optimize_short_combined.f90
!! - optimize_short_combinedpair.f90
!!
      subroutine sorteshorterror(npoints,energyerror_list,lusee)
!!
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer npoints
      integer i1,i2
      integer istruct(nblock)                                                ! internal 
!!      integer select1
      real*8 select1
!!
      real*8 energyerror_list(nblock)                                        ! in
      real*8 energyerror_temp(npoints)                                       ! in
      real*8 energyerror_copy(npoints)                                       ! in
      real*8 temp                                                            ! internal
      real*8 ethres                                                          ! out
!!
      logical lusee(nblock)                                                  ! out 
!!
!! initializations
      lusee(:)=.false.
      do i1=1,npoints
        istruct(i1)=i1
        energyerror_temp(i1)=energyerror_list(i1)
      enddo
!!
!! determine the threshold for update
      temp=(1.d0-worste)*dble(npoints)
      i2=int(temp)
!!      write(ounit,*)'find point ',i2
      energyerror_copy(:)=energyerror_temp(:)
!! caution: select changes order of energyerror_copy
      ethres=select1(i2,npoints,energyerror_copy)
!!      write(ounit,*)'ethres ',i1,ethres
!! debug
!!      do i1=1,npoints
!!        if(energyerror_copy(i1).ge.ethres)then
!!          write(ounit,*)i1,energyerror_copy(i1)
!!        endif
!!      enddo
!!
!!
!! set lusee array
      do i1=1,npoints
        if(energyerror_temp(i1).gt.ethres)then
          lusee(istruct(i1))=.true.
        endif
      enddo
!!
!! debug
!!      write(ounit,*)'lusee array'
!!      do i1=1,npoints
!!        write(ounit,'(i6,f14.8,l)')i1,&
!!          energyerror_list(i1),lusee(i1)
!!      enddo
!!
!!
!!      stop
!!
      return
      end
