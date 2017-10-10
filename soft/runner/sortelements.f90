!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! determine elementindex and pairindex arrays

!! called by
!! - readinput.f90
!! - getdimensions.f90
!!
      subroutine sortelements()
!!
      use fileunits
      use globaloptions
      use nnflags
!!
      implicit none
!!
      integer counter                          ! internal
      integer ztemp                            ! internal
      integer i1,i2,i3                         ! internal
!!
      character*2 elementtemp                  ! internal
!!
!!
      elementindex(:)=0
!!
      if(nelem.gt.1)then
 10     continue
        do i1=1,nelem-1
          if(nucelem(i1).gt.nucelem(i1+1))then
            ztemp=nucelem(i1)
            elementtemp=element(i1)
            nucelem(i1)=nucelem(i1+1)
            element(i1)=element(i1+1)
            nucelem(i1+1)=ztemp
            element(i1+1)=elementtemp
            goto 10
          endif
        enddo
      endif
!!
      do i1=1,102
        do i2=1,nelem
          if(nucelem(i2).eq.i1)then
            elementindex(i1)=i2
          endif 
        enddo
      enddo
!!
      pairindex(:,:)=0
      counter       =0
      do i1=1,nelem
        do i2=1,nelem
         if(nucelem(i2).ge.nucelem(i1))then
            counter = counter +1
            pairindex(nucelem(i1),nucelem(i2)) = counter
            pairindex(nucelem(i2),nucelem(i1)) = counter
         endif
        enddo
      enddo
!!
      if(lnntb.and.nntb_flag(3))then      
        tripletindex(:,:,:)=0 ! is done after allocation already
        counter            =0
!! loop over all unique pairs
        do i1=1,nelem
          do i2=1,nelem
!! loop over all possible neighbor elements
            if(nucelem(i2).ge.nucelem(i1))then
              do i3=1,nelem
                counter = counter +1
                tripletindex(nucelem(i1),nucelem(i2),nucelem(i3)) = counter
                tripletindex(nucelem(i2),nucelem(i1),nucelem(i3)) = counter
              enddo
            endif
          enddo
        enddo
!! debugging
!!      do i1=1,3
!!        do i2=1,3
!!          do i3=1,3
!!            write(ounit,'(3i3,i4)')nucelem(i1),nucelem(i2),nucelem(i3),&
!!              tripletindex(nucelem(i1),nucelem(i2),nucelem(i3))
!!          enddo
!!        enddo
!!      enddo
!!      stop
      endif
!!
      return
      end
