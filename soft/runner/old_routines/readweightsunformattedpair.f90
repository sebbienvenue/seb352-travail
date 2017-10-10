!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - fitting.f90 (2x)
!!            - predictionpair.f90 (2x)
!!
      subroutine readweightsunformattedpair(nelem,npairs,nucelem,&
                 maxnum_weightsp,num_weightsp,weightsp,&
                 lewald,lshort)
!!
      use fileunits
!!
      implicit none
!!
      integer nelem                         ! in
      integer npairs                        ! in
      integer maxnum_weightsp               ! in
      integer num_weightsp(npairs)          ! in
      integer i1,i2,i3                         ! internal
      integer nucelem(nelem)                ! in
      integer paircount                     ! internal
!! 
      real*8 weightsp(maxnum_weightsp,npairs)  ! out
!!
      character*20 filename
!!
      logical lewald
      logical lshort
      logical lexist
!!
      paircount = 0
      if(lshort)then
        do i1=1,nelem
          do i2=i1,nelem
             paircount = paircount + 1
             filename='weights.000.000.data'
             if(nucelem(i1).gt.99)then
               write(filename(9:11),'(i3)')nucelem(i1)
             elseif(nucelem(i1).gt.9)then
               write(filename(10:11),'(i2)')nucelem(i1)
             else
               write(filename(11:11),'(i1)')nucelem(i1)
             endif

             if(nucelem(i2).gt.99)then
               write(filename(13:15),'(i3)')nucelem(i2)
             elseif(nucelem(i2).gt.9)then
               write(filename(14:15),'(i2)')nucelem(i2)
             else
               write(filename(15:15),'(i1)')nucelem(i2)
             endif

            inquire(file=filename,exist=lexist)
            if(.not.lexist) then
              write(*,*)'Error: file not found ',filename
              stop
            endif

            open(wunit,file=filename,form='unformatted',status='old')
            rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
            do i3=1,num_weightsp(paircount)
              read(wunit)weightsp(i3,paircount)
            enddo

            close(wunit)

          enddo ! i2
        enddo ! i1
      endif ! lshort
!!
!!
      return
      end
