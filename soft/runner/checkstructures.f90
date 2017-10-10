!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - initialization.f90
!!
      subroutine checkstructures(ielem,lelement)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer ielem                   ! out
      integer i1
!!
      logical lelement(102)           ! out
!!
!! initialization
      lelement(:)=.false.
!!
      open(dataunit,file='input.data',form='formatted')
      rewind(dataunit)
!!
!! loop over all structures in input.data
!! CAUTION: in case of mode 2 totnum_structures is not determined before, because we don't want to use input.data in mode 2
!! => ielem is always 0 after this subroutine
      do i1=1,totnum_structures
        call checkonestructure(i1,lelement)
      enddo ! i1
!!
      if(ldebug)then
!!      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'============================================================='
      endif !'
      close(dataunit)
!!
!! determine the number of elements
!! Caution: not all elements need to be present in each structure
      ielem=0
      do i1=1,102
        if(lelement(i1)) ielem=ielem+1
      enddo
!!
      return
      end
