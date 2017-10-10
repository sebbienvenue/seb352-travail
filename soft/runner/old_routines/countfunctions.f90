!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - initialization.f90
!!
      subroutine countfunctions(iunit,ounit,num_functions,num_funcvalues,&
                 nelem,ldebug)
!!
      implicit none
!!
      integer iunit            ! in
      integer ounit            ! in
      integer num_functions    ! out
      integer num_funcvalues   ! out
      integer function_type    ! internal
      integer nelem            ! in
      integer i1,i2
!!
      logical ldebug
!!
!! initialization
      num_functions = 0
      num_funcvalues = 0
!!
 10   continue
      read(iunit,*,err=20,end=20)function_type
      num_functions=num_functions+1
!!
!! determine the number of structure function values including cross terms
      if(function_type.eq.1)then ! radial function
        num_funcvalues=num_funcvalues+nelem
      elseif(function_type.eq.2)then ! radial function
        num_funcvalues=num_funcvalues+nelem
      elseif(function_type.eq.3)then ! angular function
        num_funcvalues=num_funcvalues+nelem ! 2 neighbors of the same element
!! add cross terms
        if(nelem.gt.1)then
          do i1=1,nelem-1
            num_funcvalues=num_funcvalues+i1    ! 2 neighbors of different element
          enddo
        endif
      elseif(function_type.eq.4)then ! radial function
        num_funcvalues=num_funcvalues+nelem
      elseif(function_type.eq.5)then ! just Cartesian coordinates 
        num_funcvalues=num_funcvalues+1
      elseif(function_type.eq.6)then ! just Cartesian coordinates 
        num_funcvalues=num_funcvalues+1
      else
        write(ounit,*)'Error: Undefined structure function ',function_type
        stop
      endif
!!
!!      if(ldebug)then
!!        write(*,*)'type ',function_type
!!        write(*,*)'total values ',num_funcvalues
!!      endif
!!
      goto 10
 20   continue
!!
!!      if(ldebug)then
!!        write(*,*)'Number of functions found ',num_functions
!!      endif
!!
      return
      end
