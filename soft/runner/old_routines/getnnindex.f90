!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getsymmetryfunctions
!! - getallshortforces
!! - prediction
!! - getallelectrostatic
!! - fitforcesshort
!!
      subroutine getnnindex(ounit,num_functions,&
        nelem,nucelem,function_type,nnindex,ldebug)
!!
      implicit none
!!
      integer ounit                                 ! in
      integer num_functions                         ! in
      integer nelem                                 ! in
      integer nucelem(nelem)                        ! in
      integer nnindex(num_functions,nelem,nelem)    ! out
      integer function_type(num_functions)          ! in
      integer i1,i2,i3
      integer ncount
!!
      logical ldebug
!!
      ncount=1
      do i1=1,num_functions
        if(function_type(i1).eq.1)then ! radial function
          do i2=1,nelem
            nnindex(i1,i2,:)=ncount
            ncount=ncount+1
          enddo
        elseif(function_type(i1).eq.2)then
          do i2=1,nelem
            nnindex(i1,i2,:)=ncount
            ncount=ncount+1
          enddo
        elseif(function_type(i1).eq.3)then
          do i2=1,nelem
            do i3=1,nelem
              if(nucelem(i3).ge.nucelem(i2))then
                nnindex(i1,i2,i3)=ncount
                nnindex(i1,i3,i2)=ncount
                ncount=ncount+1
              endif
            enddo
          enddo
        elseif(function_type(i1).eq.4)then
          do i2=1,nelem
            nnindex(i1,i2,:)=ncount
            ncount=ncount+1
          enddo
        elseif(function_type(i1).eq.5)then
          nnindex(i1,:,:)=ncount
          ncount=ncount+1
        elseif(function_type(i1).eq.6)then
          nnindex(i1,:,:)=ncount
          ncount=ncount+1
        else
          write(*,*)'function type not implemented ',i1,function_type(i1)
          stop
        endif
!!
      enddo ! i1
!!
!!      if(ldebug)then
!!        write(ounit,*)'nnindex array'
!!        do i1=1,num_functions
!!        do i2=1,nelem
!!          do i3=1,nelem
!!            write(ounit,'(3i5,4x,i5)')i1,i2,i3,nnindex(i1,i2,i3)
!!          enddo
!!        enddo
!!        enddo
!!        write(ounit,*)'-------------------------------------------------------------'
!!      endif
!!
      return
      end
