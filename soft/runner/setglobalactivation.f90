!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine
!!
!! called by:
!! - readinput.f90
!!
      subroutine setglobalactivation(ndim,counter,&
        maxnum_layers_local,maxnodes_local,nodes_local,actfunc_local,&
        actfunc_localdummy,keyword)
!!
      use fileunits
!!
      implicit none
!!
      integer i,i1,i2,i3
      integer ndim                                                  ! in
      integer maxnum_layers_local                                   ! in
      integer maxnodes_local                                        ! in
      integer nodes_local(0:maxnum_layers_local,ndim)               ! in
      integer counter                                               ! in/out
!!
      character*1 actfunc_localdummy(maxnum_layers_local)                 ! out 
      character*1 actfunc_local(maxnodes_local,maxnum_layers_local,ndim)  ! out
      character*40 dummy                                                  ! internal
      character*40 keyword                              ! in
!! 
      counter=counter+1
!!
      backspace(nnunit)
      read(nnunit,*,ERR=99)dummy,(actfunc_localdummy(i),i=1,maxnum_layers_local)
      do i1=1,maxnum_layers_local
        do i3=1,ndim
          do i2=1,nodes_local(i1,i3)
            actfunc_local(i2,i1,i3)=actfunc_localdummy(i1)
          enddo ! i2
!! fill empty part of actfunc_local
          if(nodes_local(i1,i3).lt.maxnodes_local)then
            do i2=nodes_local(i1,i3)+1,maxnodes_local
              actfunc_local(i2,i1,i3)=' '
            enddo
          endif
        enddo ! i3
      enddo ! i1
!!
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop

      end
