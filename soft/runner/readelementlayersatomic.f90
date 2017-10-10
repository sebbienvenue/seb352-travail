!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90
!!
      subroutine readelementlayersatomic(maxnodes_local,&
        maxnum_layers_local,num_layers_local,nodes_local,actfunc_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1
      integer ztemp                                     ! internal
      integer maxnum_layers_local                       ! in 
      integer num_layers_local(nelem)                   ! out
      integer nodes_local(maxnum_layers_local,nelem)    ! out
      integer maxnodes_local                            ! in
!!
      character*40 dummy                                ! internal
      character*2 elementtemp                           ! internal
      character*1 actfunc_local(maxnodes_local,maxnum_layers_local,nelem) ! in/out

!!
      backspace(nnunit)
      read(nnunit,*,ERR=99)dummy,elementtemp
      call checkelement(elementtemp)
      call nuccharge(elementtemp,ztemp)
      backspace(nnunit)
      read(nnunit,*,ERR=99)dummy,elementtemp,num_layers_local(elementindex(ztemp))
      num_layers_local(elementindex(ztemp))=num_layers_local(elementindex(ztemp))+1
      if(num_layers_local(elementindex(ztemp)).gt.maxnum_layers_local)then
        write(ounit,*)'Error: element ',ztemp,' has too many hidden layers'
        stop !'
      endif
!! set number of nodes in new output layer to 1
      nodes_local(num_layers_local(elementindex(ztemp)),elementindex(ztemp))=1
!! delete activation functions for other output nodes
      do i1=2,maxnodes_local
        actfunc_local(i1,num_layers_local(elementindex(ztemp)),elementindex(ztemp))=' '
      enddo
!!
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
 99   continue
      write(ounit,*)'Error: keyword ',dummy
      write(ounit,*)'is missing arguments '
      stop

      end
