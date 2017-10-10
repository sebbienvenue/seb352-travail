!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90
!!
      subroutine readelementlayerspair(maxnodes_local,&
        maxnum_layers_local,num_layers_local,nodes_local,actfunc_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2
      integer ztemp1                                    ! internal
      integer ztemp2                                    ! internal
      integer maxnum_layers_local                       ! in 
      integer num_layers_local(npairs)                  ! out
      integer nodes_local(maxnum_layers_local,npairs)   ! out
      integer maxnodes_local                            ! in
      integer icount                                    ! internal
      integer jcount                                    ! internal
!!
      character*40 dummy                                ! internal
      character*2 elementtemp1                          ! internal
      character*2 elementtemp2                          ! internal
      character*1 actfunc_local(maxnodes_local,maxnum_layers_local,npairs) ! in/out

!!
      backspace(nnunit)
      read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2
      call checkelement(elementtemp1)
      call checkelement(elementtemp2)
      call nuccharge(elementtemp1,ztemp1)
      call nuccharge(elementtemp2,ztemp2)
      backspace(nnunit)
      icount=0
      jcount=0
      do i1=1,nelem
        do i2=i1,nelem
          jcount=jcount+1
            if((ztemp1.eq.elempair(jcount,1)).and.(ztemp2.eq.elempair(jcount,2)))then
              icount=jcount
            elseif((ztemp2.eq.elempair(jcount,1)).and.(ztemp1.eq.elempair(jcount,2)))then
              icount=jcount
            endif
          enddo ! i2
        enddo ! i1
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,num_layers_local(icount)
        num_layers_local(icount)=num_layers_local(icount)+1
        if(num_layers_local(icount).gt.maxnum_layers_local)then
          write(ounit,*)'Error: pair ',ztemp1,ztemp2,' has too many hidden layers'
          stop !'
        endif
!! set number of nodes in new output layer to 1
        nodes_local(num_layers_local(icount),icount)=1
!! delete activation functions for other output nodes
        do i1=2,maxnodes_local
          actfunc_local(i1,num_layers_local(icount),icount)=' '
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

