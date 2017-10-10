!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90

      subroutine getvdwparams()
!!
      use fileunits
      use nnflags 
      use globaloptions
!!
      implicit none
!!
      integer i1
      integer ztemp1
      integer ztemp2
!!
      character*2 element1
      character*2 element2
      character*40 keyword              ! internal
!!
      logical lfound(npairs)
!!
      lfound(:)=.false.
!!
      if(nn_type_vdw.eq.1)then
        allocate(vdw_param(npairs,5))
!!
        open(nnunit,file='input.nn',form='formatted',status='old')
        rewind(nnunit)
 10     continue
        read(nnunit,*,END=20)keyword
        if(keyword.eq.'vdw_param')then
          backspace(nnunit)
          read(nnunit,*,ERR=99)keyword,element1,element2
          call nuccharge(element1,ztemp1)
          call nuccharge(element2,ztemp2)
          lfound(pairindex(ztemp1,ztemp2))=.true.
          backspace(nnunit)
          read(nnunit,*,ERR=99)keyword,element1,element2,&
          vdw_param(pairindex(ztemp1,ztemp2),1),&
          vdw_param(pairindex(ztemp1,ztemp2),2),&
          vdw_param(pairindex(ztemp1,ztemp2),3),&
          vdw_param(pairindex(ztemp1,ztemp2),4)
        endif !(keyword.eq.reference)then
        goto 10
!!
 20     continue
        close(nnunit)
      endif ! nn_type_vdw
!!
!! check if all parameters have been found
      do i1=1,npairs

!!
!! CHANGE ANDI: GFORTRAN: gfortran forbids .eq. for logicals, use .eqv. instead.
!!                        ifort allows both.
!!
       !if(lfound(i1).eq..false.)then
        if(lfound(i1).eqv..false.)then
!! END CHANGE

          write(ounit,*)'ERROR: not all vdW parameters have been found in getvdwparams'
          stop
        endif
      enddo

!! write parameters
      write(ounit,'(a)')' vdW parameters for Grimme scheme:'
      write(ounit,'(a)')'             Pair     ',&
        's6    C6ij   d   Rr   '
      do i1=1,npairs
        write(ounit,'(a,a2,x,a2,4f14.6)')' vdW parameters ',&
          element(elempair(i1,1)),element(elempair(i1,2)),&
          vdw_param(i1,1),vdw_param(i1,2),vdw_param(i1,3),&
          vdw_param(i1,4)
      enddo
!!
      return
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop
!!
      end
