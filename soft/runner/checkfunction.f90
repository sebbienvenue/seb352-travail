!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine checkfunction()
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer num,nume
      integer idummy
      integer num_atoms
      integer num_pairs
      integer i1

      real*8 rdummy

      num=0
      nume=0

!! count structures in function.data
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)
!!
 11   continue
      if(nn_type_short.eq.1)then
        read(symunit,*,END=10)num_atoms
        do i1=1,num_atoms
          read(symunit,*)idummy
        enddo
      elseif(nn_type_short.eq.2)then
        read(symunit,*,END=10)num_atoms,num_pairs
        do i1=1,num_pairs
          read(symunit,*)idummy
        enddo
      else
        write(ounit,*)'ERROR in checkfunction, unknown nn_type_short'
        stop
      endif
      read(symunit,*)rdummy
      num=num+1
      goto 11
 10   continue
!!
      close(symunit)
!!
!! count structures in function.data
      open(symeunit,file='functione.data',form='formatted',status='old')
      rewind(symeunit)
!!
 21   continue
      read(symeunit,*,END=20)num_atoms
      do i1=1,num_atoms
        read(symeunit,*)idummy
      enddo
      read(symeunit,*)rdummy
      nume=nume+1
      goto 21
 20   continue
!!
      close(symeunit)
!!
      if(num.ne.nume)then
        write(ounit,*)'ERROR: You tried to fit Eshort and Charges simultaneously,'
        write(ounit,*)'but the numbers of structures are different in '
        write(ounit,*)'function.data and functione.data: ',num,nume
        stop
      endif
!!
      return
      end
