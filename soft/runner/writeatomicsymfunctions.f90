!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calcfunctions.f90 
!!
      subroutine writeatomicsymfunctions(unit_local,i1,&
        maxnum_funcvalues_local,num_funcvalues_local,symfunction_local)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
      integer num_funcvalues_local(nelem)
      integer maxnum_funcvalues_local
      integer unit_local
      integer i1,i2,i3
!!
      real*8 symfunction_local(maxnum_funcvalues_local,max_num_atoms,nblock)
!!
!
      write(unit_local,'(i6)')num_atoms_list(i1)
      do i2=1,num_atoms_list(i1)
!! check if number of digits in write format is sufficient
        do i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2)))
          if(symfunction_local(i3,i2,i1).ge.1000.d0)then
            write(ounit,*)'ERROR: symfunction is too large for write statement ',symfunction_local(i3,i2,i1)
            stop !'
          endif
          if(symfunction_local(i3,i2,i1).le.-100.d0)then
            write(ounit,*)'ERROR: symfunction is too small for write statement ',symfunction_local(i3,i2,i1)
            stop !'
          endif
        enddo
        write(unit_local,'(i3,x,1500f15.10)')zelem_list(i1,i2),&
          (symfunction_local(i3,i2,i1),i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
      enddo ! i2
      write(unit_local,'(4f20.10)')totalcharge_list(i1)/num_atoms_list(i1),&
        totalenergy_list(i1)/num_atoms_list(i1),&
        shortenergy_list(i1)/num_atoms_list(i1),&
        elecenergy_list(i1)/num_atoms_list(i1)
!! 
      return
!!
      end
