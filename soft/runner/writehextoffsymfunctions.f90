!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calcfunctions.f90 
!!
      subroutine writehextoffsymfunctions(unit_local,i1,&
        maxnum_funcvalues_local,num_funcvalues_local,symfunction_local)
!!
      use fileunits
      use globaloptions
      use structures
      use basismod
!!
      implicit none
      integer num_funcvalues_local
      integer maxnum_funcvalues_local
      integer unit_local
      integer i1,i2,i3,i4,i5
      real*8, dimension(:), allocatable :: hextoff_out
!!
      real*8 symfunction_local(maxnum_funcvalues_local,nblock)
!!
!
!
!      write(*,*) 'CMH TEST num_atoms_list'
!      write(*,*) i1
!      write(*,*) num_atoms_list
!      write(*,*) num_atoms_list(i1)
!      write(*,*) totalcharge_list(i1)
!      write(*,*) num_atoms_list(i1)
!      write(*,*) totalenergy_list(i1)
!      write(*,*) shortenergy_list(i1)
!      write(*,*) elecenergy_list(i1)
!      PAUSE
!
!
!!      write(*,*) num_atoms_list(i1)
!!      write(unit_local,*) 'test'     
      i2 = tripletindex(hextoff_training_triplet(1),hextoff_training_triplet(2),hextoff_training_triplet(3))
!!      write(*,*) i1,i2, num_funcvalues_local(i2)
!!      write(*,*) symfunction_local
!!      do i2=1,num_atoms_list(i1)
!! check if number of digits in write format is sufficient
        do i3=1,num_funcvalues_local
          if(symfunction_local(i3,i1).ge.1000.d0)then
            write(ounit,*)'ERROR: symfunction is too large for write statement ',symfunction_local(i3,i1)
            stop !'
          endif
          if(symfunction_local(i3,i1).le.-100.d0)then
            write(ounit,*)'ERROR: symfunction is too small for write statement ',symfunction_local(i3,i1)
            stop !'
          endif
        enddo
!!        write(unit_local,*) 'test2'
        write(unit_local,'(x,500f15.10)')&
            (symfunction_local(i3,i1),i3=1,num_funcvalues_local)
        allocate(hextoff_out(num_basis(elementindex(hextoff_training_triplet(1)))*&
                             num_basis(elementindex(hextoff_training_triplet(2)))))
        i5 = 0
        do i3 = 1,num_basis(elementindex(hextoff_training_triplet(1)))
          do i4 = 1,num_basis(elementindex(hextoff_training_triplet(2)))
            i5 = i5 + 1
            hextoff_out(i5) = hextoff_list(i1,i3,i4)
!!            write(*,*) hextoff_out(i5)
          enddo
        enddo
        write(unit_local,'(x,500f15.10)') (hextoff_out(i3),i3=1,(num_basis(elementindex(hextoff_training_triplet(1)))*&
                             num_basis(elementindex(hextoff_training_triplet(2)))))       
        !!write(unit_local,'()')
!!      enddo ! i2
!!      write(unit_local,'(4f20.10)')totalcharge_list(i1)/num_atoms_list(i1),&
!!        totalenergy_list(i1)/num_atoms_list(i1),&
!!        shortenergy_list(i1)/num_atoms_list(i1),&
!!        elecenergy_list(i1)/num_atoms_list(i1)
      
      return
!!
      end
