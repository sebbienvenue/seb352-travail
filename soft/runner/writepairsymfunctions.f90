!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calcpairfunctions.f90 
!!
      subroutine writepairsymfunctions(unit_local,i1,&
        pairs_charge_list,&
        maxnum_funcvalues_local,num_funcvalues_local,symfunctionp_local)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
      integer num_funcvalues_local(npairs)                                ! in
      integer unit_local                                                  ! in
      integer maxnum_funcvalues_local                                     ! in
      integer i1                                                          ! in
      integer i2,i3                                                       ! internal
      integer pairs_charge_list(2,listdim,nblock)                          ! in 
!!
      real*8 symfunctionp_local(maxnum_funcvalues_local,max_num_atoms,nblock)
!!
      write(unit_local,'(i6,x,i6)')num_atoms_list(i1),num_pairs_list(i1)
      do i2=1,num_pairs_list(i1)
!! check if number of digits in write format is sufficient
        do i3=1,num_funcvalues_local(pairindex(pairs_charge_list(1,i2,i1),pairs_charge_list(2,i2,i1)))
          if(symfunctionp_local(i3,i2,i1).ge.1000.d0)then
            write(ounit,*)'ERROR: symfunctionp is too large for write statement ',symfunctionp_local(i3,i2,i1)
            stop !'
          endif
          if(symfunctionp_local(i3,i2,i1).le.-100.d0)then
            write(ounit,*)'ERROR: symfunctionp is too small for write statement ',symfunctionp_local(i3,i2,i1)
            stop !'
          endif
        enddo
        write(unit_local,'(i3,x,i3,x,500f15.10)')pairs_charge_list(1,i2,i1),pairs_charge_list(2,i2,i1),&
          (symfunctionp_local(i3,i2,i1),i3=1,num_funcvalues_local(pairindex(pairs_charge_list(1,i2,i1),pairs_charge_list(2,i2,i1))))   
      enddo ! i2
      write(unit_local,'(4f20.10)')totalcharge_list(i1)/num_atoms_list(i1),&
        totalenergy_list(i1)/num_atoms_list(i1),&
        shortenergy_list(i1)/num_atoms_list(i1),&
        elecenergy_list(i1)/num_atoms_list(i1)
!! 
      return
!!
      end
