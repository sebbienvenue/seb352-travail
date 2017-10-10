!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - geterror.f90
!! - precondition.f90
!! - fitting_batch.f90
!! - fittingpair.f90 
!! - geterrorpair.f90
!! - preconditionpair.f90
!!
      subroutine readfunctions(iswitch,unit,npoints,ndim,&
         max_num,maxnum_funcvalues_local,num_funcvalues_local,&
         symfunction_list_local)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!     
      integer unit                                                  ! internal
      integer ndim                                                  ! in
      integer iswitch                                               ! in
      integer npoints                                               ! in
      integer max_num                                               ! in
      integer maxnum_funcvalues_local                                     ! in
      integer num_funcvalues_local(ndim)                                  ! in
      integer i1,i2,i3                                              ! internal
!!
      real*8 symfunction_list_local(maxnum_funcvalues_local,max_num,nblock)     ! out
!!
      do i1=1,npoints
        if(iswitch.eq.1)then    ! atomic NN
          read(unit,*)num_atoms_list(i1)
          do i2=1,num_atoms_list(i1)
            read(unit,*)zelem_list(i1,i2),&
              (symfunction_list_local(i3,i2,i1),i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
          enddo ! i2
        elseif(iswitch.eq.2)then ! pair NN
          read(unit,*)num_atoms_list(i1),num_pairs_list(i1)
          do i2=1,num_pairs_list(i1)
            read(unit,*)zelemp_list(1,i1,i2),zelemp_list(2,i1,i2),&
              (symfunction_list_local(i3,i2,i1),i3=1,num_funcvalues_local(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))))
          enddo ! i2
        else
          write(ounit,*)'ERROR: unknown iswitch in readfunctions ',iswitch
          stop
        endif
        read(unit,*) totalcharge_list(i1),totalenergy_list(i1),&
          shortenergy_list(i1),elecenergy_list(i1) 
      enddo ! i1
!!
      return
      end
