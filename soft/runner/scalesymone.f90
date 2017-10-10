!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - calcfunctions.f90
!! - calcpairfunctions.f90
!!
      subroutine scalesymone(ndim,&
         maxnum_funcvalues_local,num_funcvalues_local,num_atoms,&
         zelem,symfunction_local,&
         minvalue_local,maxvalue_local,avvalue_local,&
         scmin_local,scmax_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer ndim
      integer maxnum_funcvalues_local
      integer num_funcvalues_local(ndim)
      integer num_atoms
      integer zelem(max_num_atoms)
      integer i2,i3
!!
      real*8 symfunction_local(maxnum_funcvalues_local,max_num_atoms)  ! in/out
      real*8 minvalue_local(ndim,maxnum_funcvalues_local)
      real*8 maxvalue_local(ndim,maxnum_funcvalues_local)
      real*8 avvalue_local(ndim,maxnum_funcvalues_local)
      real*8 scmin_local                                               ! in
      real*8 scmax_local                                               ! in
!!
!!
        do i2=1,num_atoms
          do i3=1,num_funcvalues_local(elementindex(zelem(i2)))
            if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
              symfunction_local(i3,i2)=symfunction_local(i3,i2) &
              -avvalue_local(elementindex(zelem(i2)),i3)
            elseif(lscalesym.and..not.lcentersym)then
              symfunction_local(i3,i2)=&
             (symfunction_local(i3,i2)&
             -minvalue_local(elementindex(zelem(i2)),i3))/ &
             (maxvalue_local(elementindex(zelem(i2)),i3)-&
              minvalue_local(elementindex(zelem(i2)),i3))&
              *(scmax_local-scmin_local) + scmin_local
            elseif(lscalesym.and.lcentersym)then
              symfunction_local(i3,i2)=&
             (symfunction_local(i3,i2)-avvalue_local(elementindex(zelem(i2)),i3))&
             / &
             (maxvalue_local(elementindex(zelem(i2)),i3)-&
              minvalue_local(elementindex(zelem(i2)),i3))
            else
            endif
          enddo ! i3
        enddo ! i2
!!
      return
      end
