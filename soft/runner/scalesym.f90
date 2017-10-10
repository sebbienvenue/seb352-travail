!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! This is a multipurpose subroutine to scale general symfunctions

!! called by:
!! - getshortenergies_para.f90
!! - ewaldenergies_para.f90
!! - scalesymfit_para.f90 
!! - precondition.f90
!!
      subroutine scalesym(ndim1,ndim2,npoints,&
         maxnum_funcvalues_local,num_funcvalues_local,num_atoms_local,&
         zelem_local,symfunction_local,&
         minvalue_local,maxvalue_local,avvalue_local,&
         scmin_local,scmax_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer ndim1                                                          ! in    ! number of NNs
      integer ndim2                                                          ! in   ! number of structures
      integer npoints                                                        ! in
      integer maxnum_funcvalues_local                                        ! in
      integer num_funcvalues_local(ndim1)                                    ! in
      integer num_atoms_local(ndim2)                                         ! in
      integer zelem_local(ndim2,max_num_atoms)                               ! in
      integer i1,i2,i3                                                       ! internal
      integer itemp                                                          ! internal
!!
      real*8 symfunction_local(maxnum_funcvalues_local,max_num_atoms,ndim2)  ! in/out
      real*8 minvalue_local(ndim1,maxnum_funcvalues_local)                   ! in
      real*8 maxvalue_local(ndim1,maxnum_funcvalues_local)                   ! in
      real*8 avvalue_local(ndim1,maxnum_funcvalues_local)                    ! in
      real*8 scmin_local                                                     ! in
      real*8 scmax_local                                                     ! in
!!
      do i1=1,npoints
        do i2=1,num_atoms_local(i1)
          itemp=elementindex(zelem_local(i1,i2))
          do i3=1,num_funcvalues_local(itemp)
            if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
              symfunction_local(i3,i2,i1)=symfunction_local(i3,i2,i1) &
              -avvalue_local(itemp,i3)
            elseif(lscalesym.and..not.lcentersym)then
!! Scale each symmetry function value for the respective element
              symfunction_local(i3,i2,i1)=&
             (symfunction_local(i3,i2,i1)-minvalue_local(itemp,i3))/ &
             (maxvalue_local(itemp,i3)-minvalue_local(itemp,i3))&
              *(scmax_local-scmin_local) + scmin_local
            elseif(lscalesym.and.lcentersym)then
              symfunction_local(i3,i2,i1)=&
             (symfunction_local(i3,i2,i1)-avvalue_local(itemp,i3))&
             / (maxvalue_local(itemp,i3)-minvalue_local(itemp,i3))
            else
            endif
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
      return
      end
