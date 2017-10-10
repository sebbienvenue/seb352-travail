!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - prediction.f90
!! - predictionpair.f90
!!
      subroutine scalesym_para(natoms,atomindex,&
         ndim1,ndim2,npoints,&
         maxnum_funcvalues_local,num_funcvalues_local,&
         zelem_local,symfunction_local,&
         minvalue_local,maxvalue_local,avvalue_local,&
         scmin_local,scmax_local)
!!
      use fileunits
      use globaloptions 
!!
      implicit none
!!
      integer ndim1                                                   ! in    ! number of NNs
      integer ndim2                                                   ! in   ! number of structures
      integer npoints                                                 ! in
      integer maxnum_funcvalues_local                                 ! in 
      integer num_funcvalues_local(ndim1)                             ! in
      integer zelem_local(ndim2,max_num_atoms)                        ! in 
      integer i1,i2,i3                                                ! internal
      integer natoms                                                  ! in
      integer atomindex(natoms)                                       ! in
!!
      real*8 symfunction_local(maxnum_funcvalues_local,natoms,ndim2)  ! in/out
      real*8 minvalue_local(ndim1,maxnum_funcvalues_local)            ! in
      real*8 maxvalue_local(ndim1,maxnum_funcvalues_local)            ! in
      real*8 avvalue_local(ndim1,maxnum_funcvalues_local)             ! in
      real*8 scmin_local                                              ! in
      real*8 scmax_local                                              ! in
!!
!!
!!
      do i1=1,npoints
        do i2=1,natoms
          do i3=1,num_funcvalues_local(elementindex(zelem_local(i1,atomindex(i2))))
            if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
              symfunction_local(i3,i2,i1)=symfunction_local(i3,i2,i1) &
              -avvalue_local(elementindex(zelem_local(i1,atomindex(i2))),i3)
            elseif(lscalesym.and..not.lcentersym)then
!! Scale each symmetry function value for the respective element
              symfunction_local(i3,i2,i1)=&
             (symfunction_local(i3,i2,i1)&
             -minvalue_local(elementindex(zelem_local(i1,atomindex(i2))),i3))/ &
             (maxvalue_local(elementindex(zelem_local(i1,atomindex(i2))),i3)-&
              minvalue_local(elementindex(zelem_local(i1,atomindex(i2))),i3))&
              *(scmax_local-scmin_local) + scmin_local
            elseif(lscalesym.and.lcentersym)then
              symfunction_local(i3,i2,i1)=&
             (symfunction_local(i3,i2,i1)&
             -avvalue_local(elementindex(zelem_local(i1,atomindex(i2))),i3))&
             / &
             (maxvalue_local(elementindex(zelem_local(i1,atomindex(i2))),i3)-&
              minvalue_local(elementindex(zelem_local(i1,atomindex(i2))),i3))
            else
            endif
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
!!
!!
      return
      end
