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
!! This subroutine checks the symmetry functions for extrapolation
!!
      subroutine checkextrapolation(natoms,atomindex,ndim,&
             maxnum_funcvalues_local,num_funcvalues_local,&
             zelem,symfunction_local,&
             minvalue_local,maxvalue_local,slabel,lextrapolation)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer maxnum_funcvalues_local            ! in
      integer num_funcvalues_local(ndim)        ! in
      integer ndim                        ! in
      integer zelem(max_num_atoms)         ! in     
      integer i1,i3                     ! internal
      integer natoms                       ! in
      integer atomindex(natoms)            ! in
!!
      real*8 symfunction_local(maxnum_funcvalues_local,natoms) ! in
      real*8 minvalue_local(ndim,maxnum_funcvalues_local)            ! in
      real*8 maxvalue_local(ndim,maxnum_funcvalues_local)            ! in
      real*8 threshold                                 ! internal
!!
      character*5 slabel                    ! in
!!
      logical lextrapolation
!!      
      threshold=1.0d-8
      lextrapolation=.false.
!!
      do i1=1,natoms
        do i3=1,num_funcvalues_local(elementindex(zelem(atomindex(i1))))
          if((symfunction_local(i3,i1)-maxvalue_local(elementindex(zelem(atomindex(i1))),i3)).gt.threshold)then
            write(ounit,'(a,a5,x,2i5,x,a,2f18.8)')'### EXTRAPOLATION WARNING ### ',slabel,&
            atomindex(i1),i3,'too large ',symfunction_local(i3,i1),&
            maxvalue_local(elementindex(zelem(atomindex(i1))),i3)
            lextrapolation=.true.
          elseif((-symfunction_local(i3,i1)+minvalue_local(elementindex(zelem(atomindex(i1))),i3)).gt.threshold)then
            write(ounit,'(a,a5,x,2i5,x,a,2f18.8)')'### EXTRAPOLATION WARNING ### ',slabel,&
            atomindex(i1),i3,'too small ',symfunction_local(i3,i1),&
            minvalue_local(elementindex(zelem(atomindex(i1))),i3)
            lextrapolation=.true.
          endif
        enddo ! i3
      enddo ! i1
!!
      return
      end
