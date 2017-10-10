!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - geteshort.f90
!!
      subroutine calconeshort(num_atoms,&
        zelem,symfunction,&
        nneshort,nnatomenergy)
!!
      use fileunits
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer num_atoms                 ! in
      integer zelem(max_num_atoms)
!!
      integer i1
!!
      real*8 symfunction(maxnum_funcvalues_short_atomic,max_num_atoms)    ! in
      real*8 symfunction_atom(maxnum_funcvalues_short_atomic)             ! internal
      real*8 weights(maxnum_weights_short_atomic)                    ! internal
!! CAUTION: nnoutput assumes just one output node here
      real*8 nnoutput                                     ! internal
      real*8 nneshort                                       ! out
      real*8 nodes_values_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic) ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic)    ! just dummy in this routine
      real*8 nnatomenergy(max_num_atoms)                  ! out
!!
!!
      nneshort=0.0d0
      nnatomenergy(:)=0.0d0
!!
!!#################################################################
!! serial original:
!!
      do i1=1,num_atoms
!!
!!
        symfunction_atom(:)=symfunction(:,i1)
        weights(:)=weights_short_atomic(:,elementindex(zelem(i1)))
!!
!! calculate nodes_values and nnoutput for atom i1
        call calconenn(1,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic(elementindex(zelem(i1))),&
          maxnum_weights_short_atomic,nodes_short_atomic(0,elementindex(zelem(i1))),&
          symfunction_atom,weights,nodes_values_dummy,nodes_sum_dummy,&
          nnoutput,actfunc_short_atomic(1,1,elementindex(zelem(i1))))
!!
        nnatomenergy(i1)=nnoutput
        nneshort=nneshort+nnoutput
!!
      enddo ! i1
!!
!!
      return
      end
