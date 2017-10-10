!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!!            - prediction.f90
!!            - getdshortdw.f90
!!            - optimize_short_combined.f90
!!
      subroutine calconeshort_para(istruct,natoms,atomindex,&
        zelem,symfunction,&
        nnatomenergy,nntotalenergy)
!!
      use mpi_mod
      use fileunits
      use globaloptions
      use nnshort_atomic
      use saturation
!!
      implicit none
!!
      integer zelem(max_num_atoms)
      integer natoms                                      ! in
      integer atomindex(natoms)                           ! in
      integer istruct                                     ! in
!!
      integer i1,i2,i3
!!
      real*8 symfunction(maxnum_funcvalues_short_atomic,natoms)           ! in
!! CAUTION: just one output node is assumed here
      real*8 nnoutput       ! internal
      real*8 nodes_values_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic) ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic)    ! just dummy in this routine
      real*8 nnatomenergy(max_num_atoms)                  ! out
      real*8 nntotalenergy                                ! out
!!
!!
      do i1=1,natoms
!!
!! calculate nnoutput for atom i1
        if(ldetect_saturation)then
          saturation_element=elementindex(zelem(atomindex(i1)))
        endif
        call calconenn(1,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic(elementindex(zelem(atomindex(i1)))),&
          maxnum_weights_short_atomic,nodes_short_atomic(0,elementindex(zelem(atomindex(i1)))),&
          symfunction(1,i1),weights_short_atomic(1,elementindex(zelem(atomindex(i1)))),nodes_values_dummy,nodes_sum_dummy,&
          nnoutput,actfunc_short_atomic(1,1,elementindex(zelem(atomindex(i1)))))
        if(ldetect_saturation)then
          saturation_element=0 ! to avoid calculation of nodes_saturation in other parts of RuNNer
        endif
!!
        if(ldebug)then
!!          write(ounit,'(a,i6,a)')'nodes values for atom ',i1,' :'
          write(ounit,'(a)')'            struct  atom layer  node           nodes_sum        nodes_values'
          do i2=1,num_layers_short_atomic(elementindex(zelem(atomindex(i1)))) !'
            do i3=1,nodes_short_atomic(i2,elementindex(zelem(atomindex(i1))))
              write(ounit,'(a12,4i6,2f20.10)')'NODES SHORT ',&
                istruct,i1,i2,i3,nodes_sum_dummy(i2,i3),nodes_values_dummy(i2,i3)
            enddo ! i3
          enddo ! i2
          write(ounit,*)'-------------------------------------------------------------'
        endif ! ldebug
!!
        nnatomenergy(atomindex(i1))=nnoutput
        nntotalenergy = nntotalenergy+nnoutput
!!
      enddo ! i1
!!
      return
      end
