!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!!
      subroutine getatom_para(natoms,atomindex,&
        zelem,symfunction,&
        nnatomenergy,nnatomcharge)
!!
      use mpi_mod
      use fileunits
      use nnflags 
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer zelem(max_num_atoms)
      integer natoms                                      ! in
      integer atomindex(natoms)                           ! in
      integer nndim                                       ! internal
!!
      integer i1,i2,i3
!!
      real*8 symfunction(maxnum_funcvalues_short_atomic,natoms)                   ! in
      real*8, dimension(:)  , allocatable :: nnoutput                             ! internal
      real*8 nodes_values_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic) ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic)    ! just dummy in this routine
      real*8 nnatomenergy(max_num_atoms)                                          ! out
      real*8 nnatomcharge(max_num_atoms)                                          ! out
!!
!!============================================
!! determine number of output nodes
!!============================================
      if(lelec.and.(nn_type_elec.eq.2))then
        nndim=2
      else
        nndim=1
      endif
      allocate(nnoutput(nndim))
!!
      do i1=1,natoms
!!
!! calculate nnoutput for atom i1
        call calconenn(1,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic(elementindex(zelem(atomindex(i1)))),&
          maxnum_weights_short_atomic,nodes_short_atomic(0,elementindex(zelem(atomindex(i1)))),&
          symfunction(1,i1),weights_short_atomic(1,elementindex(zelem(atomindex(i1)))),nodes_values_dummy,nodes_sum_dummy,&
          nnoutput,actfunc_short_atomic(1,1,elementindex(zelem(atomindex(i1)))))
!!
        nnatomenergy(atomindex(i1))=nnoutput(1)
        if(lelec.and.(nn_type_elec.eq.2))then
          nnatomcharge(atomindex(i1))=nnoutput(2)
        endif
!!
      enddo ! i1
!!
      deallocate(nnoutput)
!!
      return
      end
