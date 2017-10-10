!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - prediction.f90
!! - predictionpair.f90
!!
      subroutine calconecharge_para(natoms,atomindex,&
        zelem,symfunctione,nnatomcharge)
!!
      use fileunits
      use globaloptions
      use nnewald
!!
      implicit none
!!
      integer zelem(max_num_atoms)
      integer natoms                                       ! in
      integer atomindex(natoms)                            ! in
!!
      integer i1
!!
      real*8 symfunctione(maxnum_funcvalues_elec,natoms)       ! in
      real*8 nnoutput                                      ! internal
      real*8 nnatomcharge(max_num_atoms)                   ! out, CAUTION, just some elements are calculated here
      real*8 nodes_values_dummy(maxnum_layers_elec,maxnodes_elec) ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_elec,maxnodes_elec)    ! just dummy in this routine
!!
!!
!! loop over all atoms
      do i1=1,natoms
!!
        call calconenn(1,maxnum_funcvalues_elec,maxnodes_elec,&
          maxnum_layers_elec,num_layers_elec(elementindex(zelem(atomindex(i1)))),&
          maxnum_weights_elec,nodes_elec(0,elementindex(zelem(atomindex(i1)))),&
          symfunctione(1,i1),weights_elec(1,elementindex(zelem(atomindex(i1)))),&
          nodes_values_dummy,nodes_sum_dummy,&
          nnoutput,actfunc_elec)
!!
        nnatomcharge(atomindex(i1))=nnoutput
!!
      enddo ! i1
!!
      return
      end
