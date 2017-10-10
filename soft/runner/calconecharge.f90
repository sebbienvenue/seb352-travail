!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getcharges.f90
!! - calcfunctions.f90
!! - calcpairfunctions.f90
!!
      subroutine calconecharge(num_atoms,&
        zelem,symfunctione,nnatomcharge)
!!
      use fileunits
      use globaloptions
      use nnewald
!!
      implicit none
!!
      integer num_atoms
      integer zelem(max_num_atoms)
!!
      integer i1
!!
      real*8 symfunctione(maxnum_funcvalues_elec,max_num_atoms)
      real*8 nnoutput
      real*8 nnatomcharge(max_num_atoms)                  ! out
      real*8 nodes_values_dummy(maxnum_layers_elec,maxnodes_elec) ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_elec,maxnodes_elec)    ! just dummy in this routine
!!
!!
!! loop over all atoms
      do i1=1,num_atoms
!!
          call calconenn(1,maxnum_funcvalues_elec,maxnodes_elec,&
            maxnum_layers_elec,num_layers_elec(elementindex(zelem(i1))),&
            maxnum_weights_elec,nodes_elec(0,elementindex(zelem(i1))),&
            symfunctione(1,i1),weights_elec(1,elementindex(zelem(i1))),&
            nodes_values_dummy,nodes_sum_dummy,&
            nnoutput,actfunc_elec(1,1,elementindex(zelem(i1))))
!!
          nnatomcharge(i1)=nnoutput
!!
      enddo ! i1
!!
      return
      end
