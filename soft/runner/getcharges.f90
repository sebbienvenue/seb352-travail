!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - ewaldenergies_para.f90 
!! - optimize_ewald.f90
!!
!! very similar to geteshort.f90
!!
      subroutine getcharges(ndim1,npoints,&
        zelem_local,num_atoms_local,&
        symfunctione_local,nnatomcharge_local)
!!
      use fileunits
      use globaloptions
      use nnewald
!!
      implicit none
!!
      integer npoints
      integer ndim1                                                    ! in
      integer zelem_local(ndim1,max_num_atoms) 
      integer zelem(max_num_atoms)
      integer num_atoms_local(ndim1)
      integer i1
!! 
      real*8 symfunctione_local(maxnum_funcvalues_elec,max_num_atoms,ndim1)
      real*8 nnatomcharge_local(ndim1,max_num_atoms)                   ! out
      real*8 nnatomcharge(max_num_atoms)                               ! internal
!!
!!
      do i1=1,npoints
!!
!! calculate the atomic charges 
        zelem(:)=zelem_local(i1,:)
!!
        nnatomcharge(:)=0.0d0
        call calconecharge(num_atoms_local(i1),&
          zelem,symfunctione_local(1,1,i1),&
          nnatomcharge)
!!
        nnatomcharge_local(i1,:)=nnatomcharge(:)
!!
      enddo ! i1
!!
      return
      end
