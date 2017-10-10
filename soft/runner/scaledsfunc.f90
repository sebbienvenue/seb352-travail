!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine
!! Purpose: scale dsfuncdxyz and strs arrays

!! called by: 
!! - getallelectrostatic.f90
!! - getallshortforces.f90
!! - optimize_short_combined.f90
!! - optimize_short_combinedpair.f90
!!
      subroutine scaledsfunc(max_num_neighbors_local,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        ndim,num_atoms,minvalue_local,maxvalue_local,&
        scmin_local,scmax_local,&
        zelem,dsfuncdxyz_local,strs_local)
!!
      use fileunits
      use globaloptions
!!  
      implicit none
!!
      integer maxnum_funcvalues_local                                   ! in
      integer num_funcvalues_local(ndim)                                ! in
      integer ndim                                                      ! in
      integer num_atoms                                                 ! in
      integer zelem(max_num_atoms)                                      ! in
      integer i0,i1,i2,i3                                               ! internal
      integer max_num_neighbors_local                                   ! in
!!
      real*8 dsfuncdxyz_local(maxnum_funcvalues_local,max_num_atoms,0:max_num_neighbors_local,3)   ! in/out
      real*8 minvalue_local(ndim,maxnum_funcvalues_local)                             ! in 
      real*8 maxvalue_local(ndim,maxnum_funcvalues_local)                             ! in 
      real*8 strs_local(3,3,maxnum_funcvalues_local,max_num_atoms)                     ! in/out
      real*8 scmin_local                                                   ! in
      real*8 scmax_local                                                   ! in
!!
!!
      do i0=1,3
        do i1=0,max_num_neighbors_local
          do i2=1,num_atoms
            do i3=1,num_funcvalues_local(elementindex(zelem(i2)))
!! scale each symmetry function derivative value for the respective element
          dsfuncdxyz_local(i3,i2,i1,i0)=dsfuncdxyz_local(i3,i2,i1,i0)/&
          (maxvalue_local(elementindex(zelem(i2)),i3)&
          -minvalue_local(elementindex(zelem(i2)),i3))&
          *(scmax_local-scmin_local)
            enddo ! i3
          enddo ! i2
        enddo ! i1
      enddo ! i0
!!
!! scale stress components, FIXME: CHECK IF THIS IS RIGHT!!!
      do i2=1,num_atoms
        do i3=1,num_funcvalues_local(elementindex(zelem(i2)))
          strs_local(:,:,i3,i2)=strs_local(:,:,i3,i2)/&
          (maxvalue_local(elementindex(zelem(i2)),i3)&
          -minvalue_local(elementindex(zelem(i2)),i3))&
          *(scmax_local-scmin_local)
        enddo ! i3
      enddo ! i2
!!
      return
      end
