!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: rescale dsfuncdxyz and strs arrays

!! multipurpose subroutine

!! called by: - prediction.f90
!!            - getdshortdw.f90
!!            - optimize_short_combined.f90
!!            - predictionpair.f90
!!
      subroutine scaledsfunc_para(natoms,atomindex,max_num_neighbors_local,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        ndim,minvalue_local,maxvalue_local,&
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
      integer zelem(max_num_atoms)                                      ! in
      integer i1,i2,i3,i4                                               ! internal
      integer natoms                                                    ! in
      integer atomindex(natoms)                                         ! in
      integer max_num_neighbors_local                                   ! in
!!
      real*8 dsfuncdxyz_local(maxnum_funcvalues_local,natoms,0:max_num_neighbors_local,3) ! in/out
      real*8 minvalue_local(ndim,maxnum_funcvalues_local)               ! in 
      real*8 maxvalue_local(ndim,maxnum_funcvalues_local)               ! in 
      real*8 strs_local(3,3,maxnum_funcvalues_local,natoms)             ! in/out
      real*8 scmin_local                                                ! in
      real*8 scmax_local                                                ! in
!!
!!
      do i1 = 1,3
        do i2 = 0,max_num_neighbors_local, 1
          do i3 = 1,natoms
            do i4 = 1,num_funcvalues_local(elementindex(zelem(atomindex(i3))))
!! scale each symmetry function derivative value for the respective element
              dsfuncdxyz_local(i4,i3,i2,i1)=dsfuncdxyz_local(i4,i3,i2,i1)/&
              (maxvalue_local(elementindex(zelem(atomindex(i3))),i4)&
              -minvalue_local(elementindex(zelem(atomindex(i3))),i4))&
              *(scmax_local-scmin_local)
            enddo ! i4
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
!! scale stress components, CHECK IF THIS IS RIGHT!!!
      if(ldostress)then
        do i2=1,natoms
          do i3=1,num_funcvalues_local(elementindex(zelem(atomindex(i2))))
            strs_local(:,:,i3,i2)=strs_local(:,:,i3,i2)/&
            (maxvalue_local(elementindex(zelem(atomindex(i2))),i3)&
            -minvalue_local(elementindex(zelem(atomindex(i2))),i3))&
            *(scmax_local-scmin_local)
          enddo ! i3
        enddo ! i2
      endif
!!
      return
      end

