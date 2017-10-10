!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!!
      subroutine getneighborsatomic(&
        num_atoms,num_neighbors_atomic_local,zelem,&
        max_num_neighbors_atomic_local,&
        lsta,lstc,lste,&
        maxcutoff_local,lattice,xyzstruct,lstb,lperiodic)
!!
      use mpi_mod
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1
      integer icount
      integer max_num_neighbors_atomic_local                 ! out
      integer num_neighbors_atomic_local(num_atoms)          ! out 
      integer num_atoms                                      ! in
      integer zelem(max_num_atoms)                           ! in
      integer lsta(2,max_num_atoms)                          ! out, numbers of neighbors
      integer lstc(listdim)                                  ! out, identification of atom
      integer lste(listdim)                                  ! out, nuclear charge of atom
!!
      real*8 lstb(listdim,4)                                 ! out, xyz and r_ij 
      real*8 maxcutoff_local                                 ! in
      real*8 lattice(3,3)                                    ! in
      real*8 xyzstruct(3,max_num_atoms)                      ! in
!!
      logical lperiodic                                      ! in
!!
!!=======================================================
!! initializations
!!=======================================================
      lsta(:,:)          =0
      lstb(:,:)          =0.0d0
      lstc(:)            =0
      lste(:)            =0
      num_neighbors_atomic_local(:) =0
      max_num_neighbors_atomic_local=0
!!
!!=======================================================
!! get neighbor lists lsta,lstb,lstc and lste 
!!=======================================================
      call neighbor(num_atoms,zelem,&
        lsta,lstb,lstc,lste,&
        maxcutoff_local,lattice,xyzstruct,lperiodic)
!!
!!=======================================================
!! if max number of neighbors is given in input.nn use that value
!!=======================================================
      if(lenforcemaxnumneighborsatomic)then
        max_num_neighbors_atomic_local=max_num_neighbors_atomic_input
      endif
!!
!!=======================================================
!! determine number of neighbors for all atoms
!!=======================================================
      do i1=1,num_atoms
        num_neighbors_atomic_local(i1)=lsta(2,i1)-lsta(1,i1)+1
      enddo
!!
!!=======================================================
!! determine max_num_neighbors_atomic_local if not given in input.nn
!!=======================================================
      do i1=1,num_atoms
        if(.not.lenforcemaxnumneighborsatomic)then
          max_num_neighbors_atomic_local=max(max_num_neighbors_atomic_local,lsta(2,i1)-lsta(1,i1)+1)
        else
!!        Just for checking. Could be removed to save some time
          if(max_num_neighbors_atomic_local.lt.(lsta(2,i1)-lsta(1,i1)))then
            write(ounit,*)'ERROR: max_num_neighbors_atomic in input is too small'
            stop
          endif
        endif
      enddo
!!
      return
      end
