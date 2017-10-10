!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! CAUTION: THIS SUBROUTINE CANNOT WORK, because neighbor_para fills lst arrays always from beginning

!! called by:
!!
      subroutine getneighborsatomic_para(n_start,n_end,natoms,&
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
      integer n_start                                        ! in
      integer n_end                                          ! in
      integer natoms                                         ! in
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
!! get neighbor lists lsta,lstb,lstc and lste for atoms n_start to n_end
!!=======================================================
      call neighbor_para(n_start,n_end,&
        num_atoms,zelem,lsta,lstb,lstc,lste,&
        maxcutoff_local,lattice,xyzstruct,lperiodic)
!!
!!=======================================================
!! combine neighbor lists from all processes
!!=======================================================
      call mpi_allreduce(mpi_in_place,lsta,2*max_num_atoms,&
        mpi_integer,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,lstb,listdim*4,&
        mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,lstc,listdim,&
        mpi_integer,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,lste,listdim,&
        mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
!!=======================================================
!! if max number of neighbors is given in input.nn use that value
!!=======================================================
      if(lenforcemaxnumneighborsatomic)then
        max_num_neighbors_atomic_local=max_num_neighbors_atomic_input
      endif
!!
!!=======================================================
!! determine number of neighbors for atoms n_start to n_end
!!=======================================================
      do i1=n_start,n_end
        num_neighbors_atomic_local(i1)=lsta(2,i1)-lsta(1,i1)+1
      enddo
!!
!!=======================================================
!! combine num_neighbors_atomic from all processes
!!=======================================================
      call mpi_allreduce(mpi_in_place,num_neighbors_atomic_local,num_atoms,&
        mpi_integer,mpi_sum,mpi_comm_world,mpierror)
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
