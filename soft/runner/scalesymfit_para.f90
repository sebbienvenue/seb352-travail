!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - fitting.f90
!! - fitting_batch.f90
!! - fittingpair.f90
!! - getsensitivity.f90
!!
      subroutine scalesymfit_para(iswitch,ndim,nstruct,n_start,n_end,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        minvalue_local,maxvalue_local,avvalue_local,symfunction_list_local,&
        scmin_local,scmax_local)
!!
      use mpi_mod
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!
      integer iswitch                     ! in
      integer i1                          ! internal
      integer ndim                        ! in
      integer nstruct                     ! in
      integer n_start                     ! in
      integer n_end                       ! in
      integer icount                      ! internal
      integer maxnum_funcvalues_local           ! in
      integer num_funcvalues_local(ndim)       ! in
      integer, dimension(:), allocatable :: num_atoms_mpi  ! internal
      integer, dimension(:), allocatable :: num_pairs_mpi  ! internal
      integer, dimension(:,:), allocatable :: zelem_mpi    ! internal
      integer, dimension(:,:,:), allocatable :: zelemp_mpi    ! internal
!!
      real*8 minvalue_local(ndim,maxnum_funcvalues_local)                         ! in
      real*8 maxvalue_local(ndim,maxnum_funcvalues_local)                         ! in
      real*8 avvalue_local(ndim,maxnum_funcvalues_local)                          ! in
      real*8 symfunction_list_local(maxnum_funcvalues_local,max_num_atoms,nblock)  ! in
      real*8, dimension(:,:,:), allocatable :: symfunction_mpi         ! internal
      real*8 scmin_local                                               ! in
      real*8 scmax_local                                               ! in
!!
!!
!! allocate mpi arrays for scalesym
        if((iswitch.eq.1).or.(iswitch.eq.3))then
          allocate(zelem_mpi(nstruct,max_num_atoms))
          allocate(num_atoms_mpi(nstruct))
          allocate(symfunction_mpi(maxnum_funcvalues_local,max_num_atoms,nstruct))
        elseif(iswitch.eq.2)then
          allocate(zelemp_mpi(2,nstruct,max_num_pairs))
          allocate(num_pairs_mpi(nstruct))
          allocate(symfunction_mpi(maxnum_funcvalues_local,max_num_pairs,nstruct))
        endif
!!
!! copy symmetry functions to mpi array
        icount=0
        do i1=n_start,n_end
          icount=icount+1
          symfunction_mpi(:,:,icount) = symfunction_list_local(:,:,i1)
          if((iswitch.eq.1).or.(iswitch.eq.3))then
            num_atoms_mpi(icount)       = num_atoms_list(i1)
            zelem_mpi(icount,:)         = zelem_list(i1,:)
          elseif(iswitch.eq.2)then
            num_pairs_mpi(icount)       = num_pairs_list(i1)
            zelemp_mpi(:,icount,:)      = zelemp_list(:,i1,:)
          endif
        enddo
!!
!! scale symmetry functions for the short-range interaction
!! the normal version of scalesym here also works for the parallel code
        if((iswitch.eq.1).or.(iswitch.eq.3))then
          call scalesym(ndim,nstruct,nstruct,&
            maxnum_funcvalues_local,num_funcvalues_local,num_atoms_mpi,&
            zelem_mpi,symfunction_mpi,&
            minvalue_local,maxvalue_local,avvalue_local,&
            scmin_local,scmax_local)
        elseif(iswitch.eq.2)then
          call scalesympair(ndim,nstruct,&
            num_pairs_mpi,zelemp_mpi,symfunction_mpi,&
            minvalue_local,maxvalue_local,avvalue_local)
        endif
!!
!! copy mpi array back to full array
        icount=0
        symfunction_list_local(:,:,:)=0.0d0
        do i1=n_start,n_end
          icount=icount+1
          symfunction_list_local(:,:,i1) = symfunction_mpi(:,:,icount)
        enddo
!!
!! merge all scaled arrys
        if((iswitch.eq.1).or.(iswitch.eq.3))then
          call mpi_allreduce(mpi_in_place,symfunction_list_local,&
            nblock*max_num_atoms*maxnum_funcvalues_local,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        elseif(iswitch.eq.2)then
          call mpi_allreduce(mpi_in_place,symfunction_list_local,&
            nblock*max_num_pairs*maxnum_funcvalues_local,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif
!!
!! deallocate mpi arrays
        if(allocated(num_atoms_mpi))deallocate(num_atoms_mpi)
        if(allocated(num_pairs_mpi))deallocate(num_pairs_mpi)
        if(allocated(zelem_mpi))deallocate(zelem_mpi)
        if(allocated(zelemp_mpi))deallocate(zelemp_mpi)
        deallocate(symfunction_mpi)
!!
!!
      return
      end
