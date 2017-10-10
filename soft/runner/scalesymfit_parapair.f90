!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine scalesymfit_parapair(&
        nstruct,n_start,n_end,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)

!!
      use mpi_mod
      use fileunits
      use globaloptions
      use nnshort_pair
      use structures
!!
      implicit none
!!
      integer i1                          ! internal
      integer nstruct                     ! in
      integer n_start                     ! in
      integer n_end                       ! in
      integer icount                      ! internal
      integer, dimension(:), allocatable :: num_pairs_mpi  ! internal
      integer, dimension(:,:,:), allocatable :: zelemp_mpi    ! internal
!!
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)                         ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)                         ! in
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)                          ! in
      real*8, dimension(:,:,:), allocatable :: symfunctionp_mpi      ! internal
!!
!!
!! allocate mpi arrays for scalesym
        allocate(num_pairs_mpi(nstruct))
        allocate(zelemp_mpi(2,nstruct,max_num_pairs))
        allocate(symfunctionp_mpi(maxnum_funcvalues_short_pair,max_num_pairs,nstruct))
!!
!! copy symmetry functions to mpi array
        icount=0
        do i1=n_start,n_end
          icount=icount+1
          num_pairs_mpi(icount)       = num_pairs_list(i1)
          zelemp_mpi(:,icount,:)      = zelemp_list(:,i1,:)
          symfunctionp_mpi(:,:,icount)= symfunction_short_pair_list(:,:,i1)
        enddo
!!
!! scale symmetry functions for the short-range interaction
!! the normal version of scalesym here also works for the parallel code

       call scalesympair(nstruct,nstruct,&
           num_pairs_mpi,&
           zelemp_mpi,symfunctionp_mpi,&
           minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)
!!
!! copy mpi array back to full array
        icount=0
        symfunction_short_pair_list(:,:,:)=0.0d0
        do i1=n_start,n_end
          icount=icount+1
          symfunction_short_pair_list(:,:,i1) = symfunctionp_mpi(:,:,icount)
        enddo
!!
!! merge all scaled arrys
        call mpi_allreduce(mpi_in_place,symfunction_short_pair_list,&
          nblock*max_num_atoms*maxnum_funcvalues_short_pair,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! deallocate mpi arrays
        deallocate(num_pairs_mpi)
        deallocate(zelemp_mpi)
        deallocate(symfunctionp_mpi)
!!
!!
      return
      end
