!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - geterrorpair.f90
!! - preconditionpair.f90
!!
      subroutine getshortenergies_parapair(nstruct,ndone,&
        nenergies,imaxerroreshort,&
        num_atoms_mpi,num_pairs_mpi,&
        zelem_mpi,zelemp_mpi,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        symfunctionp_mpi,lattice_mpi,&
        nnshortforce_mpi,xyzstruct_mpi,nneshort_mpi,&
        rmse_short,mad_short,maxerrorshort,shortenergy_mpi,&
        lperiodic_mpi)
!!
      use fileunits
      use fittingoptions
      use globaloptions
      use symfunctions
      use nnshort_pair
!!
      implicit none
!!
      integer nstruct                                    ! in
      integer ndone                                      ! in
      integer nenergies                                  ! out
      integer num_atoms_mpi(nstruct)                     ! in
      integer num_pairs_mpi(nstruct)                     ! in
      integer zelem_mpi(nstruct,max_num_atoms)           ! in
      integer zelemp_mpi(2,nstruct,max_num_pairs)        ! in
      integer imaxerroreshort                            ! in/out
!!
      real*8 symfunctionp_mpi(maxnum_funcvalues_short_pair,max_num_pairs,nstruct) ! in/out
      real*8 nneshort_mpi(nstruct)                       ! out
      real*8 nnshortforce_mpi(3,max_num_atoms,nstruct)     ! out
      real*8 lattice_mpi(3,3,nstruct)                    ! in
      real*8 xyzstruct_mpi(3,max_num_atoms,nstruct)      ! in
      real*8 rmse_short                                  ! out
      real*8 mad_short                                   ! out
      real*8 shortenergy_mpi(nstruct)                    ! in 
!! symmetry function parameters
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)        ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)        ! in
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)         ! in
      real*8 maxerrorshort                               ! in/out
!!
      logical lperiodic_mpi(nstruct)                     ! in
!!
!!
!! scale symmetry functions for the short-range interaction
        call scalesympair(nstruct,nstruct,&
          num_pairs_mpi,zelemp_mpi,symfunctionp_mpi,&
          minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)
!!
!! predict the short range NN output for npoint data sets
        call geteshortpair(nstruct,&
          zelemp_mpi,num_atoms_mpi,num_pairs_mpi,&
          symfunctionp_mpi,nneshort_mpi)
!!
!! predict the short range NN forces for the test points here
        if(luseforces)then
          call getallshortforcespair(nstruct,nstruct,&
            num_atoms_mpi,zelem_mpi,zelemp_mpi,&
            symfunctionp_mpi,nnshortforce_mpi,lattice_mpi,xyzstruct_mpi,&
            minvalue_short_pair,maxvalue_short_pair,lperiodic_mpi)
        endif ! luseforces
!!
!! calculate rmse_short: in/out rmse_short
        call calcrmse_energy(nstruct,nstruct,nenergies,&
             ndone,imaxerroreshort,&
             rmse_short,mad_short,maxerrorshort,&
             maxenergy,shortenergy_mpi,nneshort_mpi)
!!
!!
      return
      end
