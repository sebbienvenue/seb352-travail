!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate nneshort_mpi and nnshortforces_mpi for nstruct structures

!! called by:
!! - geterror.f90
!! - precondition.f90
!!
      subroutine getshortenergies_para(nstruct,ndone,&
        nenergies,imaxerroreshort,&
        num_atoms_mpi,zelem_mpi,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,symfunction_mpi,lattice_mpi,&
        nnshortforce_mpi,xyzstruct_mpi,nneshort_mpi,&
        rmse_short,mad_short,maxerrorshort,shortenergy_mpi,&
        lperiodic_mpi)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
      use symfunctions
      use nnshort_atomic
!!
      implicit none
!!
      integer nstruct                                    ! in
      integer ndone                                      ! in
      integer num_atoms_mpi(nstruct)                     ! in
      integer zelem_mpi(nstruct,max_num_atoms)           ! in
      integer nenergies                                  ! in/out
      integer imaxerroreshort                            ! in/out
!!
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)              ! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)              ! in
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)               ! in
      real*8 symfunction_mpi(maxnum_funcvalues_short_atomic,max_num_atoms,nstruct)     ! in/out
      real*8 nneshort_mpi(nstruct)                       ! out
      real*8 nnshortforce_mpi(3,max_num_atoms,nstruct)   ! out
      real*8 lattice_mpi(3,3,nstruct)                    ! in
      real*8 xyzstruct_mpi(3,max_num_atoms,nstruct)      ! in
      real*8 rmse_short                                  ! in/out
      real*8 mad_short                                   ! in/out
      real*8 shortenergy_mpi(nstruct)                    ! in 
!! symmetry function parameters
      real*8 maxerrorshort                               ! in/out
!!
      logical lperiodic_mpi(nstruct)                     ! in
!!
!!
!! scale symmetry functions for the short-range interaction
        call scalesym(nelem,nstruct,nstruct,&
          maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,num_atoms_mpi,&
          zelem_mpi,symfunction_mpi,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
          scmin_short_atomic,scmax_short_atomic)
!!
!! predict the short range NN output for npoint data sets
        call geteshort(nstruct,nstruct,&
          zelem_mpi,num_atoms_mpi,&
          symfunction_mpi,nneshort_mpi)
!!
!! predict the short range NN forces for the test points here
!! change T.M. : not called for preconditioner to save time
        if(luseforces.and.(.not.lprecond))then
          call getallshortforces(nstruct,nstruct,&
            num_atoms_mpi,zelem_mpi,&
            symfunction_mpi,nnshortforce_mpi,&
            lattice_mpi,xyzstruct_mpi,minvalue_short_atomic,maxvalue_short_atomic,&
            lperiodic_mpi)
        endif ! luseforces
!!
!! calculate rmse_short: in/out rmse_short
        call calcrmse_energy(nstruct,nstruct,nenergies,&
          ndone,imaxerroreshort,&
          rmse_short,mad_short,maxerrorshort,&
          maxenergy,shortenergy_mpi,nneshort_mpi)
!!

      return
      end
