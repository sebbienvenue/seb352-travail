!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - geterror.f90
!! - geterrorpair.f90
!!
      subroutine ewaldenergies_para(nstruct,ndone,&
        imaxerror_elec,&
        num_atoms_mpi,zelem_mpi,ncharges,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        symfunctione_mpi,lattice_mpi,nnelecforce_mpi,xyzstruct_mpi,&
        nnewald_mpi,&
        rmse_charge,rmse_totalcharge,rmse_ewald,&
        mad_charge,mad_totalcharge,mad_ewald,maxerrorewald,&
        elecenergy_mpi,atomcharge_mpi,nnatomcharge_mpi,&
        totalcharge_mpi,nnchargesum_mpi,&
        lperiodic_mpi)
!!
      use fileunits
      use globaloptions 
      use symfunctions
      use nnewald
!!
      implicit none
!!
      integer nstruct                                    ! in
      integer ndone                                      ! in
      integer num_atoms_mpi(nstruct)                     ! in
      integer zelem_mpi(nstruct,max_num_atoms)           ! in
      integer ncharges                                   ! out
      integer ndummy                                     ! internal
      integer imaxerror_elec                             ! in/out
!!
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec) ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec) ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)  ! in
      real*8 symfunctione_mpi(maxnum_funcvalues_elec,max_num_atoms,nstruct)   ! in/out
      real*8 nnewald_mpi(nstruct)                        ! out
      real*8 nnelecforce_mpi(3,max_num_atoms,nstruct)    ! out
      real*8 lattice_mpi(3,3,nstruct)                    ! in
      real*8 xyzstruct_mpi(3,max_num_atoms,nstruct)      ! in
      real*8 rmse_charge                                 ! out
      real*8 rmse_totalcharge                            ! out
      real*8 rmse_ewald                                  ! out
      real*8 mad_charge                                  ! out
      real*8 mad_totalcharge                             ! out
      real*8 mad_ewald                                   ! out
      real*8 elecenergy_mpi(nstruct)                    ! in
      real*8 nnatomcharge_mpi(nstruct,max_num_atoms)     ! out
      real*8 totalcharge_mpi(nstruct)                    ! in
      real*8 atomcharge_mpi(nstruct,max_num_atoms)       ! in
      real*8 nnchargesum_mpi(nstruct)                    ! out
!! symmetry function parameters
      real*8 edummy                                      ! internal
      real*8 maxerrorewald                               ! in/out
!!
      logical lperiodic_mpi(nstruct)                     ! in
!!
!! initializations
      edummy = 1.d12
      ndummy = 0
!!
!! scale the symmetry functions for the charge prediction
        call scalesym(nelem,nstruct,nstruct,&
          maxnum_funcvalues_elec,num_funcvalues_elec,num_atoms_mpi,&
          zelem_mpi,symfunctione_mpi,&
          minvalue_elec,maxvalue_elec,avvalue_elec,&
          scmin_elec,scmax_elec)
!!
!! calculate the charges on the atoms
        call getcharges(nstruct,nstruct,&
          zelem_mpi,num_atoms_mpi,&
          symfunctione_mpi,&
          nnatomcharge_mpi)
!!
!! calculate the RMSE for the atomic charges and the total charge
        call calcrmse_charge(nstruct,nstruct,ncharges,&
          zelem_mpi,num_atoms_mpi,rmse_charge,mad_charge,&
          totalcharge_mpi,rmse_totalcharge,mad_totalcharge,&
          atomcharge_mpi,nnatomcharge_mpi,nnchargesum_mpi)
!!
!! calculate the electrostatic energy for nstruct points
        call getallelectrostatic(nstruct,nstruct,&
          num_atoms_mpi,&
          zelem_mpi,lattice_mpi,minvalue_elec,maxvalue_elec,&
          nnatomcharge_mpi,xyzstruct_mpi,nnewald_mpi,&
          symfunctione_mpi,nnelecforce_mpi,&
          lperiodic_mpi)
!!
!! calculate the RMSE for the electrostatic energy: in/out rmse_ewald 
        call calcrmse_energy(nstruct,nstruct,ndummy,&
          ndone,imaxerror_elec,&
          rmse_ewald,mad_ewald,maxerrorewald,&
          edummy,elecenergy_mpi,nnewald_mpi)
!!
      return
      end
