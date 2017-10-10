!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate nneshort_mpi and nnshortforces_mpi for nstruct structures

!! called by:
!!
      subroutine getatomicoutput_para(nstruct,ndone,&
        nenergies,ncharges,&
        imaxerror_eshort,imaxerror_elec,imaxerror_etot,&
        num_atoms_mpi,zelem_mpi,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        symfunction_mpi,lattice_mpi,xyzstruct_mpi,&
        rmse_short,rmse_elec,rmse_charge,rmse_totalcharge,&
        mad_short,mad_elec,mad_charge,mad_totalcharge,&
        maxerror_eshort,maxerror_elec,maxerror_etot,&
        shortenergy_mpi,elecenergy_mpi,atomcharge_mpi,totalcharge_mpi,&
        nneshort_mpi,nnshortforce_mpi,nnatomcharge_mpi,nnchargesum_mpi,&
        nnelec_mpi,nnelecforce_mpi,&
        lperiodic_mpi)
!!
      use fileunits
      use fittingoptions
      use nnflags 
      use globaloptions
      use symfunctions
      use nnshort_atomic
!!
      implicit none
!!
      integer num_atoms_mpi(nstruct)                     ! in
      integer zelem_mpi(nstruct,max_num_atoms)           ! in
      integer nenergies                                  ! in/out
      integer ncharges                                   ! in/out
      integer nstruct                                    ! in
      integer ndone                                      ! in
      integer imaxerror_eshort                           ! in/out
      integer imaxerror_elec                             ! in/out
      integer imaxerror_etot                             ! in/out
      integer ndummy                                     ! internal
!! in
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)           ! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)           ! in
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)            ! in
      real*8 symfunction_mpi(maxnum_funcvalues_short_atomic,max_num_atoms,nstruct) ! in/out
      real*8 lattice_mpi(3,3,nstruct)                    ! in
      real*8 xyzstruct_mpi(3,max_num_atoms,nstruct)      ! in
      real*8 shortenergy_mpi(nstruct)                    ! in 
      real*8 elecenergy_mpi(nstruct)                     ! in 
      real*8 totalcharge_mpi(nstruct)                    ! in
      real*8 atomcharge_mpi(nstruct,max_num_atoms)       ! in
!! output
      real*8 nneshort_mpi(nstruct)                       ! out
      real*8 nnshortforce_mpi(3,max_num_atoms,nstruct)   ! out
      real*8 nnatomcharge_mpi(nstruct,max_num_atoms)     ! out
      real*8 nnchargesum_mpi(nstruct)                    ! out
      real*8 nnelec_mpi(nstruct)                         ! out
      real*8 nnelecforce_mpi(3,max_num_atoms,nstruct)    ! out
!! errors
      real*8 rmse_short                                  ! in/out
      real*8 rmse_elec                                   ! in/out
      real*8 rmse_charge                                 ! in/out
      real*8 rmse_totalcharge                            ! in/out
      real*8 mad_short                                   ! in/out
      real*8 mad_elec                                    ! in/out
      real*8 mad_charge                                  ! in/out
      real*8 mad_totalcharge                             ! in/out
!! symmetry function parameters
      real*8 maxerror_eshort                             ! in/out
      real*8 maxerror_elec                               ! in/out
      real*8 maxerror_etot                               ! in/out
      real*8 edummy                                      ! internal
!!
      logical lperiodic_mpi(nstruct)                     ! in
!!
!!
!!========================================================
!! scale symmetry functions 
!!========================================================
      call scalesym(nelem,nstruct,nstruct,&
        maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
        num_atoms_mpi,zelem_mpi,symfunction_mpi,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        scmin_short_atomic,scmax_short_atomic)
!!
!!========================================================
!! predict the NN output for npoint data sets
!!========================================================
      call getatom(nstruct,nstruct,&
        zelem_mpi,num_atoms_mpi,&
        symfunction_mpi,nneshort_mpi,&
        nnatomcharge_mpi,nnchargesum_mpi)
!!
!!========================================================
!! predict the short and electrostatic forces and electrostatic energy 
!!========================================================
      call getallforces(nstruct,nstruct,&
        num_atoms_mpi,zelem_mpi,&
        lattice_mpi,xyzstruct_mpi,&
        symfunction_mpi,&
        minvalue_short_atomic,maxvalue_short_atomic,&
        nnshortforce_mpi,nnelecforce_mpi,&
        nnelec_mpi,&
        lperiodic_mpi)
!!
!!========================================================
!! calculate rmse_short: in/out rmse_short
!!========================================================
      call calcrmse_energy(nstruct,nstruct,nenergies,&
        ndone,imaxerror_eshort,&
        rmse_short,mad_short,maxerror_eshort,&
        maxenergy,shortenergy_mpi,nneshort_mpi)
!!
      if(lelec.and.(nn_type_elec.eq.2))then
!!========================================================
!! calculate the RMSE for the atomic charges and the total charge
!!========================================================
        call calcrmse_charge(nstruct,nstruct,ncharges,&
          zelem_mpi,num_atoms_mpi,rmse_charge,mad_charge,&
          totalcharge_mpi,rmse_totalcharge,mad_totalcharge,&
          atomcharge_mpi,nnatomcharge_mpi,nnchargesum_mpi)
!!
!!========================================================
!! calculate the RMSE for the electrostatic energy: in/out rmse_ewald 
!!========================================================
        call calcrmse_energy(nstruct,nstruct,ndummy,&
          ndone,imaxerror_elec,&
          rmse_elec,mad_elec,maxerror_elec,&
          edummy,elecenergy_mpi,nnelec_mpi)
      endif ! lelec
!!
      return
      end
