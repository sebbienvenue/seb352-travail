module wrapper
    implicit none

    logical, allocatable :: extrapolated_atoms(:)
    logical :: bond_underrun

    ! Very ugly, but the only way to get it working is to use global variables.
    real(8), allocatable :: minvalue_short_atomic(:,:), maxvalue_short_atomic(:,:), avvalue_short_atomic(:,:)
    real(8) :: eshortmin, eshortmax
contains
    ! initwrapper sets everything up as close as the original RuNNer would do
    subroutine initwrapper(pred_max_num_atoms, ielem)
        use globaloptions, only : max_num_atoms, maxnum_funcvalues_short_atomic, nelem
        use mpi_mod, only : mpierror, mpirank, mpisize
        use mpi, only : mpi_comm_rank, mpi_comm_size, mpi_comm_world

        implicit none
        integer, intent(in) :: pred_max_num_atoms
        integer, intent(in) :: ielem

        integer :: iseed ! cannot use a constant as it is written at some place

        ! Calling initmode3 needs quiet a lot parameters. Fake some here to get the call right.
        ! num_pairs will be computed in the pair case (so not for me) and then needs: num_atoms, xyzstruct, lattice, lperiodic
        integer :: num_atoms_dummy, num_pairs_dummy, zelem_dummy(1)
        real(8) :: xyzstruct_dummy(1,1), lattice_dummy(1,1)
        logical :: lperiodic_dummy
        ! {max,min,av}value_{short_pair,elec} show up, but are never referenced due to wrong mode/flag/whatever
        real(8) :: minvalue_short_pair_dummy(1,1), maxvalue_short_pair_dummy(1,1), avvalue_short_pair_dummy(1,1)
        real(8) :: minvalue_elec_dummy(1,1), maxvalue_elec_dummy(1,1), avvalue_elec_dummy(1,1)
        ! Charge is never used in the short ranged part, but it's needed for initialization of the general mode.
        real(8) :: chargemin_dummy, chargemax_dummy

        ! mpi init as soon as possible
        call mpi_comm_rank(mpi_comm_world, mpirank, mpierror)
        call mpi_comm_size(mpi_comm_world, mpisize, mpierror)

        ! we don't have an input.data to get this value from, so the user has to provide it
        max_num_atoms = pred_max_num_atoms
        call initnn(iseed, ielem)

        ! now initialize the only arrays we will ever need
        allocate(minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic), &
                 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic), &
                 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic))
        allocate(extrapolated_atoms(max_num_atoms))

        ! initialize mode3
        call initmode3(num_atoms_dummy, num_pairs_dummy, zelem_dummy, minvalue_short_atomic, maxvalue_short_atomic, &
                       avvalue_short_atomic, minvalue_short_pair_dummy, maxvalue_short_pair_dummy, avvalue_short_pair_dummy, &
                       minvalue_elec_dummy, maxvalue_elec_dummy, avvalue_elec_dummy, eshortmin, eshortmax, chargemin_dummy, &
                       chargemax_dummy, lattice_dummy, xyzstruct_dummy, lperiodic_dummy)
    end subroutine initwrapper

    subroutine getenergyandforces(num_atoms, xyz, zelem, lattice, lperiodic, energy, forces, atomic_energies, ex_atoms, min_bond_underrun)
        use globaloptions, only : elementindex, lremoveatomenergies, max_num_atoms, maxnum_funcvalues_short_atomic, nelem
        use mpi_mod, only : mpierror
        use mpi, only : mpi_comm_world, mpi_logical, mpi_lor
        implicit none

        ! interface parameters
        integer, intent(in) :: num_atoms
        real(8), intent(in) :: xyz(3,num_atoms)
        integer, intent(in) :: zelem(num_atoms)
        real(8), intent(in) :: lattice(3,3)
        logical, intent(in) :: lperiodic
        real(8), intent(out) :: energy, forces(3,num_atoms), atomic_energies(num_atoms)
        logical, intent(out) :: ex_atoms(num_atoms)
        logical, intent(out) :: min_bond_underrun

        ! Some more dummy arguments for the call to predictionshortatomic.
        real(8) :: nnshortenergy_dummy, nnstress_short_dummy(3,3), &
                   sens_dummy(nelem,maxnum_funcvalues_short_atomic)
        ! Parameters we're really interested in.
        integer :: num_atoms_element(nelem)
        real(8) :: atomenergysum, nnatomenergy(max_num_atoms), nntotalenergy, nnshortforce(3,max_num_atoms)

        ! internal variables
        integer :: i

        extrapolated_atoms = .False.
        bond_underrun = .False.
        call predictionshortatomic(num_atoms, num_atoms_element, zelem, lattice, xyz, minvalue_short_atomic, &
                                   maxvalue_short_atomic, avvalue_short_atomic, eshortmin, eshortmax, nntotalenergy, nnshortforce, &
                                   nnatomenergy, nnshortenergy_dummy, nnstress_short_dummy, atomenergysum, sens_dummy, &
                                   lperiodic)

        call mpi_reduce(extrapolated_atoms, ex_atoms, num_atoms, mpi_logical, mpi_lor, 0, mpi_comm_world, mpierror)
        call mpi_reduce(bond_underrun, min_bond_underrun, 1, mpi_logical, mpi_lor, 0, mpi_comm_world, mpierror)
        forces = nnshortforce(:,:num_atoms)
        atomic_energies = nnatomenergy(:num_atoms)

        if (lremoveatomenergies) then
            num_atoms_element = 0
            ! fill num_atoms_element as this should be done by reading the structure which is skipped by the wrapper
            do i = 1, num_atoms, 1
                num_atoms_element(elementindex(zelem(i))) = num_atoms_element(elementindex(zelem(i))) + 1
            end do
            call addatoms(num_atoms, zelem, num_atoms_element, atomenergysum, nnatomenergy)
            energy = nntotalenergy + atomenergysum
        else
            energy = nntotalenergy
        end if
    end subroutine getenergyandforces
end module wrapper

! vim: set et ft=fortran sw=4 ts=4 tw=132:
