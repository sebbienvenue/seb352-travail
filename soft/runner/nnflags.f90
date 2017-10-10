      module nnflags 
      implicit none

!! mode = runner mode (1=calculate symfunctions, 2=fitting, 3=prediction)
      integer mode

!! nn_type_short: type of neural network: nn_type_short=1 => atom-based NN, nn_type_short=2 => pair-based NN
      integer nn_type_short
      integer nn_type_nntb
!! determine if separate charge NN is used (nn_type_elec=1) 
!! or if charges are additional NN output (nn_type_elec=2)
!! of if fixed charges are used instead of NN output (nn_type_elec=3)
      integer nn_type_elec
      integer nn_type_vdw

      logical lshort
!! calculate electrostatic energy (and forces, stress etc.) for given charges (according to keyword nn_type_elec) 
      logical lelec
      logical lnntb
      integer originatom_id
      integer zatom_id

      contains
      subroutine distribute_nnflags()
      use mpi_mod
      implicit none
!!
      call mpi_bcast(mode,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nn_type_short,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nn_type_nntb,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nn_type_elec,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nn_type_vdw,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(mode,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lelec,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lnntb,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lshort,1,mpi_logical,0,mpi_comm_world,mpierror)

      end subroutine distribute_nnflags
      end module nnflags 
