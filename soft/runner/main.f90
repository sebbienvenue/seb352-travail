!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
      program RuNNer
!!
      use fileunits 
      use nnflags
      use globaloptions
      use mpi_mod
      use timings
      use runneripi
!!
      implicit none
!!
      integer iseed ! seed for random numbers
!!
!!======================================================================
!! prepare mpi
!!======================================================================
      call mpi_init(mpierror)
      if(mpierror.ne.0)then
        write(ounit,*)'Error in mpi_init ',mpierror
        stop
      endif 
!! get number of processes mpisize
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!!======================================================================
!! prepare timer
!!======================================================================
      call zerotime(runtimeday,runtimestart,runtimeend)
      call abstime(runtimestart,runtimeday)
!!
!!======================================================================
!! initialization
!!======================================================================
      call initnn(iseed)
!!
!!======================================================================
!! mode 1: Generate symmetry functions
!!======================================================================
      if(mode.eq.1)then
        call mode1(iseed)
!!
!!======================================================================
!! mode 2: Fitting 
!!======================================================================
      elseif(mode.eq.2)then
        call mode2(iseed)
!!
!!======================================================================
!! mode 3: Prediction 
!!======================================================================
      elseif(mode.eq.3)then
        if(luseipi)then
          call predict_ipi()
        else
          call predict()
        endif
!!
      else
        write(ounit,*)'Error: unknown RuNNer mode ',mode
        stop
      endif
!!
!!======================================================================
!! final cleanup, deallocations, print summary
!!======================================================================
      call abstime(runtimeend,runtimeday)
      call cleanup()
!!
!!======================================================================
!! shutdown mpi
!!======================================================================
      call mpi_finalize(mpierror)
      if(mpierror.ne.0)then
        write(ounit,*)'mpierror finalize ',mpierror
      endif
      end
