!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine readsymfunction(ounit,iunit,nelem,&
        maxnum_funcvalues,maxnum_funcvaluese,&  
        num_funcvalues,num_funcvaluese,&  
        nucelem,&
        function_type,function_typee,&
        funccutoff,funccutoffe,maxcutoff,maxcutoffe,&
        rshift,rshifte,eta,etae,zeta,zetae,lambda,lambdae,&
        lshort,lewald,ldebug)
!!
      use mpi_mod
!!
      implicit none
!!
      integer i1,i2,i3,i4                          ! internal
      integer ounit                                ! in
      integer iunit                                ! in
      integer nelem                                ! in
      integer maxnum_funcvalues                    ! in
      integer maxnum_funcvaluese                   ! in
      integer num_funcvalues(nelem)                ! in
      integer num_funcvaluese(nelem)               ! in
      integer function_type(num_functions)         ! out 
      integer function_typee(num_functionse)       ! out 
      integer nucelem(nelem)                       ! in
!!
      real*8 funccutoff(num_functions)             ! out 
      real*8 funccutoffe(num_functionse)           ! out 
      real*8 maxcutoff                             ! out
      real*8 maxcutoffe                            ! out
      real*8 rshift(num_functions)                 ! out 
      real*8 rshifte(num_functionse)               ! out 
      real*8 eta(num_functions)                    ! out 
      real*8 etae(num_functionse)                  ! out 
      real*8 zeta(num_functions)                   ! out 
      real*8 zetae(num_functionse)                 ! out 
      real*8 lambda(num_functions)                 ! out 
      real*8 lambdae(num_functionse)               ! out 
!!
      logical lshort              ! in
      logical lewald              ! in
      logical ldebug              ! in
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!! read the symmetry functions just for process id 0
      if(mpirank.eq.0)then
!! read information on short range symmetry functions '
      if(lshort)then
        open(iunit,file='input.sym',form='formatted',status='old')
        rewind(iunit)
        do i1=1,num_functions
          read(iunit,*)function_type(i1)
          backspace(iunit)
          if(function_type(i1).eq.1)then ! pure cutoff function
            read(iunit,*)function_type(i1),funccutoff(i1)
          elseif(function_type(i1).eq.2)then ! radial function
            read(iunit,*)function_type(i1),eta(i1),rshift(i1),funccutoff(i1)
          elseif(function_type(i1).eq.3)then ! angular function
            read(iunit,*)function_type(i1),eta(i1),lambda(i1),zeta(i1),funccutoff(i1)
          elseif(function_type(i1).eq.4)then ! radial function
            read(iunit,*)function_type(i1),eta(i1),funccutoff(i1)
          elseif(function_type(i1).eq.5)then ! Just Cartesian Coordinate 
            read(iunit,*)function_type(i1),eta(i1)
          elseif(function_type(i1).eq.6)then ! Just Cartesian Coordinate 
            read(iunit,*)function_type(i1),funccutoff(i1)
          else
            write(ounit,*)'Error: Function not implemented ',function_type(i1)
            stop
          endif
          maxcutoff=max(maxcutoff,funccutoff(i1))
        enddo ! i1
        close(iunit)
      endif ! lshort
!!
!! read information on electrostatic symmetry functions
      if(lewald)then
        open(iunit,file='inpute.sym',form='formatted',status='old')
        rewind(iunit)
        do i1=1,num_functionse
          read(iunit,*)function_typee(i1)
          backspace(iunit)
          if(function_typee(i1).eq.1)then ! pure cutoff function
            read(iunit,*)function_typee(i1),funccutoffe(i1)
          elseif(function_typee(i1).eq.2)then ! radial function
            read(iunit,*)function_typee(i1),etae(i1),rshifte(i1),funccutoffe(i1)
          elseif(function_typee(i1).eq.3)then ! angular function
            read(iunit,*)function_typee(i1),etae(i1),lambdae(i1),zetae(i1),funccutoffe(i1)
          elseif(function_typee(i1).eq.4)then ! radial function
            read(iunit,*)function_typee(i1),etae(i1),funccutoffe(i1)
          elseif(function_typee(i1).eq.5)then ! just Cartesian Coordinate
            read(iunit,*)function_typee(i1),etae(i1)
          elseif(function_typee(i1).eq.5)then ! just Cartesian Coordinate
            read(iunit,*)function_typee(i1),funccutoffe(i1)
          else
            write(ounit,*)'Error: Function not implemented ',function_typee(i1)
            stop
          endif
          maxcutoffe=max(maxcutoffe,funccutoffe(i1))
        enddo ! i1
        close(iunit)
      endif ! lewald
      endif ! mpirank.eq.0
!!
!! distribute symmetry function data to all nodes
      if(lshort)then
        call mpi_bcast(function_type,num_functions,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff,num_functions,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta,num_functions,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda,num_functions,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta,num_functions,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift,num_functions,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxcutoff,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lewald)then
        call mpi_bcast(function_typee,num_functionse,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoffe,num_functionse,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(etae,num_functionse,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambdae,num_functionse,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zetae,num_functionse,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshifte,num_functionse,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxcutoffe,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!! calculate nnindex and nnindexe just for process id 0
      if(mpirank.eq.0)then
!! short range NN
        if(lshort)then
          call getnnindex(ounit,num_functions,&
            nelem,nucelem,function_type,nnindex,ldebug)
        endif ! lshort
!! electrostatic NN
        if(lewald)then
          call getnnindex(ounit,num_functionse,&
            nelem,nucelem,function_typee,nnindexe,ldebug)
        endif ! lewald
      endif ! mpirank.eq.0
!!
!! distribute the nnindex and nnindexe arrays
      if(lshort)then
        call mpi_bcast(nnindex,num_functions*nelem*nelem,&
          mpi_integer,0,mpi_comm_world,mpierror)
      endif ! lshort
      if(lewald)then
        call mpi_bcast(nnindexe,num_functionse*nelem*nelem,&
          mpi_integer,0,mpi_comm_world,mpierror)
      endif ! lewald
!!
      return
      end
