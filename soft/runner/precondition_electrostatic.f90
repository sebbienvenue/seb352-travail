!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting_electrostatic.f90 
!!
      subroutine precondition_electrostatic(ntrain,trainelem,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        avcharge,stddevcharge)
!!
      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions 
      use fittingoptions
      use symfunctions
      use structures
      use nnewald
!!
      implicit none
!!
      integer i1,i2         ! internal
      integer ntrain        ! in
      integer ncount        ! internal
      integer npoints       ! internal
      integer nstruct       ! internal
      integer n_start       ! internal
      integer n_end         ! internal
      integer, dimension(:), allocatable :: num_atoms_mpi  ! internal
      integer, dimension(:,:), allocatable :: zelem_mpi    ! internal
      integer icount          ! internal
      integer ncharge(nelem)                 ! internal
      integer trainelem(nelem)               ! internal
!!
      real*8 prefactor                       ! internal
      real*8 stddevcharge(nelem)             ! in
      real*8 avcharge(nelem)                 ! in
      real*8 nnstddevq(nelem)                ! internal
      real*8 rmse_dummy                      ! dummy
      real*8 mad_dummy                       ! dummy
      real*8 maxerror_dummy                  ! dummy
      real*8, dimension(:), allocatable :: elecenergy_mpi       ! internal
      real*8, dimension(:,:,:), allocatable :: xyzstruct_mpi     ! internal
      real*8, dimension(:,:), allocatable :: atomcharge_mpi      ! internal
      real*8, dimension(:), allocatable :: totalcharge_mpi       ! internal
      real*8, dimension(:,:,:), allocatable :: symfunctione_mpi  ! internal
      real*8, dimension(:), allocatable ::  nnewald_mpi          ! internal
      real*8, dimension(:,:), allocatable :: nnatomcharge_mpi    ! internal
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)         ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)         ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)          ! in
      real*8 nnewald_list(nblock)                     ! internal
      real*8 nnatomcharge_list(nblock,max_num_atoms)  ! internal
      real*8 avnncharge(nelem)                        ! internal
      real*8 bias                                         ! internal
!!
      logical, dimension(:), allocatable :: lperiodic_mpi    ! internal
!!
      if(mpisize.gt.1)then
        write(*,*)'ERROR: preconditioning does not yet work in parallel runs'
        stop !'
      endif
!!
!! write general header
      if(mpirank.eq.0)then
        write(ounit,*)'============================================================='
        write(ounit,*)'Weight Preconditioner:'     !'
        write(ounit,*)'----------------------'
      endif ! mpirank
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Preconditioner for the connecting weights to the output node and 
!! the bias of the output node 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the initial training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! initializations
      ncount                   = ntrain
      nnewald_list(:)          = 0.0d0
      nnatomcharge_list(:,:)   = 0.0d0 
      nnstddevq(:)             = 0.0d0
      ncharge(:)               = 0
      avnncharge(:)            = 0.0d0
!!
!! open files
      if(mpirank.eq.0)then
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
        open(symeunit,file='functione.data',form='formatted',status='old')
        rewind(symeunit) !'
      endif ! mpirank.eq.0
!!
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! determine which structures of this block should be calculated by this process
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! do all file reading for training error here at one place to allow for parallelization
!!--------------------------------------------------------------------------------------
!!
      if(mpirank.eq.0)then
!! read the symmetry functions for the charge prediction
        call readfunctions(1,symeunit,npoints,nelem,&
          max_num_atoms,maxnum_funcvalues_elec,num_funcvalues_elec,&
          symfunction_elec_list)
!! read the structures needed for the calculation of the electrostatic energy
        call getstructures(trainstructunit,npoints)
      endif ! mpirank.eq.0
!!
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalcharge_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(elecenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_elec_list,nblock*max_num_atoms*maxnum_funcvalues_elec,&
           mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomcharge_list,nblock*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
!!
!! end of file reading for training error
!!--------------------------------------------------------------------
!!
!! prepare local input arrays for this process
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
      allocate(elecenergy_mpi(nstruct))
      allocate(xyzstruct_mpi(3,max_num_atoms,nstruct))
      allocate(atomcharge_mpi(nstruct,max_num_atoms))
      allocate(lperiodic_mpi(nstruct))
      allocate(totalcharge_mpi(nstruct))
!! prepare local input/output arrays for this process
      allocate(symfunctione_mpi(maxnum_funcvalues_elec,max_num_atoms,nstruct))
!! prepare local output arrays for this process
      allocate(nnewald_mpi(nstruct))
      nnewald_mpi(:)=0.0d0
      allocate(nnatomcharge_mpi(nstruct,max_num_atoms))
      nnatomcharge_mpi(:,:)=0.0d0
!!
!!    copy local arrays
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        elecenergy_mpi(icount)      = elecenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
        atomcharge_mpi(icount,:)    = atomcharge_list(i1,:)
        totalcharge_mpi(icount)     = totalcharge_list(i1)
      enddo ! i1
!!
!! get the charges for this block of points 
!!
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        symfunctione_mpi(:,:,icount)=symfunction_elec_list(:,:,i1)
      enddo
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
        symfunctione_mpi,nnatomcharge_mpi)
!!
!! copy results back on full array
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        nnatomcharge_list(i1,:) = nnatomcharge_mpi(icount,:)
      enddo
!!
!! distribute results to all processes
      call mpi_allreduce(mpi_in_place,nnatomcharge_list,nblock*max_num_atoms,&
        mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!!    sum up charges 
      if(mpirank.eq.0)then
        do i1=1,npoints
          do i2=1,num_atoms_list(i1)
            ncharge(elementindex(zelem_list(i1,i2)))=&
              ncharge(elementindex(zelem_list(i1,i2)))+1
            avnncharge(elementindex(zelem_list(i1,i2)))&
              =avnncharge(elementindex(zelem_list(i1,i2)))+nnatomcharge_list(i1,i2)
          enddo ! i2
        enddo ! i1
      endif
!!
!!
      deallocate(num_atoms_mpi)
      deallocate(zelem_mpi)
      deallocate(elecenergy_mpi)
      deallocate(symfunctione_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(atomcharge_mpi)
      deallocate(lperiodic_mpi)
      deallocate(nnewald_mpi)
      deallocate(nnatomcharge_mpi)
      deallocate(totalcharge_mpi)
!!
      if(ncount.gt.0) goto 10
!! end block of training points
!!
!! close files
      if(mpirank.eq.0)then
        close(trainstructunit)
        close(symeunit)
      endif ! mpirank.eq.0
!!
      call mpi_allreduce(mpi_in_place,avnncharge,nelem,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,ncharge,nelem,&
           mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
!! now get averages
      do i1=1,nelem
        avnncharge(i1)=avnncharge(i1)/dble(ncharge(i1))
      enddo
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! Now get the standard deviations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!
      ncount                   = ntrain
!!
!! open files and write headers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mpirank.eq.0)then
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
        open(symeunit,file='functione.data',form='formatted',status='old')
        rewind(symeunit) !'
      endif ! mpirank.eq.0
!!
 20   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! determine which structures of this block should be calculated by this process
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! do all file reading for training error here at one place to allow for parallelization
!!--------------------------------------------------------------------------------------
!!
      if(mpirank.eq.0)then
!!
!! read the symmetry functions for the charge prediction
        call readfunctions(1,symeunit,npoints,nelem,&
          max_num_atoms,maxnum_funcvalues_elec,num_funcvalues_elec,&
          symfunction_elec_list)
!!
!! read the structures needed for the calculation of the electrostatic energy
        call getstructures(trainstructunit,npoints)
!!
      endif ! mpirank.eq.0
!!
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalcharge_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(elecenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_elec_list,nblock*max_num_atoms*maxnum_funcvalues_elec,&
           mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomcharge_list,nblock*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
!!
!! end of file reading for training error
!!--------------------------------------------------------------------
!!
!! prepare local input arrays for this process
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
      allocate(elecenergy_mpi(nstruct))
      allocate(xyzstruct_mpi(3,max_num_atoms,nstruct))
      allocate(atomcharge_mpi(nstruct,max_num_atoms))
      allocate(lperiodic_mpi(nstruct))
      allocate(totalcharge_mpi(nstruct))
!! prepare local input/output arrays for this process
      allocate(symfunctione_mpi(maxnum_funcvalues_elec,max_num_atoms,nstruct))
!! prepare local output arrays for this process
      allocate(nnewald_mpi(nstruct))
      nnewald_mpi(:)=0.0d0
      allocate(nnatomcharge_mpi(nstruct,max_num_atoms))
      nnatomcharge_mpi(:,:)=0.0d0
!!
!!    copy local arrays
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        elecenergy_mpi(icount)      = elecenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
        atomcharge_mpi(icount,:)    = atomcharge_list(i1,:)
        totalcharge_mpi(icount)     = totalcharge_list(i1)
      enddo ! i1
!!
!! get the charges for this block of points
!!
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        symfunctione_mpi(:,:,icount)=symfunction_elec_list(:,:,i1)
      enddo
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
        symfunctione_mpi,nnatomcharge_mpi)
!!
!! copy results back on full array
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        nnatomcharge_list(i1,:) = nnatomcharge_mpi(icount,:)
      enddo
!!
!! distribute results to all processes
      call mpi_allreduce(mpi_in_place,nnatomcharge_list,nblock*max_num_atoms,&
        mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! calculate part of the standard deviation   
      if(mpirank.eq.0)then
        do i1=1,npoints
          do i2=1,num_atoms_list(i1)
            nnstddevq(elementindex(zelem_list(i1,i2)))&
              =nnstddevq(elementindex(zelem_list(i1,i2)))&
              +(nnatomcharge_list(i1,i2)-avnncharge(elementindex(zelem_list(i1,i2))))**2.d0
          enddo
        enddo
      endif ! mpirank
!!
      deallocate(num_atoms_mpi)
      deallocate(zelem_mpi)
      deallocate(elecenergy_mpi)
      deallocate(symfunctione_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(atomcharge_mpi)
      deallocate(lperiodic_mpi)
      deallocate(nnewald_mpi)
      deallocate(nnatomcharge_mpi)
      deallocate(totalcharge_mpi)
!!
      if(ncount.gt.0) goto 20
!! end block of training points
!!
!! close files
      if(mpirank.eq.0)then
        close(trainstructunit)
        close(symeunit)
      endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a)')' Final preconditioning of the output values:' 
        write(ounit,*)    '--------------------------------------------'
!!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        do i1=1,nelem
          nnstddevq(i1)=nnstddevq(i1)/dble(ncharge(i1))
          nnstddevq(i1)=dsqrt(nnstddevq(i1))
        enddo
!! shift the charge bias weights
!! CHECK: According to Tobias this maybe should be done after scaling the connecting weights?
!          do i1=1,nelem
!            weights_ewald(num_weightsewald(i1),i1)&
!              =weights_ewald(num_weightsewald(i1),i1)-(avnncharge(i1)-avcharge(i1))   
!          enddo
!! scale the charge connecting weights between final hidden layer and output layer
        do i1=1,nelem
          prefactor=stddevcharge(i1)/nnstddevq(i1)
          write(ounit,'(a,a2,x,f14.6)')' Average NN charge ',element(i1),avnncharge(i1)
          write(ounit,'(a,a2,x,f14.6)')' Stddev NN charge  ',element(i1),nnstddevq(i1)
          write(ounit,'(a,a2,x,f14.6)')' Factor for connecting charge weights: ',element(i1),prefactor
          do i2=windex_elec(2*num_layers_elec(i1)-1,i1),windex_elec(2*num_layers_elec(i1),i1)
            weights_elec(i2,i1)=weights_elec(i2,i1)*prefactor
          enddo !'
        enddo
!! shift the charge bias weights, CHECK: should this be done here?
        do i1=1,nelem
          prefactor=stddevcharge(i1)/nnstddevq(i1)
          weights_elec(num_weights_elec(i1),i1)&
            =weights_elec(num_weights_elec(i1),i1)-(prefactor*avnncharge(i1)-avcharge(i1))   
        enddo
      endif ! mpirank
!!
      if(mpirank.eq.0)then
        write(ounit,*)'============================================================='
      endif ! mpirank
!!
      return
      end
