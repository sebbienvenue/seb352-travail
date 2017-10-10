!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine geterror_elec(iswitch,countepoch,ntrain,&
        imaxerror_elec,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_charge,rmse_totalcharge,&
        rmse_elec,&
        mad_charge,mad_totalcharge,&
        mad_elec,maxerror_elec)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use symfunctions
      use nnewald
      use structures
!!
      implicit none
!!
      integer i1,i2,i3      ! internal
      integer iswitch       ! in 0=train, 1=test
      integer ncharges      ! internal
      integer ncharges_sum  ! internal
      integer ntrain        ! in
      integer ndone         ! internal
      integer ndonepara     ! internal
      integer ncount        ! internal
      integer npoints       ! internal
      integer nstruct       ! internal
      integer n_start       ! internal
      integer n_end         ! internal
      integer imaxerror_elec                  ! out
      integer, dimension(:), allocatable :: num_atoms_mpi  ! internal
      integer, dimension(:,:), allocatable :: zelem_mpi    ! internal
      integer, dimension(:), allocatable :: imaxerrortemp  ! internal
      integer icount                         ! internal
      integer countepoch                     ! in
      integer nenergies                      ! internal
      integer tempunit                       ! internal
      integer pxunit                         ! internal
      integer qxunit                         ! internal
!!
      real*8, dimension(:), allocatable :: elecenergy_mpi        ! internal
      real*8, dimension(:,:,:), allocatable :: xyzstruct_mpi     ! internal
      real*8, dimension(:,:), allocatable :: atomcharge_mpi      ! internal
      real*8, dimension(:), allocatable :: totalcharge_mpi       ! internal
      real*8, dimension(:,:,:), allocatable :: symfunctione_mpi  ! internal
      real*8, dimension(:), allocatable ::  nnewald_mpi          ! internal
      real*8, dimension(:,:,:), allocatable :: nnelecforce_mpi   ! internal
      real*8, dimension(:,:), allocatable :: nnatomcharge_mpi    ! internal
      real*8, dimension(:), allocatable :: nnchargesum_mpi       ! internal
      real*8, dimension(:), allocatable :: maxerrortemp          ! internal
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)         ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)         ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)          ! in
      real*8 rmse_charge                                         ! out
      real*8 rmse_totalcharge                                    ! out
      real*8 rmse_elec                                          ! out
      real*8 mad_charge                                          ! out
      real*8 mad_totalcharge                                     ! out
      real*8 mad_elec                                           ! out
      real*8 maxerror_elec                                       ! out
      real*8 nnewald_list(nblock)                     ! internal
      real*8 nnatomcharge_list(nblock,max_num_atoms)  ! internal
      real*8 nnchargesum_list(nblock)                 ! internal
      real*8 maxforce_dummy                           ! internal
!!
      character*25 filename                           ! internal
!!
      logical, dimension(:), allocatable :: lperiodic_mpi    ! internal
!!
!!===============================================================
!!===============================================================
!! calculate the initial training error
!!===============================================================
!!===============================================================
!! initializations
!!===============================================================
      ncharges                 = 0
      ncharges_sum             = 0
      nenergies                = 0
      ncount                   = ntrain
      ndone                    = 0
      nnewald_list(:)          = 0.0d0
      nnatomcharge_list(:,:)   = 0.0d0 
      nnchargesum_list(:)      = 0.0d0
      maxforce_dummy           = 100000.d0
!!
!!===============================================================
!! open files and write headers
!!===============================================================
      if(mpirank.eq.0)then
        if(iswitch.eq.0)then ! train error
          open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
          rewind(trainstructunit)
          open(symeunit,file='functione.data',form='formatted',status='old')
          rewind(symeunit) !'
          if(lwritetraincharges)then
            filename='traincharges.000000.out'
            if(countepoch.gt.9999)then
              write(ounit,*)'Error: too many epochs in geterror'
              write(ounit,*)'switch off lwritetraincharges'
              stop
            elseif(countepoch.gt.999)then
              write(filename(16:19),'(i4)')countepoch
            elseif(countepoch.gt.99)then
              write(filename(17:19),'(i3)')countepoch
            elseif(countepoch.gt.9)then
              write(filename(18:19),'(i2)')countepoch
            else
              write(filename(19:19),'(i1)')countepoch
            endif
            open(trainqxunit,file=filename,form='formatted',status='replace')
            rewind(trainqxunit) !'
            write(trainqxunit,'(4a)')' point   atom    element      ',&
            'Q(DFT)     Q(NN)     Delta'
          endif
!!
        elseif(iswitch.eq.1)then ! test error
          open(teststructunit,file='teststruct.data',form='formatted',status='old')
          rewind(teststructunit)
          open(tymeunit,file='testinge.data',form='formatted',status='old')
          rewind(tymeunit)
          if(lwritetraincharges)then
            filename='testcharges.000000.out'
            if(countepoch.gt.999)then
              write(filename(15:18),'(i4)')countepoch
            elseif(countepoch.gt.99)then
              write(filename(16:18),'(i3)')countepoch
            elseif(countepoch.gt.9)then
              write(filename(17:18),'(i2)')countepoch
            else
              write(filename(18:18),'(i1)')countepoch
            endif
            open(testqxunit,file=filename,form='formatted',status='replace')
            rewind(testqxunit) !'
            write(testqxunit,'(4a)')' point   atom    element      ',&
            'Q(DFT)     Q(NN)     Delta'
          endif
        endif ! iswitch
      endif ! mpirank.eq.0
!!
!!===============================================================
!! loop over all structures
!!===============================================================
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!!===============================================================
!! determine which nstruct structures of this block should be calculated by this process
!!===============================================================
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!===============================================================
!! do all file reading for error determination here at one place to allow for parallelization
!!===============================================================
      if(mpirank.eq.0)then
!!===============================================================
!! read the symmetry functions for the charge prediction (train or test)
!!===============================================================
        if(iswitch.eq.0)then ! train
          tempunit=symeunit          
        elseif(iswitch.eq.1)then ! test
          tempunit=tymeunit          
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
        call readfunctions(1,tempunit,npoints,nelem,&
          max_num_atoms,maxnum_funcvalues_elec,num_funcvalues_elec,&
          symfunction_elec_list)
!!
!!===============================================================
!! read the structures needed for the calculation of the electrostatic energy
!! is needed for lelec (structure for electrostatics)
!! must be called after readfunctions because it needs num_atoms_list
!!===============================================================
        if(iswitch.eq.0)then ! train
          tempunit=trainstructunit          
        elseif(iswitch.eq.1)then ! test
          tempunit=teststructunit          
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
        call getstructures(tempunit,npoints)
!!
      endif ! mpirank.eq.0
!!
!!===============================================================
!! distribute the data to all processes
!!===============================================================
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
!!===============================================================
!! end of file reading for error determination
!!===============================================================
!!
!!===============================================================
!! prepare local arrays for this process
!!===============================================================
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
      allocate(elecenergy_mpi(nstruct))
      allocate(xyzstruct_mpi(3,max_num_atoms,nstruct))
      allocate(atomcharge_mpi(nstruct,max_num_atoms))
      allocate(lperiodic_mpi(nstruct))
      allocate(totalcharge_mpi(nstruct))
      allocate(symfunctione_mpi(maxnum_funcvalues_elec,max_num_atoms,nstruct))
      allocate(nnewald_mpi(nstruct))
      nnewald_mpi(:)=0.0d0
      allocate(nnelecforce_mpi(3,max_num_atoms,nstruct))
      nnelecforce_mpi(:,:,:)=0.0d0
      allocate(nnatomcharge_mpi(nstruct,max_num_atoms))
      nnatomcharge_mpi(:,:)=0.0d0
      allocate(nnchargesum_mpi(nstruct))
      nnchargesum_mpi(:)=0.0d0
!!
!!===============================================================
!! fill local input arrays
!!===============================================================
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
!!===============================================================
!! get the electrostatic energies 
!!===============================================================
      ndonepara=ndone+n_start-1
      call ewaldenergies_para(nstruct,ndonepara,&
        imaxerror_elec,num_atoms_mpi,zelem_mpi,ncharges,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        symfunction_elec_list(1,1,n_start),lattice_list(1,1,n_start),&
        nnelecforce_mpi,xyzstruct_mpi,nnewald_mpi,&
        rmse_charge,rmse_totalcharge,rmse_elec,&
        mad_charge,mad_totalcharge,mad_elec,maxerror_elec,&
        elecenergy_mpi,atomcharge_mpi,nnatomcharge_mpi,&
        totalcharge_mpi,nnchargesum_mpi,&
        lperiodic_mpi)
!!
!!===============================================================
!! for parallel case: find imaxerror_elec and maxerror_elec from all processes
!!===============================================================
      if(mpisize.gt.1)then
        allocate(imaxerrortemp(mpisize))
        imaxerrortemp(:)=0
        allocate(maxerrortemp(mpisize))
        maxerrortemp(:) =0.0d0
        imaxerrortemp(mpirank+1)=imaxerror_elec
        maxerrortemp(mpirank+1) =maxerror_elec
        call mpi_allreduce(mpi_in_place,imaxerrortemp,mpisize,&
          mpi_integer,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,maxerrortemp,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        do i1=1,mpisize
          if(maxerrortemp(i1).gt.maxerror_elec)then
            maxerror_elec = maxerrortemp(i1)
            imaxerror_elec= imaxerrortemp(i1)
          endif
        enddo ! i1
        deallocate(imaxerrortemp)
        deallocate(maxerrortemp)
      endif ! mpisize
!!
!!===============================================================
!! copy results back on full array
!!===============================================================
      icount=0
      nnewald_list(:)           =0.0d0
      nnatomcharge_list(:,:)    =0.0d0
      nnchargesum_list(:)       =0.0d0
      do i1=n_start,n_end
        icount=icount+1
        nnewald_list(i1)          = nnewald_mpi(icount)
        nnatomcharge_list(i1,:)   = nnatomcharge_mpi(icount,:)
        nnchargesum_list(i1)      = nnchargesum_mpi(icount)
      enddo
!!
!!===============================================================
!! distribute results to all processes
!!===============================================================
      call mpi_allreduce(mpi_in_place,nnewald_list,nblock,&
        mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,nnatomcharge_list,nblock*max_num_atoms,&
        mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,nnchargesum_list,nblock,&
        mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!!===============================================================
!! prepare file units for detailed output 
!!===============================================================
      if(iswitch.eq.0)then ! train set
        pxunit=trainpxunit
        qxunit=trainqxunit
      elseif(iswitch.eq.1)then ! test set
        pxunit=testpxunit
        qxunit=testqxunit
      else
        write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
        stop
      endif
!!
      if(mpirank.eq.0)then
!!===============================================================
!! if requested write energy details for all points to trainpoints.XXXXXX.out
!!===============================================================
        if(lwritetrainpoints)then
          do i1=1,npoints
            write(pxunit,'(i6,x,4f14.8)')i1+ndone,&
            elecenergy_list(i1),nnewald_list(i1),&
            totalcharge_list(i1),nnchargesum_list(i1)
          enddo
        endif
!!
!!===============================================================
!! if requested write charge details for all points to traincharges.XXXXXX.out
!!===============================================================
        if(lwritetraincharges)then
          if(lelec)then
            do i1=1,npoints
              do i2=1,num_atoms_list(i1)
                write(qxunit,'(i6,x,i6,x,i3,x,3f14.8)')i1+ndone,i2,&
                  zelem_list(i1,i2),&
                  atomcharge_list(i1,i2),nnatomcharge_list(i1,i2),&
                  atomcharge_list(i1,i2)-nnatomcharge_list(i1,i2)
              enddo ! i2
            enddo ! i1
          endif ! lelec
        endif ! lwritetraincharges
!!
        ndone=ndone+npoints
!!
      endif ! mpirank.eq.0
!!
!!===============================================================
!! deallocate temporary arrays for this process
!!===============================================================
      deallocate(num_atoms_mpi)
      deallocate(zelem_mpi)
      deallocate(elecenergy_mpi)
      deallocate(symfunctione_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(atomcharge_mpi)
      deallocate(lperiodic_mpi)
      deallocate(nnewald_mpi)
      deallocate(nnelecforce_mpi)
      deallocate(nnatomcharge_mpi)
      deallocate(totalcharge_mpi)
      deallocate(nnchargesum_mpi)
!!
!!===============================================================
!! if there are structures left go to next block of structures 
!!===============================================================
      if(ncount.gt.0) goto 10
!!===============================================================
!! end block of structures 
!!===============================================================
!!
!!===============================================================
!! close files
!!===============================================================
      if(mpirank.eq.0)then
        if(iswitch.eq.0)then ! training
          close(trainstructunit)
          close(symeunit)
          if(lwritetraincharges)     close(trainqxunit)
          if(lwritetrainpoints)      close(trainpxunit)
        elseif(iswitch.eq.1)then ! testing
          close(teststructunit)
          close(tymeunit)
          if(lwritetraincharges)     close(testqxunit)
          if(lwritetrainpoints)      close(testpxunit)
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif ! iswitch
      endif ! mpirank.eq.0
!!
!!===============================================================
!! combine the partial rmse values of all processes here only to avoid double counting
!!===============================================================
      call mpi_allreduce(mpi_in_place,rmse_charge,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,rmse_totalcharge,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,rmse_elec,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,mad_charge,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,mad_totalcharge,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,mad_elec,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,ncharges,1,&
           mpi_integer,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,nenergies,1,&
           mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
      ncharges_sum=ncharges_sum+ncharges
      ncharges=0
!!
!!===============================================================
!! calculate the final RMSEs 
!!===============================================================
      call getrmse_elec(ntrain,ncharges_sum,&
        rmse_charge,rmse_totalcharge,rmse_elec,&
        mad_charge,mad_totalcharge,mad_elec)
!!
      return
      end
