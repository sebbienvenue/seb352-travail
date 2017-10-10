!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine geterror_atomic(iswitch,countepoch,ntrain,&
        imaxerror_eshort,imaxerror_elec,imaxerror_etot,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        rmse_short,rmse_charge,rmse_totalcharge,&
        rmse_elec,rmse_etot,&
        rmse_force_s,rmse_force_t,rmse_force_e,&
        mad_short,mad_charge,mad_totalcharge,&
        mad_elec,mad_etot,&
        mad_force_s,mad_force_t,mad_force_e,&
        maxerror_eshort,maxerror_elec,maxerror_etot)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use symfunctions
      use nnshort_atomic
      use structures
!!
      implicit none
!!
      integer i1,i2,i3                                     ! internal
      integer iswitch                                      ! in 0=train, 1=test
      integer ncharges                                     ! internal
      integer ncharges_sum                                 ! internal
      integer nforces_short                                ! internal
      integer nforces_elec                                 ! internal
      integer nforces_total                                ! internal
      integer ntrain                                       ! in
      integer ndone                                        ! internal
      integer ndonepara                                    ! internal
      integer ncount                                       ! internal
      integer npoints                                      ! internal
      integer nstruct                                      ! internal
      integer n_start                                      ! internal
      integer n_end                                        ! internal
      integer imaxerror_eshort                             ! out
      integer imaxerror_elec                               ! out
      integer imaxerror_etot                               ! out
      integer, dimension(:), allocatable :: num_atoms_mpi  ! internal
      integer, dimension(:,:), allocatable :: zelem_mpi    ! internal
      integer, dimension(:), allocatable :: imaxerror_temp ! internal
      integer icount                                       ! internal
      integer countepoch                                   ! in
      integer nenergies                                    ! internal
      integer netot                                        ! internal
      integer tempunit                                     ! internal
      integer pxunit                                       ! internal
      integer fxunit                                       ! internal
      integer qxunit                                       ! internal
!!
      real*8, dimension(:), allocatable :: elecenergy_mpi        ! internal
      real*8, dimension(:), allocatable :: shortenergy_mpi       ! internal
      real*8, dimension(:,:,:), allocatable :: xyzstruct_mpi     ! internal
      real*8, dimension(:,:), allocatable :: atomcharge_mpi      ! internal
      real*8, dimension(:), allocatable :: totalcharge_mpi       ! internal
      real*8, dimension(:), allocatable ::  nneshort_mpi         ! internal
      real*8, dimension(:,:,:), allocatable :: nnshortforce_mpi  ! internal
      real*8, dimension(:), allocatable ::  nnelec_mpi           ! internal
      real*8, dimension(:,:,:), allocatable :: nnelecforce_mpi   ! internal
      real*8, dimension(:,:), allocatable :: nnatomcharge_mpi    ! internal
      real*8, dimension(:), allocatable :: nnchargesum_mpi       ! internal
      real*8, dimension(:), allocatable :: maxerror_temp         ! internal
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)! in
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic) ! in
      real*8 rmse_short                               ! out
      real*8 rmse_charge                              ! out
      real*8 rmse_totalcharge                         ! out
      real*8 rmse_elec                                ! out
      real*8 rmse_etot                                ! out
      real*8 rmse_force_s                             ! out
      real*8 rmse_force_t                             ! out
      real*8 rmse_force_e                             ! out
      real*8 mad_short                                ! out
      real*8 mad_charge                               ! out
      real*8 mad_totalcharge                          ! out
      real*8 mad_elec                                 ! out
      real*8 mad_etot                                 ! out
      real*8 mad_force_s                              ! out
      real*8 mad_force_t                              ! out
      real*8 mad_force_e                              ! out
      real*8 maxerror_eshort                          ! out
      real*8 maxerror_elec                            ! out
      real*8 maxerror_etot                            ! out
      real*8 nneshort_list(nblock)                    ! internal
      real*8 nnshortforce_list(3,max_num_atoms,nblock)! internal
      real*8 nntotforce_list(3,max_num_atoms,nblock)  ! internal
      real*8 nnelec_list(nblock)                      ! internal
      real*8 nnelecforce_list(3,max_num_atoms,nblock) ! internal
      real*8 nnatomcharge_list(nblock,max_num_atoms)  ! internal
      real*8 nnchargesum_list(nblock)                 ! internal
      real*8 nnetot_list(nblock)                      ! internal
      real*8 edummy                                   ! internal
      real*8 forcesum(3)                              ! internal
      real*8 maxforce_dummy                           ! internal
!!
      character*25 filename                           ! internal
!!
      logical, dimension(:), allocatable :: lperiodic_mpi    ! internal
!!
!!==============================================
!! initializations
!!==============================================
!! counters
      nenergies                = 0
      netot                    = 0
      nforces_short            = 0
      nforces_elec             = 0
      nforces_total            = 0
      ncharges                 = 0
      ncharges_sum             = 0
      ncount                   = ntrain
      ndone                    = 0
!! reference data
      shortforce_list(:,:,:)   = 0.0d0
      elecforce_list(:,:,:)    = 0.0d0
!! data
      nnelec_list(:)           = 0.0d0
      nneshort_list(:)         = 0.0d0
      nnelecforce_list(:,:,:)  = 0.0d0
      nnshortforce_list(:,:,:) = 0.0d0
      nntotforce_list(:,:,:)   = 0.0d0
      nnatomcharge_list(:,:)   = 0.0d0 
      nnchargesum_list(:)      = 0.0d0
      nnetot_list(:)           = 0.0d0
!! others
      edummy                   = 1.d12
      maxforce_dummy           = 100000.d0
!!
!!==============================================
!! open files and write headers
!!==============================================
      if(mpirank.eq.0)then
        if(iswitch.eq.0)then
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit) !'
          if(lwritetrainpoints)then
            filename='trainpoints.000000.out'
            if(countepoch.gt.9999)then
              write(ounit,*)'Error: too many epochs in geterror'
              write(ounit,*)'switch off lwritetrainpoints'
              stop
            elseif(countepoch.gt.999)then
              write(filename(15:18),'(i4)')countepoch
            elseif(countepoch.gt.99)then
              write(filename(16:18),'(i3)')countepoch
            elseif(countepoch.gt.9)then
              write(filename(17:18),'(i2)')countepoch
            else
              write(filename(18:18),'(i1)')countepoch
            endif
            open(trainpxunit,file=filename,form='formatted',status='replace')
            rewind(trainpxunit) !'
            write(trainpxunit,'(4a)')' point      Etot(DFT)      Etot(NN)  ',&
             'E_short(DFT)   E_short(NN)  ',&
             'E_elec(DFT)    E_elec(NN)      ',&
             'Qsum(DFT)     Qsum(NN)'
          endif
          open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
          rewind(trainstructunit) !'
!!
          if(lelec)then
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
          endif ! lelec
!!
          if(luseforces)then
            if(lwritetrainforces)then
              filename='trainforces.000000.out'
              if(countepoch.gt.9999)then
                write(ounit,*)'Error: too many epochs in geterror'
                write(ounit,*)'switch off lwritetrainforces'
                stop
              elseif(countepoch.gt.999)then
                write(filename(15:18),'(i4)')countepoch
              elseif(countepoch.gt.99)then
                write(filename(16:18),'(i3)')countepoch
              elseif(countepoch.gt.9)then
                write(filename(17:18),'(i2)')countepoch
              else
                write(filename(18:18),'(i1)')countepoch
              endif
              open(trainfxunit,file=filename,form='formatted',status='replace')
              write(trainfxunit,'(4a)')' point   atom       F_s(DFT)       F_s(NN)',&
              '         Delta          ',&
              'F_e(DFT)       F_e(NN)         Delta'
            endif
            open(trainfunit,file='trainforces.data',form='formatted',status='old')
            rewind(trainfunit)
            if(lelec)then
              open(trainfeunit,file='trainforcese.data',form='formatted',status='old')
              rewind(trainfeunit)
            endif
          endif !' luseforces
!!
        elseif(iswitch.eq.1)then ! test error
          open(tymunit,file='testing.data',form='formatted',status='old')
          rewind(tymunit)
          if(lwritetrainpoints)then
            filename='testpoints.000000.out'
            if(countepoch.gt.999)then
              write(filename(16:17),'(i4)')countepoch
            elseif(countepoch.gt.99)then
              write(filename(15:17),'(i3)')countepoch
            elseif(countepoch.gt.9)then
              write(filename(16:17),'(i2)')countepoch
            else
              write(filename(17:17),'(i1)')countepoch
            endif
            open(testpxunit,file=filename,form='formatted',status='replace')
            rewind(testpxunit) !'
            write(testpxunit,'(4a)')' point      Etot(DFT)      Etot(NN)  ',&
             'E_short(DFT)   E_short(NN)  ',&
             'E_elec(DFT)    E_elec(NN)      ',&
             'Qsum(DFT)     Qsum(NN)'
          endif
          open(teststructunit,file='teststruct.data',form='formatted',status='old')
          rewind(teststructunit) !'
!!
          if(lelec)then
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
          endif ! lelec
!!
          if(luseforces)then
            if(lwritetrainforces)then
              filename='testforces.000000.out'
              if(countepoch.gt.999)then
                write(filename(14:17),'(i4)')countepoch
              elseif(countepoch.gt.99)then
                write(filename(15:17),'(i3)')countepoch
              elseif(countepoch.gt.9)then
                write(filename(16:17),'(i2)')countepoch
              else
                write(filename(17:17),'(i1)')countepoch
              endif
              open(testfxunit,file=filename,form='formatted',status='replace')
              write(testfxunit,'(4a)')' point   atom       F_s(DFT)       F_s(NN)',&
              '         Delta          ',&
              'F_e(DFT)       F_e(NN)         Delta'
            endif
            open(testfunit,file='testforces.data',form='formatted',status='old')
            rewind(testfunit)
            if(lelec)then
              open(testfeunit,file='testforcese.data',form='formatted',status='old')
              rewind(testfeunit)
            endif
          endif !' luseforces
        endif ! iswitch
      endif ! mpirank.eq.0
!!
!!==============================================
!! loop block-wise over all structures 
!!==============================================
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!!==============================================
!! determine which nstruct structures of this block should be calculated by this process
!!==============================================
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!==============================================
!!==============================================
!! do all file reading here for nstruct structures at one place to allow for parallelization
!!==============================================
!!==============================================
      if(mpirank.eq.0)then
!!==============================================
!! read npoint short range symmetry function sets (train or test)
!!==============================================
        if(iswitch.eq.0)then ! train
          tempunit=symunit          
        elseif(iswitch.eq.1)then ! test
          tempunit=tymunit          
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
        call readfunctions(1,tempunit,npoints,nelem,&
          max_num_atoms,maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
          symfunction_short_atomic_list)
!!
!!==============================================
!! read the structures needed for the calculation of the electrostatic energy
!! and get reference forces from DFT
!! must be called after readfunctions because it needs num_atoms_list
!!==============================================
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
!!==============================================
!! read short range forces from trainforces.data or testforces.data => shortforce_list
!!==============================================
        if(luseforces)then
          if(iswitch.eq.0)then ! train
            tempunit=trainfunit          
          elseif(iswitch.eq.1)then ! test
            tempunit=testfunit          
          else
            write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
            stop
          endif
          call readforces(tempunit,npoints,shortforce_list)
        endif
!!
!!==============================================
!! read electrostatic forces from trainforcese.data or testforcese.data => elecforce_list
!!==============================================
        if(luseforces.and.lelec)then
          if(iswitch.eq.0)then ! train
            tempunit=trainfeunit          
          elseif(iswitch.eq.1)then ! test
            tempunit=testfeunit          
          else
            write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
            stop
          endif
          call readforces(tempunit,npoints,elecforce_list)
        endif
      endif ! mpirank.eq.0
!!
!!==============================================
!! distribute the data to all processes
!!==============================================
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_short_atomic_list,nblock*max_num_atoms*maxnum_funcvalues_short_atomic,&
           mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
      if(lelec.and.(nn_type_elec.eq.2))then
        call mpi_bcast(elecenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(atomcharge_list,nblock*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(totalcharge_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(luseforces)then
        call mpi_bcast(shortforce_list,nblock*max_num_atoms*3,&
          mpi_real8,0,mpi_comm_world,mpierror)
        if(lelec.and.(nn_type_elec.eq.2))then
          call mpi_bcast(elecforce_list,nblock*max_num_atoms*3,&
            mpi_real8,0,mpi_comm_world,mpierror)
        endif
      endif
!!==============================================
!!==============================================
!! end of file reading 
!!==============================================
!!==============================================
!!
!!==============================================
!! allocate local arrays for parallel processes 
!!==============================================
!! prepare local input arrays for this process
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
      allocate(elecenergy_mpi(nstruct))
      allocate(shortenergy_mpi(nstruct))
      allocate(xyzstruct_mpi(3,max_num_atoms,nstruct))
      allocate(atomcharge_mpi(nstruct,max_num_atoms))
      allocate(lperiodic_mpi(nstruct))
      allocate(totalcharge_mpi(nstruct))
!! prepare local output arrays for this process
      allocate(nneshort_mpi(nstruct))
      nneshort_mpi(:)=0.0d0
      allocate(nnshortforce_mpi(3,max_num_atoms,nstruct))
      nnshortforce_mpi(:,:,:)=0.0d0
      allocate(nnelec_mpi(nstruct))
      nnelec_mpi(:)=0.0d0
      allocate(nnelecforce_mpi(3,max_num_atoms,nstruct))
      nnelecforce_mpi(:,:,:)=0.0d0
      allocate(nnatomcharge_mpi(nstruct,max_num_atoms))
      nnatomcharge_mpi(:,:)=0.0d0
      allocate(nnchargesum_mpi(nstruct))
      nnchargesum_mpi(:)=0.0d0
!!
!!==============================================
!! fill local arrays for structures n_start to n_end
!!==============================================
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        shortenergy_mpi(icount)     = shortenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
        if(lelec)then
          elecenergy_mpi(icount)    = elecenergy_list(i1)
          atomcharge_mpi(icount,:)  = atomcharge_list(i1,:)
          totalcharge_mpi(icount)   = totalcharge_list(i1)
        endif ! lelec
      enddo ! i1
!!
!!==============================================
!! get the energies, forces and charges 
!!==============================================
        ndonepara=ndone+n_start-1
!! FIXME: getatomicoutput_para is not completed yet
        call getatomicoutput_para(nstruct,ndonepara,&
          nenergies,ncharges,&
          imaxerror_eshort,imaxerror_elec,imaxerror_etot,&
          num_atoms_mpi,zelem_mpi,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
          symfunction_short_atomic_list(1,1,n_start),&
          lattice_list(1,1,n_start),xyzstruct_mpi,&
          rmse_short,rmse_elec,rmse_charge,rmse_totalcharge,&
          mad_short,mad_elec,mad_charge,mad_totalcharge,&
          maxerror_eshort,maxerror_elec,maxerror_etot,&
          shortenergy_mpi,elecenergy_mpi,atomcharge_mpi,totalcharge_mpi,&
          nneshort_mpi,nnshortforce_mpi,nnatomcharge_mpi,nnchargesum_mpi,&
          nnelec_mpi,nnelecforce_mpi,&
          lperiodic_mpi)
!!
!!==============================================
!! for parallel case: find imaxerror_eshort and maxerror_eshort from all processes
!!==============================================
        if(mpisize.gt.1)then
          allocate(imaxerror_temp(mpisize))
          imaxerror_temp(:)=0
          allocate(maxerror_temp(mpisize))
          maxerror_temp(:) =0.0d0
          imaxerror_temp(mpirank+1)=imaxerror_eshort
          maxerror_temp(mpirank+1) =maxerror_eshort
          call mpi_allreduce(mpi_in_place,imaxerror_temp,mpisize,&
            mpi_integer,mpi_sum,mpi_comm_world,mpierror)
          call mpi_allreduce(mpi_in_place,maxerror_temp,mpisize,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
          do i1=1,mpisize
            if(maxerror_temp(i1).gt.maxerror_eshort)then
              maxerror_eshort = maxerror_temp(i1)
              imaxerror_eshort= imaxerror_temp(i1)
            endif
          enddo ! i1
          deallocate(imaxerror_temp)
          deallocate(maxerror_temp)
        endif
!!==============================================
!! for parallel case: find imaxerror_elec and maxerror_elec from all processes
!!==============================================
        if(lelec)then
          if(mpisize.gt.1)then
            allocate(imaxerror_temp(mpisize))
            imaxerror_temp(:)=0
            allocate(maxerror_temp(mpisize))
            maxerror_temp(:) =0.0d0
            imaxerror_temp(mpirank+1)=imaxerror_elec
            maxerror_temp(mpirank+1) =maxerror_elec
            call mpi_allreduce(mpi_in_place,imaxerror_temp,mpisize,&
              mpi_integer,mpi_sum,mpi_comm_world,mpierror)
            call mpi_allreduce(mpi_in_place,maxerror_temp,mpisize,&
              mpi_real8,mpi_sum,mpi_comm_world,mpierror)
            do i1=1,mpisize
              if(maxerror_temp(i1).gt.maxerror_elec)then
                maxerror_elec = maxerror_temp(i1)
                imaxerror_elec= imaxerror_temp(i1)
              endif
            enddo ! i1
            deallocate(imaxerror_temp)
            deallocate(maxerror_temp)
          endif
        endif ! lelec
!!
!!==============================================
!! for parallel case: find imaxerror_etot and maxerror_etot from all processes
!!==============================================
        if(lelec)then
          if(mpisize.gt.1)then
            allocate(imaxerror_temp(mpisize))
            imaxerror_temp(:)=0
            allocate(maxerror_temp(mpisize))
            maxerror_temp(:) =0.0d0
            imaxerror_temp(mpirank+1)=imaxerror_etot
            maxerror_temp(mpirank+1) =maxerror_etot
            call mpi_allreduce(mpi_in_place,imaxerror_temp,mpisize,&
              mpi_integer,mpi_sum,mpi_comm_world,mpierror)
            call mpi_allreduce(mpi_in_place,maxerror_temp,mpisize,&
              mpi_real8,mpi_sum,mpi_comm_world,mpierror)
            do i1=1,mpisize
              if(maxerror_temp(i1).gt.maxerror_etot)then
                maxerror_etot = maxerror_temp(i1)
                imaxerror_etot= imaxerror_temp(i1)
              endif
            enddo ! i1
            deallocate(imaxerror_temp)
            deallocate(maxerror_temp)
          endif
        endif ! lelec
!!
!!==============================================
!! copy results back on full array
!!==============================================
        icount=0
        nneshort_list(:)         =0.0d0
        nnshortforce_list(:,:,:) =0.0d0
        if(lelec)then
          nnelec_list(:)            =0.0d0
          nnelecforce_list(:,:,:)   =0.0d0
          nnatomcharge_list(:,:)    =0.0d0
          nnchargesum_list(:)       =0.0d0
        endif ! lelec
        do i1=n_start,n_end
          icount=icount+1
          nneshort_list(i1)         = nneshort_mpi(icount)
          nnshortforce_list(:,:,i1) = nnshortforce_mpi(:,:,icount)
          if(lelec)then
            nnelec_list(i1)           = nnelec_mpi(icount)
            nnelecforce_list(:,:,i1)  = nnelecforce_mpi(:,:,icount)
            nnatomcharge_list(i1,:)   = nnatomcharge_mpi(icount,:)
            nnchargesum_list(i1)      = nnchargesum_mpi(icount)
          endif
        enddo
!!
!!==============================================
!! distribute results to all processes
!!==============================================
        call mpi_allreduce(mpi_in_place,nneshort_list,nblock,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,nnshortforce_list,nblock*max_num_atoms*3,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(lelec)then
          call mpi_allreduce(mpi_in_place,nnelec_list,nblock,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
          call mpi_allreduce(mpi_in_place,nnelecforce_list,nblock*max_num_atoms*3,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
          call mpi_allreduce(mpi_in_place,nnatomcharge_list,nblock*max_num_atoms,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
          call mpi_allreduce(mpi_in_place,nnchargesum_list,nblock,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif ! lelec
!!
!!==============================================
!! check NN forces if requested
!!==============================================
        if(mpirank.eq.0)then
          if(lcheckf)then
            do i1=1,npoints
              forcesum(:)=0.0d0
              do i2=1,num_atoms_list(i1)
                do i3=1,3
                  forcesum(i3)=forcesum(i3)+nnshortforce_list(i3,i2,i1)
                enddo ! i3 
              enddo ! i2
              do i2=1,3
                if(abs(forcesum(i2)).gt.0.000001d0)then
                  write(ounit,*)'Error in force sum ',i1,i2,forcesum(i2)
                  stop
                endif
              enddo ! i2
            enddo ! i1
          endif ! lcheckf
        endif ! mpirank.eq.0
!!
!!==============================================
!! get the total energy from short range and electrostatic contributions
!!==============================================
      nnetot_list(:)=nneshort_list(:)+nnelec_list(:)
!!
!!==============================================
!! calculate the RMSE for the total energy: in/out is rmse_etot
!!==============================================
      ndonepara=ndone+n_start-1
      call calcrmse_energy(nblock,npoints,netot,&
        ndonepara,imaxerror_etot,&
        rmse_etot,mad_etot,maxerror_etot,&
        edummy,totalenergy_list,nnetot_list)
!!
!!==============================================
!! calculate the RMSE of the total forces here: in/out rmse_force_s + nforces_short
!!==============================================
      if(luseforces.and.lshort)then
        call calcrmse_forces(npoints,&
          nforces_short,rmse_force_s,mad_force_s,&
          shortforce_list,nnshortforce_list,maxforce)
      endif
!!
!!==============================================
!! calculate the RMSE of the electrostatic forces here: in/out rmse_force_e + nforces_elec
!!==============================================
      if(luseforces.and.lelec)then
        call calcrmse_forces(npoints,&
          nforces_elec,rmse_force_e,mad_force_e,&
          elecforce_list,nnelecforce_list,maxforce_dummy)
      endif
!!
!!==============================================
!! calculate the RMSE of the total forces here: in/out rmse_force_t + nforces_total
!!==============================================
      if(luseforces)then
        if((lshort).and.(.not.lelec))then
          totforce_list(:,:,:)  =shortforce_list(:,:,:)
          nntotforce_list(:,:,:)=nnshortforce_list(:,:,:)
          call calcrmse_forces(npoints,&
            nforces_total,rmse_force_t,mad_force_t,&
            totforce_list,nntotforce_list,maxforce_dummy)
        elseif((.not.lshort).and.(lelec))then
          totforce_list(:,:,:)  =  elecforce_list(:,:,:)
          nntotforce_list(:,:,:)=nnelecforce_list(:,:,:)
          call calcrmse_forces(npoints,&
            nforces_total,rmse_force_t,mad_force_t,&
            totforce_list,nntotforce_list,maxforce_dummy)
        elseif(lshort.and.lelec)then
          totforce_list(:,:,:)  =  shortforce_list(:,:,:)+  elecforce_list(:,:,:)
          nntotforce_list(:,:,:)=nnshortforce_list(:,:,:)+nnelecforce_list(:,:,:)
          call calcrmse_forces(npoints,&
            nforces_total,rmse_force_t,mad_force_t,&
            totforce_list,nntotforce_list,maxforce_dummy)
        else ! must not be possible:
          write(ounit,*)'Error in geterror_atomic'
          stop
        endif
      endif
!!
!!==============================================
!! prepare file units for output
!!==============================================
      if(mpirank.eq.0)then
        if(iswitch.eq.0)then ! train set
          pxunit=trainpxunit
          qxunit=trainqxunit
          fxunit=trainfxunit
        elseif(iswitch.eq.1)then ! test set
          pxunit=testpxunit
          qxunit=testqxunit
          fxunit=testfxunit
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
!!
!!==============================================
!! if requested write energy details for all points to trainpoints.XXXXXX.out
!!==============================================
        if(lwritetrainpoints)then
          do i1=1,npoints
            write(pxunit,'(i6,x,8f14.8)')i1+ndone,totalenergy_list(i1),nnetot_list(i1),&
            shortenergy_list(i1),nneshort_list(i1),elecenergy_list(i1),nnelec_list(i1),&
            totalcharge_list(i1),nnchargesum_list(i1)
          enddo
        endif
!!
!!==============================================
!! if requested write charge details for all points to traincharges.XXXXXX.out
!!==============================================
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
!!==============================================
!! write force details for all points to trainforces.XXXXXX.out
!!==============================================
        if(luseforces.and.lwritetrainforces)then
          do i1=1,npoints
            do i2=1,num_atoms_list(i1)
            write(fxunit,'(2i5,a4,2f14.6,2x,f14.8,2x,2f14.6,2x,f14.8)')&
              i1+ndone,i2,' Fx ',&
              shortforce_list(1,i2,i1),nnshortforce_list(1,i2,i1),&
              shortforce_list(1,i2,i1)-nnshortforce_list(1,i2,i1),&
              elecforce_list(1,i2,i1),nnelecforce_list(1,i2,i1),&
              elecforce_list(1,i2,i1)-nnelecforce_list(1,i2,i1)
            write(fxunit,'(2i5,a4,2f14.6,2x,f14.8,2x,2f14.6,2x,f14.8)')&
              i1+ndone,i2,' Fy ',&
              shortforce_list(2,i2,i1),nnshortforce_list(2,i2,i1),&
              shortforce_list(2,i2,i1)-nnshortforce_list(2,i2,i1),&
              elecforce_list(2,i2,i1),nnelecforce_list(2,i2,i1),&
              elecforce_list(2,i2,i1)-nnelecforce_list(2,i2,i1)
            write(fxunit,'(2i5,a4,2f14.6,2x,f14.8,2x,2f14.6,2x,f14.8)')&
              i1+ndone,i2,' Fz ',&
              shortforce_list(3,i2,i1),nnshortforce_list(3,i2,i1),&
              shortforce_list(3,i2,i1)-nnshortforce_list(3,i2,i1),&
              elecforce_list(3,i2,i1),nnelecforce_list(3,i2,i1),&
              elecforce_list(3,i2,i1)-nnelecforce_list(3,i2,i1)
            enddo ! i2
          enddo ! i1
        endif ! luseforces
!!
        ndone=ndone+npoints
!!
      endif ! mpirank.eq.0
!!
!!==============================================
!! cleanup
!!==============================================
      deallocate(num_atoms_mpi)
      deallocate(zelem_mpi)
      deallocate(elecenergy_mpi)
      deallocate(shortenergy_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(atomcharge_mpi)
      deallocate(lperiodic_mpi)
      deallocate(nneshort_mpi)
      deallocate(nnshortforce_mpi)
      deallocate(nnelec_mpi)
      deallocate(nnelecforce_mpi)
      deallocate(nnatomcharge_mpi)
      deallocate(totalcharge_mpi)
      deallocate(nnchargesum_mpi)
!!
!!==============================================
!! if there are structures left to calculate go to next block of structures
!!==============================================
      if(ncount.gt.0) goto 10
!!==============================================
!! looping over structures done
!!==============================================
!!
!!==============================================
!! close files
!!==============================================
      if(mpirank.eq.0)then
        if(iswitch.eq.0)then
          close(symunit)
          close(trainstructunit)
          if(luseforces)                       close(trainfunit)
          if(lelec.and.luseforces)             close(trainfeunit)
          if(lelec.and.lwritetraincharges)     close(trainqxunit)
          if(lwritetrainpoints)                close(trainpxunit)
          if(luseforces.and.lwritetrainforces) close(trainfxunit)
        elseif(iswitch.eq.1)then
          close(tymunit)
          close(teststructunit)
          if(luseforces)                       close(testfunit)
          if(lelec.and.luseforces)             close(testfeunit)
          if(lelec.and.lwritetraincharges)     close(testqxunit)
          if(lwritetrainpoints)                close(testpxunit)
          if(luseforces.and.lwritetrainforces) close(testfxunit)
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif ! iswitch
      endif ! mpirank.eq.0
!!
!!==============================================
!! combine the partial rmse values of all processes here only to avoid double counting
!!==============================================
      call mpi_allreduce(mpi_in_place,rmse_short,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,rmse_charge,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,rmse_totalcharge,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,rmse_elec,1,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,mad_short,1,&
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
!!==============================================
!! calculate the final RMSEs
!!==============================================
      call getrmse(ntrain,ncharges_sum,nenergies,netot,&
        nforces_short,nforces_elec,nforces_total,&
        rmse_short,rmse_charge,rmse_totalcharge,rmse_elec,rmse_etot,&
        rmse_force_s,rmse_force_e,rmse_force_t,&
        mad_short,mad_charge,mad_totalcharge,mad_elec,mad_etot,&
        mad_force_s,mad_force_e,mad_force_t)
!!
      return
      end
