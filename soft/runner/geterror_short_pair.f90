!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine geterror_short_pair(iswitch,countepoch,ntrain,&
        imaxerror_eshort,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        rmse_short,rmse_force_s,mad_short,&
        mad_force_s,maxerroreshort)
!!
      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions 
      use fittingoptions 
      use symfunctions
      use nnshort_pair
      use structures
!!
      implicit none
!!
      integer i1,i2,i3                                           ! internal
      integer imaxerror_eshort                                   ! out
      integer iswitch                                            ! in 0=train, 1=test
      integer nforces                                            ! internal
      integer ntrain                                             ! in
      integer ndone                                              ! internal
      integer ndonepara                                          ! internal
      integer ncount                                             ! internal
      integer npoints                                            ! internal
      integer nstruct                                            ! internal
      integer n_start                                            ! internal
      integer n_end                                              ! internal
      integer nenergies                                          ! internal
      integer, dimension(:), allocatable :: num_atoms_mpi        ! internal
      integer, dimension(:), allocatable :: num_pairs_mpi        ! internal
      integer, dimension(:,:), allocatable :: zelem_mpi          ! internal
      integer, dimension(:,:,:), allocatable :: zelemp_mpi       ! internal
      integer, dimension(:), allocatable :: imaxerrortemp        ! internal
      integer icount                                             ! internal
      integer countepoch                                         ! in
      integer tempunit                                           ! internal
      integer pxunit                                             ! internal
      integer fxunit                                             ! internal
!!
      real*8, dimension(:), allocatable :: shortenergy_mpi       ! internal
      real*8, dimension(:,:,:), allocatable :: xyzstruct_mpi     ! internal
      real*8, dimension(:,:,:), allocatable :: symfunctionp_mpi  ! internal
      real*8, dimension(:), allocatable ::  nneshort_mpi         ! internal
      real*8, dimension(:,:,:), allocatable :: nnshortforce_mpi  ! internal
      real*8, dimension(:), allocatable :: maxerrortemp          ! internal
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)    ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)    ! in
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)     ! in
      real*8 rmse_short                                          ! out
      real*8 rmse_force_s                                        ! out
      real*8 mad_short                                           ! out
      real*8 mad_force_s                                         ! out
      real*8 maxerroreshort                                      ! out
      real*8 nneshort_list(nblock)                               ! internal
      real*8 nnshortforce_list(3,max_num_atoms,nblock)           ! internal
      real*8 forcesum(3)                                         ! internal
!!
      character*25 filename                                      ! internal
!!
      logical, dimension(:), allocatable :: lperiodic_mpi        ! internal
!!
!!============================================================
!! initializations
!!============================================================
      nforces                  = 0
      ncount                   = ntrain
      ndone                    = 0
      nneshort_list(:)         = 0.0d0
      shortforce_list(:,:,:)   = 0.0d0
      nnshortforce_list(:,:,:) = 0.0d0
      nenergies                = 0
!!
!!============================================================
!! open files and write headers
!!============================================================
      if(mpirank.eq.0)then
        if(iswitch.eq.0)then
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit)
          if(lwritetrainpoints)then
            filename='trainpoints.000000.out'               ! File Name
            if(countepoch.gt.9999)then
              write(ounit,*)'Error: too many epochs in geterrorpair'
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
            write(trainpxunit,'(a)')' point      E_short(DFT)      E_short(NN)  '
          endif
          open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
          rewind(trainstructunit)
          if(luseforces)then
            if(lwritetrainforces)then
              filename='trainforces.000000.out'            ! File Name
              if(countepoch.gt.9999)then
                write(ounit,*)'Error: too many epochs in geterrorpair'
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
              write(trainfxunit,'(2a)')' point   atom       F_s(DFT)       F_s(NN)',&
              '         Delta          '
            endif
            open(trainfunit,file='trainforces.data',form='formatted',status='old')
            rewind(trainfunit)
          endif !' luseforces
!!
        elseif(iswitch.eq.1)then ! test error
          open(tymunit,file='testing.data',form='formatted',status='old')
          rewind(tymunit)
          if(lwritetrainpoints)then
            filename='testpoints.000000.out'          ! File Name
            if(countepoch.gt.999)then
              write(filename(14:17),'(i4)')countepoch
            elseif(countepoch.gt.99)then
              write(filename(15:17),'(i3)')countepoch
            elseif(countepoch.gt.9)then
              write(filename(16:17),'(i2)')countepoch
            else
              write(filename(17:17),'(i1)')countepoch
            endif
            open(testpxunit,file=filename,form='formatted',status='replace')
            rewind(testpxunit) !'
            write(testpxunit,'(a)')' point      E_short(DFT)      E_short(NN)  '
          endif
          open(teststructunit,file='teststruct.data',form='formatted',status='old')
          rewind(teststructunit)
          if(luseforces)then
            if(lwritetrainforces)then
              filename='testforces.000000.out'       ! File Name
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
              write(testfxunit,'(2a)')' point   atom       F_s(DFT)       F_s(NN)',&
              '         Delta          '
            endif
            open(testfunit,file='testforces.data',form='formatted',status='old')
            rewind(testfunit)
          endif !' luseforces
        endif ! iswitch
      endif ! mpirank.eq.0
!!
!!============================================================
!! loop block-wise over structures
!!============================================================
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!!============================================================
!! determine which nstruct structures of this block should be calculated by this process
!!============================================================
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!============================================================
!! do all file reading for training error here at one place to allow for parallelization
!!============================================================
!!
      if(mpirank.eq.0)then
!!============================================================
!! read the short range symmetry functions 
!!============================================================
        if(iswitch.eq.0)then ! train
          tempunit=symunit
        elseif(iswitch.eq.1)then ! test
          tempunit=tymunit
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
        call readfunctions(2,tempunit,npoints,npairs,&
          max_num_pairs,maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
          symfunction_short_pair_list)
!!
!!============================================================
!! read the structures needed for the calculation of the forces
!! must be called after readfunctions because it needs num_atoms_list
!!============================================================
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
!!============================================================
!! read short range forces from trainforces.data or testforces.data => shortforce_list
!!============================================================
        if(luseforces.and.lshort)then
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
      endif ! mpirank.eq.0
!!
!!============================================================
!! distribute the data to all processes
!!============================================================
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelemp_list,nblock*max_num_pairs,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_short_pair_list,nblock*max_num_pairs*maxnum_funcvalues_short_pair,&
        mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortforce_list,nblock*max_num_atoms*3,&
        mpi_real8,0,mpi_comm_world,mpierror)
!!
!!============================================================
!! end of file reading for training error
!!============================================================
!!
!!============================================================
!! allocate local arrays for this process
!!============================================================
      allocate(num_atoms_mpi(nstruct))
      allocate(num_pairs_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
      allocate(zelemp_mpi(2,nstruct,max_num_pairs))
      allocate(shortenergy_mpi(nstruct))
      allocate(xyzstruct_mpi(3,max_num_atoms,nstruct))
      allocate(lperiodic_mpi(nstruct))
      allocate(symfunctionp_mpi(maxnum_funcvalues_short_pair,max_num_pairs,nstruct))
      allocate(nneshort_mpi(nstruct))
      nneshort_mpi(:)=0.0d0
      allocate(nnshortforce_mpi(3,max_num_atoms,nstruct))
      nnshortforce_mpi(:,:,:)=0.0d0
!!
!!============================================================
!! fill local arrays for this process
!!============================================================
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        num_pairs_mpi(icount)       = num_pairs_list(i1)
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        zelemp_mpi(:,icount,:)      = zelemp_list(:,i1,:)
        shortenergy_mpi(icount)     = shortenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
      enddo ! i1
!!
!!============================================================
!! get the short range energies 
!!============================================================
       ndonepara=ndone+n_start-1
       call getshortenergies_parapair(nstruct,ndonepara,&
          nenergies,imaxerror_eshort,&
          num_atoms_mpi,num_pairs_mpi,&
          zelem_mpi,zelemp_mpi,&
          minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
          symfunction_short_pair_list(1,1,n_start),lattice_list(1,1,n_start),&
          nnshortforce_mpi,xyzstruct_mpi,nneshort_mpi,&
          rmse_short,mad_short,maxerroreshort,shortenergy_mpi,&
          lperiodic_mpi)
!!
!!============================================================
!! for parallel case: find imaxerror_eshort and maxerroreshort from all processes
!!============================================================
        if(mpisize.gt.1)then
          allocate(imaxerrortemp(mpisize))
          imaxerrortemp(:)=0
          allocate(maxerrortemp(mpisize))
          maxerrortemp(:) =0.0d0
          imaxerrortemp(mpirank+1)=imaxerror_eshort
          maxerrortemp(mpirank+1) =maxerroreshort
          call mpi_allreduce(mpi_in_place,imaxerrortemp,mpisize,&
            mpi_integer,mpi_sum,mpi_comm_world,mpierror)
          call mpi_allreduce(mpi_in_place,maxerrortemp,mpisize,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
          do i1=1,mpisize
            if(maxerrortemp(i1).gt.maxerroreshort)then
              maxerroreshort = maxerrortemp(i1)
              imaxerror_eshort= imaxerrortemp(i1)
            endif
          enddo ! i1
          deallocate(imaxerrortemp)
          deallocate(maxerrortemp)
        endif
!!
!!============================================================
!! copy results back on full array
!!============================================================
        icount=0
        nneshort_list(:)         =0.0d0
        nnshortforce_list(:,:,:) =0.0d0
        do i1=n_start,n_end
          icount=icount+1
          nneshort_list(i1)         = nneshort_mpi(icount)
          nnshortforce_list(:,:,i1) = nnshortforce_mpi(:,:,icount)
        enddo
!!
!!============================================================
!! distribute results to all processes
!!============================================================
        call mpi_allreduce(nneshort_list,nneshort_list,nblock,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(nnshortforce_list,nnshortforce_list,nblock*max_num_atoms*3,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!!============================================================
!! check NN forces if requested
!!============================================================
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
!!============================================================
!! calculate the RMSE of the short range forces 
!!============================================================
      if(luseforces)then
        call calcrmse_forces(npoints,&
          nforces,rmse_force_s,mad_force_s,&
          shortforce_list,nnshortforce_list,maxforce)
      endif
!!
!!============================================================
!! set file units depending on training or test set 
!!============================================================
      if(iswitch.eq.0)then ! train set
        pxunit=trainpxunit
        fxunit=trainfxunit
      elseif(iswitch.eq.1)then ! test set
        pxunit=testpxunit
        fxunit=testfxunit
      else
        write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
        stop
      endif
!!
!!============================================================
!! if requested write energy details for all points to trainpoints.out
!!============================================================
      if(mpirank.eq.0)then
        if(lwritetrainpoints)then
          do i1=1,npoints
            write(pxunit,'(i6,x,2f14.8)')i1+ndone,&
            shortenergy_list(i1),nneshort_list(i1)
          enddo
        endif
!!
!!============================================================
!! if requested write force details for all points to trainforces.out
!!============================================================
        if(luseforces.and.lwritetrainforces)then
          do i1=1,npoints
            do i2=1,num_atoms_list(i1)
            write(fxunit,'(2i5,a4,2f14.6,2x,f14.8)')&
              i1+ndone,i2,' Fx ',&
              shortforce_list(1,i2,i1),nnshortforce_list(1,i2,i1),&
              shortforce_list(1,i2,i1)-nnshortforce_list(1,i2,i1)
            write(fxunit,'(2i5,a4,2f14.6,2x,f14.8)')&
              i1+ndone,i2,' Fy ',&
              shortforce_list(2,i2,i1),nnshortforce_list(2,i2,i1),&
              shortforce_list(2,i2,i1)-nnshortforce_list(2,i2,i1)
            write(fxunit,'(2i5,a4,2f14.6,2x,f14.8)')&
              i1+ndone,i2,' Fz ',&
              shortforce_list(3,i2,i1),nnshortforce_list(3,i2,i1),&
              shortforce_list(3,i2,i1)-nnshortforce_list(3,i2,i1)
            enddo ! i2
          enddo ! i1
        endif ! luseforces
!!
      endif ! mpirank.eq.0
!!
      ndone=ndone+npoints
!!
!!============================================================
!! deallocate local arrays 
!!============================================================
      deallocate(num_pairs_mpi)
      deallocate(num_atoms_mpi)
      deallocate(zelem_mpi)
      deallocate(zelemp_mpi)
      deallocate(shortenergy_mpi)
      deallocate(symfunctionp_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(lperiodic_mpi)
      deallocate(nneshort_mpi)
      deallocate(nnshortforce_mpi)
!!
!!============================================================
!! if there are structures left go to next block of structures 
!!============================================================
      if(ncount.gt.0) goto 10
!!============================================================
!! end of block wise loop 
!!============================================================
!!
!!============================================================
!! close files
!!============================================================
      if(mpirank.eq.0)then
        if(iswitch.eq.0)then
          close(symunit)
          close(trainstructunit)
          if(luseforces)                       close(trainfunit)
          if(lwritetrainpoints)                close(trainpxunit)
          if(luseforces.and.lwritetrainforces) close(trainfxunit)
        elseif(iswitch.eq.1)then
          close(tymunit)
          close(teststructunit)
          if(luseforces)                       close(testfunit)
          if(lwritetrainpoints)                close(testpxunit)
          if(luseforces.and.lwritetrainforces) close(testfxunit)
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif ! iswitch
      endif ! mpirank.eq.0
!!
!!============================================================
!! combine the partial rmse values of all processes here only to avoid double counting
!!============================================================
      call mpi_allreduce(mpi_in_place,rmse_short,1,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,mad_short,1,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,nenergies,1,mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
!!============================================================
!! calculate the final RMSEs for the training set
!!============================================================
      call getrmse_short(nenergies,&
        nforces,rmse_short,rmse_force_s,mad_short,mad_force_s)
!!
      return
      end
