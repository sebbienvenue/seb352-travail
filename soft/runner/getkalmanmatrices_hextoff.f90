!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine
!! this is not yet generalized to include pair case
!! this subroutine should be generalized to replace also getkalmanmatricesunformatted

!! called by:
!!
      subroutine getkalmanmatrices_hextoff(ndim,&
         iseed,kseed,lseed,mseed,&
         corrdim,&
         maxcorrdim,num_weightsfree_local,&
         corrmatrix_list,&
         kalmanlambda_local)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer i0,i1,i2,i3,i4         ! internal
      integer ndim               ! in
      integer num_weightsfree_local(ndim)    ! in
      integer itemp               ! internal
      integer n_start             ! internal
      integer n_end               ! internal
      integer icount
      integer jcount
      integer, dimension(:)  , allocatable :: kaldim_list
      integer, dimension(:)  , allocatable :: displs
      integer, dimension(:)  , allocatable :: recvcnts
!! Kalman matrix dimensions:
      integer corrdim(ndim)             ! in
      integer corrfdim(ndim)            ! in
      integer kaldim(ndim)              ! in
      integer maxcorrdim                 ! in
      integer maxcorrfdim                ! in
      integer isum                       ! internal
      integer iseed
      integer kseed
      integer lseed
      integer mseed
!!
      real*8 corrmatrix_list(maxcorrdim,ndim)                       ! out
      real*8 kalmanlambda_local(ndim)                                ! out
      real*8, dimension(:)  , allocatable :: corrmatrix_full         ! internal 
      real*8, dimension(:)  , allocatable :: corrmatrixtemp_list
      real*8, dimension(:)  , allocatable :: corrmatrixfull_list
      real*8, dimension(:)  , allocatable :: corrmatrixtemp
!!
      character*21 kalfilenames(ndim)   ! internal
      character*21 ctemp21              ! internal
      character*24 kalfilenamep(ndim)   ! internal
      character*24 ctemp24              ! internal
      character*40 kalfilenamehextoff(ndim)   ! internal
!!
      logical lexist            ! internal
!!
!! initializations
      isum=0
      i0 = tripletindex(hextoff_training_triplet(1),hextoff_training_triplet(2),hextoff_training_triplet(3))
      write(kalfilenamehextoff(1),'(A,I3.3,A,I3.3,A,I3.3,A)') &
          'kalman.hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
!!
!! initialization of the correlation matrix for the Kalman filter
      if((mpisize.eq.1).or.lompmkl)then
        if(optmodehextoff.eq.1)then
          call initialcorrmatrix(num_weightsfree_local(1),corrmatrix_list(1,1))
        endif
      else ! parallel case
!! Hextoff element fitting
        if(optmodehextoff.eq.1)then
            allocate(corrmatrix_full(num_weightsfree_local(1)*num_weightsfree_local(1)))
            corrmatrix_full(:)=0.0d0
!! set diagonal elements to 1.0d0
            do i2=1,num_weightsfree_local(1)
              do i3=1,num_weightsfree_local(1)
                if(i2.eq.i3)then
                  itemp=(i2-1)*num_weightsfree_local(1)+i3
                  corrmatrix_full(itemp)=1.0d0
                endif
              enddo ! i3
            enddo ! i2
!! distribute correlation matrix to processes, distribute only full columns
            call mpifitdistribution(num_weightsfree_local(1),itemp,n_start,n_end)
            do i2=1,corrdim(1)
              corrmatrix_list(i2,1)&
                =corrmatrix_full(i2+(n_start-1)*num_weightsfree_local(1))
            enddo ! i2
            deallocate(corrmatrix_full)
        endif !
      endif ! mpisize
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read Kalman filter data if requested
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(lrestkalman)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read Kalman data for hextoff NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((optmodehextoff.eq.1))then
!!          do i4,1,ndim
            i4 = 1
            inquire(file=kalfilenamehextoff(i4),exist=lexist)
            if(.not.lexist)then
              write(ounit,*)'Error: file not found ',kalfilenames(i4)
              stop
            endif
            if((mpisize.eq.1).or.lompmkl)then
                if(lnntb.and.nntb_flag(3))then
                  open(kalunit,file=kalfilenamehextoff(i4),form='formatted',status='old')
                else
                  write(ounit,*)'ERROR: unknown nntb type in writekalman ', nntb_flag
                  stop
                endif  
              rewind(kalunit) !'

              if(lreadunformatted)then
                read(kalunit)kalmanlambda_local(i4)
              else
                read(kalunit,*)kalmanlambda_local(i4)
              endif
              do i1=1,num_weightsfree_local(i4)*(num_weightsfree_local(i4)+1)/2
                if(lreadunformatted)then
                  read(kalunit)corrmatrix_list(i1,i4)
                else
                  read(kalunit,*)corrmatrix_list(i1,i4)
                endif
              enddo
              if(lreadunformatted)then
                read(kalunit)iseed
                read(kalunit)kseed
                read(kalunit)lseed
                read(kalunit)mseed
              else
                read(kalunit,*)iseed
                read(kalunit,*)kseed
                read(kalunit,*)lseed
                read(kalunit,*)mseed
              endif
              close(kalunit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            else ! parallel case
              if(mpirank.eq.0)then
                if(lnntb.and.(nntb_flag(3)))then
                  open(kalunit,file=kalfilenames(i4),form='formatted',status='old')
                elseif(nn_type_short.eq.2)then
                  open(kalunit,file=kalfilenamep(i4),form='formatted',status='old')
                else
                  write(ounit,*)'ERROR: unknown nntb type in writekalman ',nntb_flag
                  stop
                endif  
                rewind(kalunit)
                if(lreadunformatted)then
                  read(kalunit)kalmanlambda_local(i4) 
                else
                  read(kalunit,*)kalmanlambda_local(i4) 
                endif
              endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read and distribute correlation matrices for short range weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if(mpirank.eq.0)then
                allocate(corrmatrixfull_list(num_weightsfree_local(i4)*num_weightsfree_local(i4)))
                allocate(corrmatrixtemp_list(num_weightsfree_local(i4)*(num_weightsfree_local(i4)+1)/2))
                allocate(corrmatrixtemp(num_weightsfree_local(i4)*num_weightsfree_local(i4)))
!!' read packed correlation matrices
                do i1=1,num_weightsfree_local(i4)*(num_weightsfree_local(i4)+1)/2
                  if(lreadunformatted)then
                    read(kalunit)corrmatrixtemp_list(i1)
                  else
                    read(kalunit,*)corrmatrixtemp_list(i1)
                  endif
                enddo
                if(lreadunformatted)then
                  read(kalunit)iseed
                  read(kalunit)kseed
                  read(kalunit)lseed
                  read(kalunit)mseed
                else
                  read(kalunit,*)iseed
                  read(kalunit,*)kseed
                  read(kalunit,*)lseed
                  read(kalunit,*)mseed
                endif
!! unpack correlation matrices
                icount=0
                jcount=0
                do i1=1,num_weightsfree_local(i4)
                  do i2=1,num_weightsfree_local(i4)
                    jcount=jcount+1
                    if(i2.ge.i1)then
                      icount=icount+1
                      corrmatrixfull_list(jcount)=corrmatrixtemp_list(icount)
                    endif
                  enddo ! i2
                enddo ! i1
                deallocate(corrmatrixtemp_list)
!! fill alsot the upper triangle
                icount=0
                jcount=0
                do i1=1,num_weightsfree_local(i4)
                  do i2=1,num_weightsfree_local(i4)
                    jcount=jcount+1
                    if(i2.lt.i1)then
                      icount=(i2-1)*num_weightsfree_local(i4)+i1
                      corrmatrixfull_list(jcount)=corrmatrixfull_list(icount)
                    endif
                  enddo ! i2
                enddo ! i1
              endif ! mpirank.eq.0
              call mpi_barrier(mpi_comm_world,mpierror)
!! determine the number of elements for each process
              allocate(kaldim_list(mpisize))
              kaldim_list(:)=0
              allocate(displs(mpisize))
              displs(:)=0
              allocate(recvcnts(mpisize))
              recvcnts(:)=0
!! determine the displacements for kaldim_list
              do i1=1,mpisize
                displs(i1)=i1-1
              enddo
              recvcnts(:)=1
!! determine the corrmatrix fragment size of all processes
              call mpi_gatherv(kaldim(i4),1,mpi_integer,&
                kaldim_list,recvcnts,displs,mpi_integer,0,mpi_comm_world,mpierror)
              kaldim_list(:)=kaldim_list(:)*num_weightsfree_local(i4)
              call mpi_bcast(kaldim_list,mpisize,mpi_integer,0,mpi_comm_world,mpierror)
! determine the displacements for corrmatrixtemp
              displs(:)=0
              do i1=2,mpisize
                do i2=1,i1-1
                  displs(i1)=displs(i1)+kaldim_list(i2)
                enddo ! i2
              enddo ! i1
!! distribute parts of the correlation matrix to all processesi
              if(mpirank.eq.0)then
                corrmatrixtemp(:)=corrmatrixfull_list(:)
              endif ! mpirank.eq.0
              call mpi_barrier(mpi_comm_world,mpierror)
              call mpi_scatterv(corrmatrixtemp,&
                kaldim_list,displs,mpi_real8,&
                corrmatrix_list(1,i4),kaldim_list(mpirank+1),mpi_real8,&
                0,mpi_comm_world,mpierror)
              deallocate(kaldim_list)
              deallocate(displs)
              deallocate(recvcnts)
              if(mpirank.eq.0)then
                deallocate(corrmatrixfull_list)
                deallocate(corrmatrixtemp)
              endif ! mpirank.eq.0
              if(mpirank.eq.0)then
                close(kalunit) 
              endif ! mpirank.eq.0
            endif ! mpisize.eq.1
!!          enddo ! i4=1,ndim
          call mpi_bcast(kalmanlambda_local,ndim,mpi_real8,0,mpi_comm_world,mpierror)
        endif ! 
!!
      write(ounit,*)'WARNING: reading correlation matrix for force fitting from file is not implemented yet'
!!
      endif ! lrestkalman

      return
      end
