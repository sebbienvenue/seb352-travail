!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! this subroutine is obsolete and should be merged with getkalmanmatrices

!! called by:
!! - fitting.f90
!! - fitting_batch.f90
!!
      subroutine getkalmanmatricesunformatted(ndim,&
         iseed,kseed,lseed,mseed,&
         kaldim,kaledim,corrdim,corrfdim,corredim,corrcdim,&
         maxcorrdim,maxcorrfdim,maxcorredim,&
         num_weightsshortfree,num_weightsewaldfree,&
         corrmatrix_list,corrmatrixf_list,corrmatrixe_list,corrmatrixc,&
         kalmanlambda_local,&
         lshort,lewald)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3,i4         ! internal
      integer ndim               ! in
      integer num_weightsshortfree(ndim)    ! in
      integer num_weightsewaldfree(nelem)    ! in
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
      integer corredim(nelem)            ! in
      integer corrcdim                   ! in
      integer kaldim(ndim)              ! in
      integer kaledim(nelem)             ! in
      integer maxcorrdim                 ! in
      integer maxcorrfdim                ! in
      integer maxcorredim                ! in
      integer isum                       ! internal
      integer iseed
      integer kseed
      integer lseed
      integer mseed
!!
      real*8 corrmatrix_list(maxcorrdim,ndim)                       ! out
      real*8 corrmatrixf_list(maxcorrfdim,ndim)                     ! out
      real*8 corrmatrixe_list(maxcorredim,nelem)                     ! out
      real*8 corrmatrixc(corrcdim)                                   ! out
      real*8 kalmanlambda_local(ndim)                                ! out
      real*8, dimension(:)  , allocatable :: corrmatrix_full         ! internal 
      real*8, dimension(:)  , allocatable :: corrmatrixtemp_list
      real*8, dimension(:)  , allocatable :: corrmatrixfull_list
      real*8, dimension(:)  , allocatable :: corrmatrixtemp
!!
      character*21 kalfilenames(ndim)  ! internal
      character*21 ctemp21              ! internal
      character*20 kalfilenamee(nelem)  ! internal
      character*20 ctemp20              ! internal
!!
      logical lshort            ! in
      logical lewald            ! in
      logical lexist            ! internal
!!
!! initializations
      isum=0
      kalfilenames(:)='kalman.short.000.data'
      kalfilenamee(:)='kalman.elec.000.data'
      do i1=1,ndim
        ctemp21=kalfilenames(i1)
        if(nucelem(i1).gt.99)then
          write(ctemp21(14:16),'(i3)')nucelem(i1)
        elseif(nucelem(i1).gt.9)then
          write(ctemp21(15:16),'(i2)')nucelem(i1)
        else
          write(ctemp21(16:16),'(i1)')nucelem(i1)
        endif
        kalfilenames(i1)=ctemp21
      enddo
      do i1=1,nelem
        ctemp20=kalfilenamee(i1)
        if(nucelem(i1).gt.99)then
          write(ctemp20(13:15),'(i3)')nucelem(i1)
        elseif(nucelem(i1).gt.9)then
          write(ctemp20(14:15),'(i2)')nucelem(i1)
        else
          write(ctemp20(15:15),'(i1)')nucelem(i1)
        endif
        kalfilenamee(i1)=ctemp20
      enddo
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!! initialization of the correlation matrix for the Kalman filter
      if((mpisize.eq.1).or.lompmkl)then
        if(lewald.and.lchargeconstraint)then
          do i1=1,nelem
            isum=isum+num_weightsewaldfree(i1)
          enddo
        endif
        do i1=1,ndim
          if(lshort.and.(optmodee.eq.1))then
            call initialcorrmatrix(num_weightsshortfree(i1),corrmatrix_list(1,i1))
          endif
          if(lshort.and.(optmodef.eq.1))then
            call initialcorrmatrix(num_weightsshortfree(i1),corrmatrixf_list(1,i1))
          endif
          if(lewald.and.(optmodeq.eq.1))then
            call initialcorrmatrix(num_weightsewaldfree(i1),corrmatrixe_list(1,i1))
            if(lchargeconstraint)then
              call initialcorrmatrix(isum,corrmatrixc)
            endif
          endif
        enddo
!!
      else ! parallel case
!! short range part for energy fitting 
        if(lshort.and.(optmodee.eq.1))then
          do i1=1,ndim
            allocate(corrmatrix_full(num_weightsshortfree(i1)*num_weightsshortfree(i1)))
            corrmatrix_full(:)=0.0d0
!! set diagonal elements to 1.0d0
            do i2=1,num_weightsshortfree(i1)
              do i3=1,num_weightsshortfree(i1)
                if(i2.eq.i3)then
                  itemp=(i2-1)*num_weightsshortfree(i1)+i3
                  corrmatrix_full(itemp)=1.0d0
                endif
              enddo ! i3
            enddo ! i2
!! distribute correlation matrix to processes, distribute only full columns
            call mpifitdistribution(num_weightsshortfree(i1),&
              itemp,n_start,n_end)
            do i2=1,corrdim(i1)
              corrmatrix_list(i2,i1)&
                =corrmatrix_full(i2+(n_start-1)*num_weightsshortfree(i1))
            enddo ! i2
            deallocate(corrmatrix_full)
          enddo ! i1
        endif ! lshort
!! short range part for force fitting 
        if(lshort.and.(optmodef.eq.1))then
          do i1=1,ndim
            allocate(corrmatrix_full(num_weightsshortfree(i1)*num_weightsshortfree(i1)))
            corrmatrix_full(:)=0.0d0
!! set diagonal elements to 1.0d0
            do i2=1,num_weightsshortfree(i1)
              do i3=1,num_weightsshortfree(i1)
                if(i2.eq.i3)then
                  itemp=(i2-1)*num_weightsshortfree(i1)+i3
                  corrmatrix_full(itemp)=1.0d0
                endif
              enddo ! i3
            enddo ! i2
!! distribute correlation matrix to processes, distribute only full columns
            call mpifitdistribution(num_weightsshortfree(i1),&
              itemp,n_start,n_end)
            do i2=1,corrfdim(i1)
              corrmatrixf_list(i2,i1)&
                =corrmatrix_full(i2+(n_start-1)*num_weightsshortfree(i1))
            enddo ! i2
            deallocate(corrmatrix_full)
          enddo ! i1
        endif ! lshort
!!
!! electrostatic part
        if(lewald.and.lchargeconstraint)then
          isum=0
          do i1=1,nelem
            isum=isum+num_weightsewaldfree(i1)
          enddo
        endif
        if(lewald.and.(optmodeq.eq.1))then
          do i1=1,nelem
            allocate(corrmatrix_full(num_weightsewaldfree(i1)*num_weightsewaldfree(i1)))
            corrmatrix_full(:)=0.0d0
!! set diagonal elements to 1.0d0
            do i2=1,num_weightsewaldfree(i1)
              do i3=1,num_weightsewaldfree(i1)
                if(i2.eq.i3)then
                  itemp=(i2-1)*num_weightsewaldfree(i1)+i3
                  corrmatrix_full(itemp)=1.0d0
                endif
              enddo ! i3
            enddo ! i2
!! distribute correlation matrix to processes, distribute only full columns
            call mpifitdistribution(num_weightsewaldfree(i1),&
              itemp,n_start,n_end)
            do i2=1,corredim(i1)
              corrmatrixe_list(i2,i1)&
                =corrmatrix_full(i2+(n_start-1)*num_weightsewaldfree(i1))
            enddo ! i2
            deallocate(corrmatrix_full)
          enddo ! i1
!! FIXME: at the moment we have the full corrmatrixc on all processes
          if(lchargeconstraint)then
            call initialcorrmatrix(isum,corrmatrixc)
          endif ! lchargeconstraint
        endif ! lewald
      endif ! mpisize
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read Kalman filter data if requested
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(lrestkalman)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read Kalman data for short range NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(lshort.and.((optmodee.eq.1).or.(optmodef.eq.1)))then
          do i4=1,ndim
            inquire(file=kalfilenames(i4),exist=lexist)
            if(.not.lexist)then
              write(ounit,*)'Error: file not found ',kalfilenames(i4)
              stop
            endif
            if((mpisize.eq.1).or.lompmkl)then
              open(kalunit,file=kalfilenames(i4),form='unformatted',status='old')
                rewind(kalunit) !'
                read(kalunit)kalmanlambda_local(i4)
                do i1=1,num_weightsshortfree(i4)*(num_weightsshortfree(i4)+1)/2
                  read(kalunit)corrmatrix_list(i1,i4)
                enddo
                read(kalunit)iseed
                read(kalunit)kseed
                read(kalunit)lseed
                read(kalunit)mseed
              close(kalunit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            else ! parallel case
              if(mpirank.eq.0)then
                open(kalunit,file=kalfilenames(i4),form='unformatted',status='old')
                read(kalunit)kalmanlambda_local(i4) 
              endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read and distribute correlation matrices for short range weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if(mpirank.eq.0)then
                allocate(corrmatrixfull_list(num_weightsshortfree(i4)*num_weightsshortfree(i4)))
                allocate(corrmatrixtemp_list(num_weightsshortfree(i4)*(num_weightsshortfree(i4)+1)/2))
                allocate(corrmatrixtemp(num_weightsshortfree(i4)*num_weightsshortfree(i4)))
!!' read packed correlation matrices
                do i1=1,num_weightsshortfree(i4)*(num_weightsshortfree(i4)+1)/2
                  read(kalunit)corrmatrixtemp_list(i1)
                enddo
                read(kalunit)iseed
                read(kalunit)kseed
                read(kalunit)lseed
                read(kalunit)mseed
!! unpack correlation matrices
                icount=0
                jcount=0
                do i1=1,num_weightsshortfree(i4)
                  do i2=1,num_weightsshortfree(i4)
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
                do i1=1,num_weightsshortfree(i4)
                  do i2=1,num_weightsshortfree(i4)
                    jcount=jcount+1
                    if(i2.lt.i1)then
                      icount=(i2-1)*num_weightsshortfree(i4)+i1
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
              kaldim_list(:)=kaldim_list(:)*num_weightsshortfree(i4)
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
          enddo ! i4=1,ndim
          call mpi_bcast(kalmanlambda_local,ndim,mpi_real8,0,mpi_comm_world,mpierror)
        endif ! lshort
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read Kalman data for electrostatic NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(lewald.and.(optmodeq.eq.1))then
          do i4=1,nelem
            inquire(file=kalfilenamee(i4),exist=lexist)
            if(.not.lexist)then
              write(ounit,*)'Error: file not found ',kalfilenamee(i4)
              stop
            endif
            if((mpisize.eq.1).or.lompmkl)then
              open(kaleunit,file=kalfilenamee(i4),form='unformatted',status='old')
              rewind(kaleunit) !'
              read(kaleunit)kalmanlambdae(i4)
              do i1=1,num_weightsewaldfree(i4)*(num_weightsewaldfree(i4)+1)/2
                read(kaleunit)corrmatrixe_list(i1,i4)
              enddo
              read(kaleunit)iseed
              read(kaleunit)kseed
              read(kaleunit)lseed
              read(kaleunit)mseed
              close(kaleunit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            else ! parallel case
              if(mpirank.eq.0)then
                open(kaleunit,file=kalfilenamee(i4),form='unformatted',status='old')
                read(kaleunit)kalmanlambdae(i4) !'
              endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read and distribute correlatione matrices for electrostatic weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if(mpirank.eq.0)then
                allocate(corrmatrixfull_list(num_weightsewaldfree(i4)*num_weightsewaldfree(i4)))
                allocate(corrmatrixtemp_list(num_weightsewaldfree(i4)*(num_weightsewaldfree(i4)+1)/2))
                allocate(corrmatrixtemp(num_weightsewaldfree(i4)*num_weightsewaldfree(i4)))
!!' read packed correlatione matrices
                do i1=1,num_weightsewaldfree(i4)*(num_weightsewaldfree(i4)+1)/2
                  read(kaleunit)corrmatrixtemp_list(i1)
                enddo
                read(kaleunit)iseed
                read(kaleunit)kseed
                read(kaleunit)lseed
                read(kaleunit)mseed
!! unpack correlation matrices
                icount=0
                jcount=0
                do i1=1,num_weightsewaldfree(i4)
                  do i2=1,num_weightsewaldfree(i4)
                    jcount=jcount+1
                    if(i2.ge.i1)then
                      icount=icount+1
                      corrmatrixfull_list(jcount)=corrmatrixtemp_list(icount)
                   endif
                  enddo ! i2
                enddo ! i1
                deallocate(corrmatrixtemp_list)
!! fill also the upper triangle
                icount=0
                jcount=0
                do i1=1,num_weightsewaldfree(i4)
                  do i2=1,num_weightsewaldfree(i4)
                    jcount=jcount+1
                    if(i2.lt.i1)then
                      icount=(i2-1)*num_weightsewaldfree(i4)+i1
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
              call mpi_gatherv(kaledim(i4),1,mpi_integer,&
                kaldim_list,recvcnts,displs,mpi_integer,0,mpi_comm_world,mpierror)
              kaldim_list(:)=kaldim_list(:)*num_weightsewaldfree(i4)
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
                corrmatrixe_list(1,i4),kaldim_list(mpirank+1),mpi_real8,&
                0,mpi_comm_world,mpierror)
              deallocate(kaldim_list)
              deallocate(displs)
              deallocate(recvcnts)
              if(mpirank.eq.0)then
                deallocate(corrmatrixfull_list)
                deallocate(corrmatrixtemp)
              endif ! mpirank.eq.0
              if(mpirank.eq.0)then
                close(kaleunit) 
              endif ! mpirank.eq.0
            endif ! mpisize.eq.1
          enddo ! i4=1,nelem
          call mpi_bcast(kalmanlambdae,nelem,mpi_real8,0,mpi_comm_world,mpierror)
        endif ! lewald
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read matrices for charge constraint
!! FIXME: so far all processes have the full constraint matrices
        if(lewald.and.(optmodeq.eq.1))then
          if(lchargeconstraint)then
            isum=0
            do i4=1,nelem
              isum=isum+num_weightsewaldfree(i4)
            enddo
            if(mpirank.eq.0)then
              inquire(file='kalmanc.data',exist=lexist)
              if(.not.lexist)then
                write(ounit,*)'kalmanc.data file not found'
                stop
              endif
              open(kalcunit,file='kalmanc.data',form='unformatted',status='old')
                rewind(kalcunit) !'
                read(kalcunit)kalmanlambdac
                do i1=1,isum*(isum+1)/2
                  read(kalcunit)corrmatrixc(i1)
                enddo
              close(kalcunit)
            endif ! mpirank.eq.0
!! distribute kalmanlambdac, corrmatrixc 
            call mpi_bcast(kalmanlambdac,1,&
              mpi_real8,0,mpi_comm_world,mpierror)
            call mpi_bcast(corrmatrixc,&
              isum*(isum+1)/2,&
              mpi_real8,0,mpi_comm_world,mpierror)
          endif ! lchargeconstraint
        endif ! lewald
!!
      write(ounit,*)'WARNING: reading correlation matrix for force fitting from file is not implemented yet'
      endif ! lrestkalman
!!
      return
      end
