!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! called by:
!!
      subroutine writekalman_short(ndim,&
        iseed,kseed,lseed,mseed,&
        maxcorrdim,kaldim,&
        num_weightsfree_local,&
        kalmanlambda_local,corrmatrix_list)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer ndim             ! in
      integer num_weightsfree_local(ndim)  ! in
      integer i1,i2,i4          ! internal
      integer icount            ! internal
      integer jcount            ! internal
      integer, dimension(:)  , allocatable :: kaldim_list 
      integer, dimension(:)  , allocatable :: displs 
      integer, dimension(:)  , allocatable :: recvcnts 
!! Kalman matrix dimensions:
      integer maxcorrdim           ! in
      integer kaldim(ndim)        ! in
      integer iseed
      integer kseed
      integer lseed
      integer mseed

      real*8 kalmanlambda_local(ndim)         ! in
!!
      real*8 corrmatrix_list(maxcorrdim,ndim)   ! in
      real*8, dimension(:)  , allocatable :: corrmatrixtemp 
      real*8, dimension(:)  , allocatable :: corrmatrixtemp_triangle 
!!
!!
      character*21 kalfilenames(ndim)  ! internal
      character*21 ctemp21              ! internal
      character*24 kalfilenamep(ndim)  ! internal
      character*24 ctemp24              ! internal
!!
!!
      write(ounit,*)'WARNING: correlation matrix for forces is not written to file'
!!'
!! initializations
      if(nn_type_short.eq.1)then
        kalfilenames(:)='kalman.short.000.data'
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
      elseif(nn_type_short.eq.2)then
        kalfilenamep(:)='kalman.pair.000.000.data'
        icount=0
        do i1=1,nelem
          do i2=i1,nelem
            icount=icount+1
            ctemp24=kalfilenamep(icount)
            if(nucelem(i1).gt.99)then
              write(ctemp24(13:15),'(i3)')nucelem(i1)
            elseif(nucelem(i1).gt.9)then
              write(ctemp24(14:15),'(i2)')nucelem(i1)
            else
              write(ctemp24(15:15),'(i1)')nucelem(i1)
            endif
            if(nucelem(i2).gt.99)then
              write(ctemp24(17:19),'(i3)')nucelem(i2)
            elseif(nucelem(i2).gt.9)then
              write(ctemp24(18:19),'(i2)')nucelem(i2)
            else
              write(ctemp24(19:19),'(i1)')nucelem(i2)
            endif
            kalfilenamep(icount)=ctemp24
          enddo
        enddo
      else
        write(ounit,*)'ERROR: unknown nn_type_short in writekalman ',nn_type_short
        stop
      endif ! nn_type_short
!!
      if(lsavekalman)then
!! write Kalman data for short range NN
          do i4=1,ndim
            if(mpirank.eq.0)then
              if(nn_type_short.eq.1)then
                open(kalunit,file=kalfilenames(i4),form='formatted',status='replace') 
              elseif(nn_type_short.eq.2)then
                open(kalunit,file=kalfilenamep(i4),form='formatted',status='replace') 
              else
                write(ounit,*)'ERROR: unknown nn_type_short in writekalman ',nn_type_short
                stop
              endif
              if(lwriteunformatted)then
                write(kalunit)kalmanlambda_local(i4) !'
              else
                write(kalunit,*)kalmanlambda_local(i4) !'
              endif
            endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write corrmatrix for short range weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for mpisize.eq.1 we do not need to compress
            if((mpisize.eq.1).or.lompmkl)then
              if(mpirank.eq.0)then
                do i1=1,num_weightsfree_local(i4)*(num_weightsfree_local(i4)+1)/2
                  if(lwriteunformatted)then
                    write(kalunit)corrmatrix_list(i1,i4)
                  else
                    write(kalunit,'(e20.12)')corrmatrix_list(i1,i4)
                  endif
                enddo
              endif ! mpirank.eq.0
              if(lwriteunformatted)then
                write(kalunit)iseed
                write(kalunit)kseed
                write(kalunit)lseed
                write(kalunit)mseed
              else
                write(kalunit,*)iseed
                write(kalunit,*)kseed
                write(kalunit,*)lseed
                write(kalunit,*)mseed
              endif
            else ! more than one core
!! collect the correlation matrices of all processes
              if(mpirank.eq.0)then
                allocate(corrmatrixtemp_triangle(num_weightsfree_local(i4)*(num_weightsfree_local(i4)+1)/2))
              endif ! mpirank.eq.0
              allocate(corrmatrixtemp(num_weightsfree_local(i4)*num_weightsfree_local(i4)))
              allocate(kaldim_list(mpisize))
              call mpi_barrier(mpi_comm_world,mpierror)
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
!!
              kaldim_list(:)=kaldim_list(:)*num_weightsfree_local(i4)
              call mpi_bcast(kaldim_list,mpisize,mpi_integer,0,mpi_comm_world,mpierror)
!! determine the displacements for corrmatrixtemp 
              displs(:)=0
              do i1=2,mpisize
                do i2=1,i1-1
                  displs(i1)=displs(i1)+kaldim_list(i2)
                enddo ! i2
              enddo ! i1
!! collect the corrmatrix fragments of all processes
              call mpi_gatherv(corrmatrix_list(1,i4),&
                kaldim_list(mpirank+1),mpi_real8,&
                corrmatrixtemp,kaldim_list,displs,mpi_real8,0,mpi_comm_world,mpierror)
!! store in full array
              call mpi_barrier(mpi_comm_world,mpierror)
              deallocate(kaldim_list)
              deallocate(displs)
              deallocate(recvcnts)
!! compress corrmatrix to lower triangle before saving'
              if(mpirank.eq.0)then
                icount=0
                jcount=0
                do i1=1,num_weightsfree_local(i4)
                  do i2=1,num_weightsfree_local(i4)
                    jcount=jcount+1
                    if(i2.ge.i1)then
                      icount=icount+1
                      corrmatrixtemp_triangle(icount)=corrmatrixtemp(jcount)
                    endif
                  enddo ! i2
                enddo ! i1
              endif ! mpirank.eq.0
              deallocate(corrmatrixtemp)
!! now save correlation matrix for short range weights
              if(mpirank.eq.0)then
                do i1=1,num_weightsfree_local(i4)*(num_weightsfree_local(i4)+1)/2
                  if(lwriteunformatted)then
                    write(kalunit)corrmatrixtemp_triangle(i1)
                  else
                    write(kalunit,'(e20.12)')corrmatrixtemp_triangle(i1)
                  endif
                enddo
                deallocate(corrmatrixtemp_triangle)
                if(lwriteunformatted)then
                  write(kalunit)iseed
                  write(kalunit)kseed
                  write(kalunit)lseed
                  write(kalunit)mseed
                else
                  write(kalunit,*)iseed
                  write(kalunit,*)kseed
                  write(kalunit,*)lseed
                  write(kalunit,*)mseed
                endif
              endif ! mpirank.eq.0
            endif ! mpisize
            call mpi_barrier(mpi_comm_world,mpierror)
            if(mpirank.eq.0)then
              close(kalunit)
            endif ! mpirank.eq.0
          enddo ! i4=1,ndim
!!
      endif ! lsavekalman
!!
      return
      end
