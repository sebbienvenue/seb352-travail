!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine writekalmanunformattedpair(&
        iseed,kseed,lseed,mseed,&
        maxcorrdim,maxcorrfdim,maxcorredim,&
        kaldim,kaledim,corrcdim,&
        num_weightspairfree,num_weightsewaldfree,&
        corrmatrix_list,corrmatrixe_list,corrmatrixc,&
        lshort,lewald)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer num_weightspairfree(npairs)  ! in
      integer num_weightsewaldfree(nelem)  ! in
      integer i1,i2,i4       ! internal
      integer icount            ! internal
      integer jcount            ! internal
      integer, dimension(:)  , allocatable :: kaldim_list 
      integer, dimension(:)  , allocatable :: displs 
      integer, dimension(:)  , allocatable :: recvcnts 
!! Kalman matrix dimensions:
      integer maxcorrdim           ! in
      integer maxcorrfdim          ! in
      integer maxcorredim          ! in
      integer corrcdim             ! in
      integer kaldim(npairs)         ! in
      integer kaledim(nelem)       ! in
      integer constraintdim        ! internal
      integer iseed
      integer kseed
      integer lseed
      integer mseed

      real*8 corrmatrix_list(maxcorrdim,npairs)   ! in
      real*8 corrmatrixe_list(maxcorredim,nelem) ! in
      real*8 corrmatrixc(corrcdim)               ! in
      real*8, dimension(:)  , allocatable :: corrmatrixtemp 
      real*8, dimension(:)  , allocatable :: corrmatrixtemp_triangle 
!!
!!
      character*21 kalfilenamesp(npairs)  ! internal
      character*21 ctemp21              ! internal
      character*20 kalfilenamee(nelem)  ! internal
      character*20 ctemp20              ! internal
!!
      logical lshort            ! in
      logical lewald            ! in
!!
      write(ounit,*)'WARNING: correlation matrix for forces is not written to file'
!!'
!! initializations
      kalfilenamesp(:)='kalman.short.000.data'
      kalfilenamee(:)='kalman.elec.000.data'
      do i1=1,nelem
        ctemp21=kalfilenamesp(i1)
        if(nucelem(i1).gt.99)then
          write(ctemp21(14:16),'(i3)')nucelem(i1)
        elseif(nucelem(i1).gt.9)then
          write(ctemp21(15:16),'(i2)')nucelem(i1)
        else
          write(ctemp21(16:16),'(i1)')nucelem(i1)
        endif
        kalfilenamesp(i1)=ctemp21
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
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!      write(ounit,*)mpirank,' mpisize ',mpisize
!!
!!      write(ounit,*)mpirank,' writekalman starts ',mpisize
!!
      if(lsavekalman)then
!! write Kalman data for short range NN
        if(lshort)then
          do i4=1,nelem
            if(mpirank.eq.0)then
              open(kalunit,file=kalfilenamesp(i4),form='unformatted',status='replace') 
              write(kalunit)kalmanlambdap(i4) !'
            endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write corrmatrix for short range weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for mpisize.eq.1 we do not need to compress
            if(mpisize.eq.1)then
              if(mpirank.eq.0)then
                do i1=1,num_weightspairfree(i4)*(num_weightspairfree(i4)+1)/2
                  write(kalunit)corrmatrix_list(i1,i4)
                enddo
              endif ! mpirank.eq.0
              write(kalunit)iseed
              write(kalunit)kseed
              write(kalunit)lseed
              write(kalunit)mseed
            else ! more than one core
!! collect the correlation matrices of all processes
              if(mpirank.eq.0)then
                allocate(corrmatrixtemp_triangle(num_weightspairfree(i4)*(num_weightspairfree(i4)+1)/2))
              endif ! mpirank.eq.0
              allocate(corrmatrixtemp(num_weightspairfree(i4)*num_weightspairfree(i4)))
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
              kaldim_list(:)=kaldim_list(:)*num_weightspairfree(i4)
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
                do i1=1,num_weightspairfree(i4)
                  do i2=1,num_weightspairfree(i4)
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
                do i1=1,num_weightspairfree(i4)*(num_weightspairfree(i4)+1)/2
                  write(kalunit)corrmatrixtemp_triangle(i1)
                enddo
                deallocate(corrmatrixtemp_triangle)
                write(kalunit)iseed
                write(kalunit)kseed
                write(kalunit)lseed
                write(kalunit)mseed
              endif ! mpirank.eq.0
            endif ! mpisize
            call mpi_barrier(mpi_comm_world,mpierror)
            if(mpirank.eq.0)then
              close(kalunit)
            endif ! mpirank.eq.0
          enddo ! i4=1,nelem
        endif ! lshort
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write Kalman data for electrostatic NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(lewald)then
          do i4=1,nelem
            if(mpirank.eq.0)then
              open(kaleunit,file=kalfilenamee(i4),form='unformatted',status='replace')
              write(kaleunit)kalmanlambdae(i4) !'
            endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write corrmatrixe for charge weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for mpisize.eq.1 we do not need to compress
            if(mpisize.eq.1)then
              if(mpirank.eq.0)then
                do i1=1,num_weightsewaldfree(i4)*(num_weightsewaldfree(i4)+1)/2
                  write(kaleunit)corrmatrixe_list(i1,i4)
                enddo
                write(kaleunit)iseed
                write(kaleunit)kseed
                write(kaleunit)lseed
                write(kaleunit)mseed
              endif ! mpirank.eq.0
            else ! more than one core
!! collect the correlation matrices of all processes
              if(mpirank.eq.0)then
                allocate(corrmatrixtemp_triangle(num_weightsewaldfree(i4)*(num_weightsewaldfree(i4)+1)/2))
              endif ! mpirank.eq.0
              allocate(corrmatrixtemp(num_weightsewaldfree(i4)*num_weightsewaldfree(i4)))
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
!! determine the corrmatrixe fragment size of all processes
              call mpi_gatherv(kaledim(i4),1,mpi_integer,&
                kaldim_list,recvcnts,displs,mpi_integer,0,mpi_comm_world,mpierror)
!!
              kaldim_list(:)=kaldim_list(:)*num_weightsewaldfree(i4)
              call mpi_bcast(kaldim_list,mpisize,mpi_integer,0,mpi_comm_world,mpierror)
!! determine the displacements for corrmatrixtemp 
              displs(:)=0
              do i1=2,mpisize
                do i2=1,i1-1
                  displs(i1)=displs(i1)+kaldim_list(i2)
                enddo ! i2
              enddo ! i1
!! collect the corrmatrixe fragments of all processes
              call mpi_gatherv(corrmatrixe_list(1,i4),&
                kaldim_list(mpirank+1),mpi_real8,&
                corrmatrixtemp,kaldim_list,displs,mpi_real8,0,mpi_comm_world,mpierror)
              call mpi_barrier(mpi_comm_world,mpierror)
              deallocate(kaldim_list)
              deallocate(displs)
              deallocate(recvcnts)
!! compress corrmatrixe to lower triangle before saving'
              if(mpirank.eq.0)then
                icount=0
                jcount=0
                do i1=1,num_weightsewaldfree(i4)
                  do i2=1,num_weightsewaldfree(i4)
                    jcount=jcount+1
                    if(i2.ge.i1)then
                      icount=icount+1
                      corrmatrixtemp_triangle(icount)=corrmatrixtemp(jcount)
                    endif
                  enddo ! i2
                enddo ! i1
              endif ! mpirank.eq.0
              deallocate(corrmatrixtemp)
!! now save correlation matrix for charge weights
              if(mpirank.eq.0)then
                do i1=1,num_weightsewaldfree(i4)*(num_weightsewaldfree(i4)+1)/2
                  write(kaleunit)corrmatrixtemp_triangle(i1)
                enddo
                deallocate(corrmatrixtemp_triangle)
                write(kaleunit)iseed
                write(kaleunit)kseed
                write(kaleunit)lseed
                write(kaleunit)mseed
              endif ! mpirank.eq.0
            endif ! mpisize
            call mpi_barrier(mpi_comm_world,mpierror)
            if(mpirank.eq.0)then
              close(kaleunit)
            endif ! mpirank.eq.0
          enddo ! i4=1,nelem
        endif ! lewald


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write charge constraint data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(lewald)then
          if(lchargeconstraint)then
            constraintdim=0
            do i4=1,nelem
              constraintdim=constraintdim+num_weightsewaldfree(i4)
            enddo ! i4
            if(mpirank.eq.0)then
              open(kalcunit,file='kalmanc.data',form='unformatted',status='replace')
                write(kalcunit)kalmanlambdac !'
                do i1=1,constraintdim*(constraintdim+1)/2
                  write(kalcunit)corrmatrixc(i1)
                enddo
              close(kalcunit)
            endif ! mpirank.eq.0
          endif ! lchargeconstraint
        endif ! lewald
      endif ! lsavekalman
!!
      return
      end
