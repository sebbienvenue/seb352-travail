!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! called by:
!!
      subroutine writekalman_elec(&
        iseed,kseed,lseed,mseed,&
        maxcorredim,kaledim,corrcdim,&
        num_weightsewaldfree,&
        corrmatrixe_list,corrmatrixc)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer num_weightsewaldfree(nelem)  ! in
      integer i1,i2,i4          ! internal
      integer icount            ! internal
      integer jcount            ! internal
      integer, dimension(:)  , allocatable :: kaldim_list 
      integer, dimension(:)  , allocatable :: displs 
      integer, dimension(:)  , allocatable :: recvcnts 
!! Kalman matrix dimensions:
      integer maxcorredim          ! in
      integer corrcdim             ! in
      integer kaledim(nelem)       ! in
      integer constraintdim        ! internal
      integer iseed
      integer kseed
      integer lseed
      integer mseed
!!
      real*8 corrmatrixe_list(maxcorredim,nelem) ! in
      real*8 corrmatrixc(corrcdim)               ! in
      real*8, dimension(:)  , allocatable :: corrmatrixtemp 
      real*8, dimension(:)  , allocatable :: corrmatrixtemp_triangle 
!!
      character*20 kalfilenamee(nelem)  ! internal
      character*20 ctemp20              ! internal
!!
!! initializations
      kalfilenamee(:)='kalman.elec.000.data'
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
      if(lsavekalman)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write Kalman data for electrostatic NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i4=1,nelem
          if(mpirank.eq.0)then
            open(kaleunit,file=kalfilenamee(i4),form='formatted',status='replace')
            if(lwriteunformatted)then
              write(kaleunit)kalmanlambdae(i4) !'
            else
              write(kaleunit,*)kalmanlambdae(i4) !'
            endif
          endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write corrmatrixe for charge weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for mpisize.eq.1 we do not need to compress
          if((mpisize.eq.1).or.lompmkl)then
            if(mpirank.eq.0)then
              do i1=1,num_weightsewaldfree(i4)*(num_weightsewaldfree(i4)+1)/2
                if(lwriteunformatted)then
                  write(kaleunit)corrmatrixe_list(i1,i4)
                else
                  write(kaleunit,*)corrmatrixe_list(i1,i4)
                endif
              enddo
              if(lwriteunformatted)then
                write(kaleunit)iseed
                write(kaleunit)kseed
                write(kaleunit)lseed
                write(kaleunit)mseed
              else
                write(kaleunit,*)iseed
                write(kaleunit,*)kseed
                write(kaleunit,*)lseed
                write(kaleunit,*)mseed
              endif
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
                if(lwriteunformatted)then
                  write(kaleunit)corrmatrixtemp_triangle(i1)
                else
                  write(kaleunit,*)corrmatrixtemp_triangle(i1)
                endif
              enddo
              deallocate(corrmatrixtemp_triangle)
              if(lwriteunformatted)then
                write(kaleunit)iseed
                write(kaleunit)kseed
                write(kaleunit)lseed
                write(kaleunit)mseed
              else
                write(kaleunit,*)iseed
                write(kaleunit,*)kseed
                write(kaleunit,*)lseed
                write(kaleunit,*)mseed
              endif
            endif ! mpirank.eq.0
          endif ! mpisize
          call mpi_barrier(mpi_comm_world,mpierror)
          if(mpirank.eq.0)then
            close(kaleunit)
          endif ! mpirank.eq.0
        enddo ! i4=1,nelem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write charge constraint data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(lchargeconstraint)then
          constraintdim=0
          do i4=1,nelem
            constraintdim=constraintdim+num_weightsewaldfree(i4)
          enddo ! i4
          if(mpirank.eq.0)then
            open(kalcunit,file='kalmanc.data',form='formatted',status='replace')
              if(lwriteunformatted)then
                write(kalcunit)kalmanlambdac !'
              else
                write(kalcunit,*)kalmanlambdac !'
              endif
              do i1=1,constraintdim*(constraintdim+1)/2
                if(lwriteunformatted)then
                  write(kalcunit)corrmatrixc(i1)
                else
                  write(kalcunit,*)corrmatrixc(i1)
                endif
              enddo
            close(kalcunit)
          endif ! mpirank.eq.0
        endif ! lchargeconstraint
      endif ! lsavekalman
!!
      return
      end
