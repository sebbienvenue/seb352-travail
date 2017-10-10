!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting_batch.f90
!!
!! This routine is called once for each block of points
!!
      subroutine updateshort(nshort,&
         maxkaldim,kaldim,corrdim,maxcorrdim,&
         num_weights_short_atomic_free,&
         wconstraintidx,&
         dshortdw,errorshort,&
         corrmatrix_list)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer num_weights_short_atomic_free(nelem)              ! in
      integer wconstraintidx(maxnum_weights_short_atomic,nelem) ! in
      integer i1,i2,i3                                 ! internal
!! Kalman matrix dimensions:
      integer maxcorrdim                               ! in
      integer maxkaldim                                ! in
      integer corrdim(nelem)                           ! in
      integer kaldim(nelem)                            ! in
      integer nshort(nelem)                            ! in
!!
      real*8 deshortdw(maxnum_weights_short_atomic,1,nelem)       ! in
      real*8 dshortdw(maxnum_weights_short_atomic,nelem)       ! in
      real*8 errorshort(nelem)                         ! in
!!
      real*8 errore                                    ! internal
      real*8, dimension(:)  , allocatable :: weights                        ! internal
      real*8, dimension(:) , allocatable :: deshortdw_temp                  ! internal
      real*8 corrmatrix_list(maxcorrdim,nelem)                              ! in/out
!!
!!
!! initializations
      errore                 = 0.0d0

!      write(ounit,*)'ounit ',ounit
!      write(ounit,*)'nshort ',nshort(1)
!      write(ounit,*)'maxkaldim ',maxkaldim
!      write(ounit,*)'kaldim ',kaldim(1)
!      write(ounit,*)'maxcorrdim ',maxcorrdim
!      write(ounit,*)'nelem ',nelem
!      write(ounit,*)'maxnum_weightsshort ',maxnum_weightsshort
!      write(ounit,*)'num_weightsshort ',num_weightsshort(1)
!      write(ounit,*)'num_weightsshortfree ',num_weightsshort(1)
!      write(ounit,*)'errorshort ',errorshort(1)
!      write(ounit,*)'kalmanlambda ',kalmanlambda(1)
!      write(ounit,*)'kalmannue ',kalmannue
!      write(ounit,*)'corrdim inside updateshort ',corrdim(1)
!!      errorshort(:)=0.01d0
!!      write(ounit,*)'setting errorshort to 0.01'
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update the weights for each element:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! normalize the error and the derivatives in case grouping has been used
!!
!! chose optimization algorithm
      if(optmodee.eq.1)then ! Kalman filter
!! update weights for each element
        do i2=1,nelem
!!
          errore =errorshort(i2)/dble(nshort(i2))
          deshortdw(:,1,i2)=dshortdw(:,i2)/dble(nshort(i2))
!!
!! prepare temporary arrays
          allocate(weights(num_weights_short_atomic_free(i2)))
          allocate(deshortdw_temp(num_weights_short_atomic_free(i2)))
!!
!! reduce weights array to free weights
          do i3=1,num_weights_short_atomic_free(i2)
            weights(i3)=weights_short_atomic(wconstraintidx(i3,i2),i2)
            deshortdw_temp(i3)=deshortdw(i3,1,i2)
!!            write(ounit,*)'deshortdw_temp ',deshortdw_temp(i3)
          enddo ! i3
!!
          if(mpisize.eq.1)then
!! updatekalman_para cannot be used for 1 process because the corrmatrix array has a different dimension then
            if(mpirank.eq.0)then
              call updatekalman(kaldim(i2),corrdim(i2),&
                num_weights_short_atomic_free(i2),kalmanlambda(i2),kalmannue,&
                weights,deshortdw_temp,&
                corrmatrix_list(1,i2),errore)
            endif ! mpirank.eq.0
            call mpi_bcast(weights,num_weights_short_atomic_free(i2),&
              mpi_real8,0,mpi_comm_world,mpierror)
!!
          else ! real parallel case
!! updatekalman cannot be used for more than 1 process because the corrmatrix array has a different dimension then
            call updatekalman_para(paramode,kaldim(i2),corrdim(i2),&
              num_weights_short_atomic_free(i2),kalmanlambda(i2),kalmannue,&
              weights,deshortdw_temp,&
              corrmatrix_list(1,i2),errore)
          endif ! mpisize.eq.1
!!
!! expand weights array back to original array
          do i3=1,num_weights_short_atomic_free(i2)
            weights_short_atomic(wconstraintidx(i3,i2),i2)=weights(i3)
          enddo ! i3
!!
          deallocate(weights)
          deallocate(deshortdw_temp)
        enddo ! i2=1,nelem
!!
      elseif(optmodee.eq.2)then ! conjugate gradient'
        write(ounit,*)'CG optimization not yet implemented in optimize_short'
        stop
!!
!! '
      elseif(optmodee.eq.3)then ! steepest descent'
!! update weights for each element
        do i2=1,nelem
          errore =errorshort(i2)/dble(nshort(i2))
          deshortdw(:,1,i2)=dshortdw(:,i2)/dble(nshort(i2))
!!
          allocate(weights(num_weights_short_atomic_free(i2)))
          allocate(deshortdw_temp(num_weights_short_atomic_free(i2)))
!!
         if(mpirank.eq.0)then
!!
!! reduce weights array to free weights
           do i3=1,num_weights_short_atomic_free(i2)
             weights(i3)=weights_short_atomic(wconstraintidx(i3,i2),i2)
             deshortdw_temp(i3)=deshortdw(i3,1,i2)
           enddo ! i3
!!
           call updatesteepest(&
             num_weights_short_atomic_free(i2),weights,&
             deshortdw_temp,errore,steepeststepe)
!!
!! expand weights array back to original array
             do i3=1,num_weights_short_atomic_free(i2)
               weights_short_atomic(wconstraintidx(i3,i2),i2)=weights(i3)
             enddo ! i3
!!
           endif ! mpirank.eq.0
           call mpi_bcast(weights_short_atomic,maxnum_weights_short_atomic*nelem,&
             mpi_real8,0,mpi_comm_world,mpierror)
!!
           deallocate(weights)
           deallocate(deshortdw_temp)
         enddo ! i2=1,nelem
       endif ! optmodee
!!
      return
      end
