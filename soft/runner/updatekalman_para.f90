!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - optimize_short_combined.f90
!! - optimize_short_combinedpair.f90
!! - optimize_ewald.f90
!! - fitforcesshort.f90
!!
      subroutine updatekalman_para(paramode,kaldim,corrdim,&
                num_weights,kalmanlambda,kalmannue,&
                weights,dedw,corrmatrix,error)
!!
      use mpi_mod
      use fileunits
!! Don't use fittingoptions, because updatekalman is used for different quantities
!!
      implicit none
!!
      integer num_weights                                         ! in
      integer i1,i2,i3,i4                                         ! internal
      integer nstruct                                             ! internal
      integer n_start                                             ! internal
      integer n_end                                               ! internal
      integer kaldim                                              ! in
      integer corrdim                                             ! in
      integer icount                                              ! internal
      integer jcount                                              ! internal
      integer paramode                                            ! in
!! 
      real*8 kalmanlambda                                         ! in/out
      real*8 kalmannue                                            ! in
      real*8 weights(num_weights)                                 ! in/out
      real*8 wtemp(num_weights)                                   ! internal
      real*8 dedw(num_weights)                                    ! in
      real*8 kalgainmat(kaldim)                                   ! internal
      real*8 corrmatrix(corrdim)                                  ! in/out
      real*8 error                                                ! in
      real*8 alpha                                                ! internal
      real*8 inverse                                              ! internal
      real*8 invlambda                                            ! internal
      real*8 coh(num_weights)                                     ! internal
      real*8 ddot                                                 ! internal
      real*8 vec(num_weights)                                     ! internal
      real*8 checksum
!!
!! initializations
      coh(:)   = 0.0d0
      invlambda= 1.d0/kalmanlambda
      inverse  = 0.0d0
!!
!! debug
!!      checksum=0.0d0
!!      do i1=1,corrdim
!!        checksum=checksum+corrmatrix(i1)
!!      enddo
!!      call mpi_allreduce(mpi_in_place,checksum,1,&
!!       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!      write(ounit,*)mpirank,' corrmatrix ',checksum
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!
!! Step1) calculation of the vector coh: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpifitdistribution(num_weights,nstruct,n_start,n_end)
!!      write(ounit,*)mpirank,'nstruct,n_start,n_end ',nstruct,n_start,n_end
!!
!! calculate part of the coh of this process 
!!      icount=0
!!      do i1=n_start,n_end
!!        do i2=1,num_weights
!!          vec(i2)=corrmatrix(icount+i2)
!!!!          write(ounit,*)mpirank,i2,vec(i2)
!!        enddo ! i2
!!        icount=icount+num_weights
!!        coh(i1)=ddot(num_weights,dedw,1,vec,1)
!!      enddo ! i1

!! alternative (both seem to work correctly)
      icount=0
      do i1=n_start,n_end
        coh(i1)=ddot(num_weights,dedw,1,corrmatrix(icount*num_weights+1),1)
        icount=icount+1
      enddo ! i1
!!
!! distribute coh matrix, needed in step 5 
      call mpi_allreduce(mpi_in_place,coh,num_weights,&
       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! debug
!!      checksum=0.0d0
!!      do i1=1,num_weights
!!        checksum=checksum+coh(i1)
!!      enddo
!!      write(ounit,*)mpirank,' coh ',checksum
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! debug
!!      if(mpirank.eq.0)then
!!        do i1=1,num_weights
!!         write(ounit,*)'coh ',i1,coh(i1)
!!        enddo
!!      endif
!!      stop
!!
!! Step2) calculation of inverse:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! parallel version
!! CAUTION: In the parallel version the result very slightly
!! depends on the number of processes because the order of the
!! summation matters. Even invisible changes after the 15th digit
!! still sum up to significant errors in long fits
!! => don't use the parallel version!!!
!!
      if(paramode.eq.1)then
!! serial version 
        inverse=ddot(num_weights,dedw,1,coh,1)
      elseif(paramode.eq.2)then
!! calculate part of the scalar product of this process 
!!        do i1=n_start,n_end
!!          inverse=inverse+dedw(i1)*coh(i1)  
!!        enddo ! i1
!! alternative (both seem to work)
        inverse=ddot(nstruct,dedw(n_start),1,coh(n_start),1)
!! combine the partial sums of all processes
        call mpi_allreduce(mpi_in_place,inverse,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      else
        write(ounit,*)'Error: unsupported paramode in updatekalman_para'
        stop
      endif

!! final step
      inverse = 1.d0/(kalmanlambda+inverse)
!!
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!
!! Step3) update of kalgainmat:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the part of kalgainmat of this process
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        kalgainmat(icount)=inverse*coh(i1)
!!        write(ounit,*)'kalgainmat ',i1,inverse,coh(i1),kalgainmat(icount)
      enddo ! i1

!! debug
!!      checksum=0.0d0
!!      do i1=1,nstruct
!!        checksum=checksum+kalgainmat(i1)
!!      enddo
!!      call mpi_allreduce(mpi_in_place,checksum,1,&
!!       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!      write(ounit,*)mpirank,' kalgainmat ',checksum
!!      call mpi_barrier(mpi_comm_world,mpierror)
!! debug
!!      if(mpirank.eq.0)then
!!      do i1=1,kaldim
!!        write(ounit,*)mpirank,' kalgainmat ',i1,kalgainmat(i1) 
!!      enddo
!!      endif
!!
!! Step4) update weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      icount=0
      wtemp(:)=0.0d0
      do i1=n_start,n_end
        icount=icount+1
        wtemp(i1)=weights(i1)+error*kalgainmat(icount)
!!        write(ounit,*)'updatekalman ',i1,error,kalgainmat(icount)
      enddo ! i1
!!
!! combine the wtemp of all processes
      call mpi_allreduce(wtemp,weights,num_weights,&
       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!!
!! debug
!!      checksum=0.0d0
!!      do i1=1,num_weights
!!        checksum=checksum+weights(i1)
!!      enddo
!!      write(ounit,*)mpirank,' weights ',checksum
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! Step5) update correlation matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      alpha=-inverse*invlambda
      icount=0
      do i1=n_start,n_end
        do i2=1,num_weights
          icount=icount+1
          corrmatrix(icount)=invlambda*corrmatrix(icount)+alpha*coh(i1)*coh(i2)
        enddo ! i2
      enddo ! i1

!! debug
!!      checksum=0.0d0
!!      do i1=1,corrdim
!!        checksum=checksum+corrmatrix(i1)
!!      enddo
!!      call mpi_allreduce(mpi_in_place,checksum,1,&
!!       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!      write(ounit,*)mpirank,' corrmatrix ',checksum
!!      call mpi_barrier(mpi_comm_world,mpierror)

!!      icount=n_start
!!      jcount=1
!!      do i1=1,corrdim
!!        corrmatrix(i1)=corrmatrix(i1)+alpha*coh(icount)*coh(jcount)
!!        icount=icount+1
!!        if(icount.gt.(kaldim+n_start)) icount=n_start
!!        jcount=jcount+1
!!        if(jcount.gt.num_weights) jcount=1
!!      enddo
!!
!! This slows down the code a lot here and thus is done already above:
!!      corrmatrix(:)=invlambda*corrmatrix(:)

!! debug
!!      if(mpirank.eq.0)then
!!        do i1=1,corrdim
!!         write(ounit,*)'corrmat ',i1,corrmatrix(i1)
!!        enddo
!!      endif
!!      stop
!!
!!
!! Step6) update kalmanlambda 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      kalmanlambda=kalmannue*kalmanlambda + 1.d0 - kalmannue
!!
!!
!!
!!
      return
      end
