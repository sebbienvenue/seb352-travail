!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine getkalmanmatricespair(&
         corrdim,corrfdim,corredim,corrcdim,&
         maxcorrdim,maxcorrfdim,maxcorredim,&
         num_weightsewaldfree,num_weightspairfree,&
         corrmatrix_list,corrmatrixf_list,corrmatrixe_list,corrmatrixc,&
         lshort,lewald)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3                ! internal
      integer num_weightsewaldfree(nelem)! in
      integer num_weightspairfree(npairs)
      integer itemp                      ! internal
      integer n_start                    ! internal
      integer n_end                      ! internal
!! Kalman matrix dimensions:
      integer corrdim(npairs)            ! in
      integer corrfdim(npairs)           ! in
      integer corredim(nelem)            ! in
      integer corrcdim                   ! in
      integer maxcorrdim                 ! in
      integer maxcorrfdim                ! in
      integer maxcorredim                ! in
      integer isum                       ! internal
!!
      real*8 corrmatrix_list(maxcorrdim,npairs)                      ! out
      real*8 corrmatrixf_list(maxcorrfdim,npairs)                    ! out
      real*8 corrmatrixe_list(maxcorredim,nelem)                     ! out
      real*8 corrmatrixc(corrcdim)                                   ! out
      real*8, dimension(:)  , allocatable :: corrmatrix_full         ! internal 
!!
      character*21 kalfilenames(nelem)  ! internal
      character*21 ctemp21              ! internal
      character*20 kalfilenamee(nelem)  ! internal
      character*20 ctemp20              ! internal
      character*24 kalfilenamep(npairs)  ! internal
      character*24 ctemp24              ! internal
!!
      logical lshort                    ! in
      logical lewald                    ! in
!!
!! initializations
      isum=0
      kalfilenames(:)='kalman.short.000.data'
      kalfilenamep(:)='kalman.pair.000.000.data'
      kalfilenamee(:)='kalman.elec.000.data'

      do i1=1,nelem
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
      if(mpisize.eq.1)then  !------------------------IF
        if(lewald.and.lchargeconstraint)then
          do i1=1,nelem
            isum=isum+num_weightsewaldfree(i1)
          enddo
        endif
        do i1=1,npairs 
          if(lshort.and.(optmodee.eq.1))then
            call initialcorrmatrix(num_weightspairfree(i1),corrmatrix_list(1,i1))
          endif
          if(lshort.and.(optmodef.eq.1))then
            call initialcorrmatrix(num_weightspairfree(i1),corrmatrixf_list(1,i1))
          endif
        enddo

        do i1=1,nelem
          if(lewald.and.(optmodeq.eq.1))then
            call initialcorrmatrix(num_weightsewaldfree(i1),corrmatrixe_list(1,i1))
            if(lchargeconstraint)then
              call initialcorrmatrix(isum,corrmatrixc)
            endif
          endif
        enddo
!!
      else ! parallel case   -------------------ELSE
!! short range part for energy fitting
        if(lshort.and.(optmodee.eq.1))then

          do i1=1,npairs 
            allocate(corrmatrix_full(num_weightspairfree(i1)*num_weightspairfree(i1)))
            corrmatrix_full(:)=0.0d0
!! set diagonal elements to 1.0d0
            do i2=1,num_weightspairfree(i1)
              do i3=1,num_weightspairfree(i1)
                if(i2.eq.i3)then
                  itemp=(i2-1)*num_weightspairfree(i1)+i3
                  corrmatrix_full(itemp)=1.0d0
                endif
              enddo ! i3
            enddo ! i2
!! distribute correlation matrix to processes, distribute only full columns
            call mpifitdistribution(num_weightspairfree(i1),itemp,n_start,n_end)
            do i2=1,corrdim(i1)
              corrmatrix_list(i2,i1)=corrmatrix_full(i2+(n_start-1)*num_weightspairfree(i1))
            enddo ! i2
            deallocate(corrmatrix_full)
          enddo ! i1
        endif ! lshort
!!
!! short range part for force fitting
        if(lshort.and.(optmodef.eq.1))then
          do i1=1,npairs
            allocate(corrmatrix_full(num_weightspairfree(i1)*num_weightspairfree(i1)))
            corrmatrix_full(:)=0.0d0
!! set diagonal elements to 1.0d0
            do i2=1,num_weightspairfree(i1)
              do i3=1,num_weightspairfree(i1)
                if(i2.eq.i3)then
                  itemp=(i2-1)*num_weightspairfree(i1)+i3
                  corrmatrix_full(itemp)=1.0d0
                endif
              enddo ! i3
            enddo ! i2
!! distribute correlation matrix to processes, distribute only full columns
            call mpifitdistribution(num_weightspairfree(i1),itemp,n_start,n_end)
            do i2=1,corrfdim(i1)
              corrmatrixf_list(i2,i1)&
                =corrmatrix_full(i2+(n_start-1)*num_weightspairfree(i1))
            enddo ! i2
            deallocate(corrmatrix_full)
          enddo ! i1
        endif ! lshort
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 lshort

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
      endif ! mpisize --------------------------------ENDIF
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read Kalman filter data if requested
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(lrestkalman)then                             !-------FOR RESTARTING THE CALCULATIONS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read Kalman data for short range NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(ounit,*)'Error: no restart for Kalman matrices implemented for nntype 2'
        stop
!!
      endif ! lrestkalman
!!
!!
      return
      end
