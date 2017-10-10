!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine readotf(countepoch)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer countepoch                                ! in

      character*40 keyword                              ! internal





!===========================================
! read input.otf file (only process 0)
!===========================================
      if(mpirank.eq.0)then
        open(otfunit,file='input.otf',form='formatted',status='old')
 50     continue
        read(otfunit,*,END=51) keyword 

        if(keyword.eq.'write_trainforces')then
          if(.not.lwritetrainforces)then
            write(ounit,*)countepoch,&
              ' OTF switching on write_trainforces'
            lwritetrainforces=.true.
          endif

        elseif(keyword.eq.'write_trainpoints')then
          if(.not.lwritetrainpoints)then
            write(ounit,*)countepoch,&
              ' OTF switching on write_trainpoints'
            lwritetrainpoints=.true.
          endif

        else
          write(ounit,*)'ERROR: invalid keyword in input.otf ',keyword
          stop

        endif

        goto 50
 51     continue
        close(otfunit)

      endif ! mpirank.eq.0

!===========================================
! distribute modified parameters to all other processes (MPI) 
!===========================================
      if(mpisize.gt.1)then
        call mpi_bcast(lwritetrainforces,1,mpi_logical,0,mpi_comm_world,mpierror)
        call mpi_bcast(lwritetrainpoints,1,mpi_logical,0,mpi_comm_world,mpierror)
      endif

      return
      end
