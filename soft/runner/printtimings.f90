!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine printtimings()
!!
      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions
      use timings
!!
      implicit none
!!
      integer i1
!!
      real*8, dimension(:), allocatable :: time_array_local 
!! mode 3
      real*8 maxtimeshort
      real*8 maxtimeallocshort
      real*8 maxtimesymshort
      real*8 maxtimeextrapolationshort
      real*8 maxtimescalesymshort
      real*8 maxtimescaledsfuncshort
      real*8 maxtimeeshort
      real*8 maxtimefshort
      real*8 maxtimesshort
      real*8 maxtimesymelec1
      real*8 maxtimesymelec2
      real*8 maxtimecharge
      real*8 maxtimecomm1
      real*8 maxtimeelec
      real*8 maxtimeeelec
!!
      allocate(time_array_local(mpisize))
!!
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Timing summary:                   t/min or  t/s'
        write(ounit,*)'-------------------------------------------------------------'
      endif ! mpirank.eq.0



      if(mode.eq.3)then
!! get maximum process times for mode 3
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timeshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimeshort=0.0d0
          do i1=1,mpisize
            maxtimeshort=max(maxtimeshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timeallocshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimeallocshort=0.0d0
          do i1=1,mpisize
            maxtimeallocshort=max(maxtimeallocshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timesymshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimesymshort=0.0d0
          do i1=1,mpisize
            maxtimesymshort=max(maxtimesymshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timeextrapolationshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimeextrapolationshort=0.0d0
          do i1=1,mpisize
            maxtimeextrapolationshort=max(maxtimeextrapolationshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timescalesymshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimescalesymshort=0.0d0
          do i1=1,mpisize
            maxtimescalesymshort=max(maxtimescalesymshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timescaledsfuncshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimescaledsfuncshort=0.0d0
          do i1=1,mpisize
            maxtimescaledsfuncshort=max(maxtimescaledsfuncshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timeeshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimeeshort=0.0d0
          do i1=1,mpisize
            maxtimeeshort=max(maxtimeeshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timefshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimefshort=0.0d0
          do i1=1,mpisize
            maxtimefshort=max(maxtimefshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timesshort
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimesshort=0.0d0
          do i1=1,mpisize
            maxtimesshort=max(maxtimesshort,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timeelec
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimeelec=0.0d0
          do i1=1,mpisize
            maxtimeelec=max(maxtimeelec,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timesymelec1
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimesymelec1=0.0d0
          do i1=1,mpisize
            maxtimesymelec1=max(maxtimesymelec1,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timecharge
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimecharge=0.0d0
          do i1=1,mpisize
            maxtimecharge=max(maxtimecharge,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timecomm1
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimecomm1=0.0d0
          do i1=1,mpisize
            maxtimecomm1=max(maxtimecomm1,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timesymelec2
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimesymelec2=0.0d0
          do i1=1,mpisize
            maxtimesymelec2=max(maxtimesymelec2,time_array_local(i1))
          enddo
        endif
!!
        time_array_local(:)=0.0d0
        time_array_local(mpirank+1)=timeeelec
        call mpi_allreduce(mpi_in_place,time_array_local,mpisize,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(mpirank.eq.0)then
          maxtimeeelec=0.0d0
          do i1=1,mpisize
            maxtimeeelec=max(maxtimeeelec,time_array_local(i1))
          enddo
        endif
!!
!! summation for mode 3
        call mpi_allreduce(mpi_in_place,timeshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timeallocshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timesymshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timeextrapolationshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timescalesymshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timescaledsfuncshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timeeshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timefshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timesshort,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timesymelec1,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timecharge,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timecomm1,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timeelec,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timesymelec2,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,timeeelec,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      endif ! mode.eq.3

      if(mpirank.eq.0)then
!!
        timeinitnn=timeinitnnend-timeinitnnstart
        write(ounit,'(a35,x,2f8.2)')' TOTTIME INITNN                 ',&
          timeinitnn/60.d0,timeinitnn
        timereadinput=timereadinputend-timereadinputstart
        write(ounit,'(a35,x,2f8.2)')'   TOTTIME READINPUT            ',&
          timereadinput/60.d0,timereadinput
!!
        if(mode.eq.1)then
          timemode1=timemode1end-timemode1start
          write(ounit,'(a35,x,2f8.2)')' TOTTIME MODE 1                 ',&
            timemode1/60.d0,timemode1
        endif
!!
        if(mode.eq.2)then
          timemode2=timemode2end-timemode2start
          write(ounit,'(a35,x,2f8.2)')' TOTTIME MODE 2                 ',&
            timemode2/60.d0,timemode2
        endif
!!
        if(mode.eq.3)then
          timemode3=timemode3end-timemode3start
          write(ounit,'(a35,x,2f8.2)')' TOTTIME MODE 3                   ',&
            timemode3/60.d0,timemode3
          write(ounit,'(a35,x,2f8.2)')'   TOTTIME SHORT RANGE PART       ',&
            timeshort/60.d0,timeshort
          write(ounit,'(a35,x,2f8.2)')'   average SHORT RANGE PART       ',&
            timeshort/(60.d0*dble(mpisize)),timeshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'   maximum SHORT RANGE PART       ',&
            maxtimeshort/(60.d0),maxtimeshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME ALLOC SHORT          ',&
            timeallocshort/60.d0,timeallocshort
          write(ounit,'(a35,x,2f8.2)')'     average ALLOC SHORT          ',&
            timeallocshort/(60.d0*dble(mpisize)),timeallocshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum ALLOC SHORT          ',&
            maxtimeallocshort/(60.d0),maxtimeallocshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME SYMFUNCTION SHORT    ',&
            timesymshort/60.d0,timesymshort
          write(ounit,'(a35,x,2f8.2)')'     average SYMFUNCTION SHORT    ',&
            timesymshort/(60.d0*dble(mpisize)),timesymshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum SYMFUNCTION SHORT    ',&
            maxtimesymshort/(60.d0),maxtimesymshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME EXTRAPOLATION SHORT  ',&
            timeextrapolationshort/60.d0,timeextrapolationshort
          write(ounit,'(a35,x,2f8.2)')'     average EXTRAPOLATION SHORT  ',&
            timeextrapolationshort/(60.d0*dble(mpisize)),timeextrapolationshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum EXTRAPOLATION SHORT  ',&
            maxtimeextrapolationshort/(60.d0),maxtimeextrapolationshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME SCALESYMSHORT        ',&
            timescalesymshort/60.d0,timescalesymshort
          write(ounit,'(a35,x,2f8.2)')'     average SCALESYMSHORT        ',&
            timescalesymshort/(60.d0*dble(mpisize)),timescalesymshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum SCALESYMSHORT        ',&
            maxtimescalesymshort/(60.d0),maxtimescalesymshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME SCALEDSFUNCSHORT     ',&
            timescaledsfuncshort/60.d0,timescaledsfuncshort
          write(ounit,'(a35,x,2f8.2)')'     average SCALEDSFUNCSHORT     ',&
            timescaledsfuncshort/(60.d0*dble(mpisize)),timescaledsfuncshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum SCALEDSFUNCSHORT     ',&
            maxtimescaledsfuncshort/(60.d0),maxtimescaledsfuncshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME ESHORT               ',&
            timeeshort/60.d0,timeeshort
          write(ounit,'(a35,x,2f8.2)')'     average ESHORT               ',&
            timeeshort/(60.d0*dble(mpisize)),timeeshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum ESHORT               ',&
            maxtimeeshort/(60.d0),maxtimeeshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME FSHORT               ',&
            timefshort/60.d0,timefshort
          write(ounit,'(a35,x,2f8.2)')'     average FSHORT               ',&
            timefshort/(60.d0*dble(mpisize)),timefshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum FSHORT               ',&
            maxtimefshort/(60.d0),maxtimefshort
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME STRESS SHORT         ',&
            timesshort/60.d0,timesshort
          write(ounit,'(a35,x,2f8.2)')'     average STRESS SHORT         ',&
            timesshort/(60.d0*dble(mpisize)),timesshort/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum STRESS SHORT         ',&
            maxtimesshort/(60.d0),maxtimesshort
          write(ounit,'(a35,x,2f8.2)')'   TOTTIME ELEC                   ',&
            timeelec/60.d0,timeelec
          write(ounit,'(a35,x,2f8.2)')'   average ELEC                   ',&
            timeelec/(60.d0*dble(mpisize)),timeelec/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'   maximum ELEC                   ',&
            maxtimeelec/(60.d0),maxtimeelec
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME SYMELEC1             ',&
            timesymelec1/60.d0,timesymelec1
          write(ounit,'(a35,x,2f8.2)')'     average SYMELEC1             ',&
            timesymelec1/(60.d0*dble(mpisize)),timesymelec1/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum SYMELEC1             ',&
            maxtimesymelec1/(60.d0),maxtimesymelec1
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME CHARGES              ',&
            timecharge/60.d0,timecharge
          write(ounit,'(a35,x,2f8.2)')'     average CHARGES              ',&
            timecharge/(60.d0*dble(mpisize)),timecharge/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum CHARGES              ',&
            maxtimecharge/(60.d0),maxtimecharge
          write(ounit,'(a35,x,2f8.2)')'   TOTTIME COMMUNICATION 1        ',&
            timecomm1/60.d0,timecomm1
          write(ounit,'(a35,x,2f8.2)')'   average COMMUNICATION 1        ',&
            timecomm1/(60.d0*dble(mpisize)),timecomm1/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'   maximum COMMUNICATION 1        ',&
            maxtimecomm1/(60.d0),maxtimecomm1
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME SYMELEC2             ',&
            timesymelec1/60.d0,timesymelec2
          write(ounit,'(a35,x,2f8.2)')'     average SYMELEC2             ',&
            timesymelec2/(60.d0*dble(mpisize)),timesymelec2/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum SYMELEC2             ',&
            maxtimesymelec2/(60.d0),maxtimesymelec2
          write(ounit,'(a35,x,2f8.2)')'     TOTTIME EELEC                ',&
            timeeelec/60.d0,timeeelec
          write(ounit,'(a35,x,2f8.2)')'     average EELEC                ',&
            timeeelec/(60.d0*dble(mpisize)),timeeelec/dble(mpisize)
          write(ounit,'(a35,x,2f8.2)')'     maximum EELEC                ',&
            maxtimeeelec/(60.d0),maxtimeeelec
        endif
!!
        timefinalize=timefinalizeend-timefinalizestart
          write(ounit,'(a35,x,2f8.2)')' TOTTIME FINALIZE                 ',&
          timefinalize/60.d0,timefinalize
!!
      endif ! mpirank.eq.0
!!
      deallocate(time_array_local)
!!
      return
      end
