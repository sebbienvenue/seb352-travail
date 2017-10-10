!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
        subroutine gettime(timesec,timemin,timehours)
        implicit none
        real*8 timetot
        real*8 timeold
        real*8 timediff
        real*8 timesec,timemin,timehours
        real*8 seconds
        integer day,dayold
        integer nhours,nminutes
        integer icon
!!
        save timeold,dayold
!!
        character*10 date,time
!!
        call date_and_time(date,time)
!!        write(*,*)'date,time', date,time
!!
        read(time,'(i2,i2,f6.3)')nhours,nminutes,seconds
!!        write(*,*)'nhours,nminutes,seconds',nhours,nminutes,seconds
!!
        timetot=3600.d0*dble(nhours)+60.d0*dble(nminutes)+seconds
!!        write(*,*)'timetot ',timetot
!!
        read(date,'(6x,i2)') day
!!        write(*,*)'day ',day
!!
        timediff = timetot - timeold
        if(day.ne.dayold) timediff = timediff + 3600.d0*24.d0
        timeold = timetot
        dayold = day
!!        write(*,*)'timediff ',timediff
        timesec=timediff
        timemin=timesec/60.d0
        timehours=timemin/60.d0
!!
        return
        end
!!
