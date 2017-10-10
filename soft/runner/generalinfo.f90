!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! -main.f90
!!
      subroutine generalinfo()
!!
      use fileunits
!!
      implicit none
!!
      integer nhours,nminutes,day,year,month

!!
!! CHANGE ANDI: GFORTRAN: need return value for gfortran hostnm and getcwd functions
!!
      integer ierr
!! END CHANGE

!!
      real*8 seconds
!!
      character*20 machname
      character*20 hostnm 
      character*60 directory 
      character*60 getcwd 
      character*20 username 
      character*10 date,time
!!
!!      call hostnm(machname,status) ! breda gets trouble with this
!!

!!
!! CHANGE ANDI: GFORTRAN: ifort and gfortran have different functions for machine name and current directory
!!                        to use these preprocessor directives compile this file with -cpp option
!!
     !machname=hostnm()
     !call getlog(username)
     !call date_and_time(date,time)
     !read(time,'(i2,i2,f6.3)')nhours,nminutes,seconds
     !read(date,'(i4,i2,i2)') year,month,day
     !directory=getcwd()
#ifdef __INTEL_COMPILER
      machname=hostnm()
      directory=getcwd()
#else
      ierr=hostnm(machname)
      ierr=getcwd(directory)
#endif

      call getlog(username)
      call date_and_time(date,time)
      read(time,'(i2,i2,f6.3)')nhours,nminutes,seconds
      read(date,'(i4,i2,i2)') year,month,day
!! END CHANGE

!!
      write(ounit,*)'General job information:'
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Executing host    : ',machname
      write(ounit,*)'User name         : ',username
      write(ounit,'(a,i2,a,i2,a,i4)')' Starting date     : ',day,'.',month,'.',year
      write(ounit,'(a,i2,a,i2,a,i2)')' Starting time     : ',nhours,' h ',nminutes,' min '
      write(ounit,'(2a)')' Working directory : ',directory
!!      write(ounit,'(4a)')' Job started on ',machname,' by user ',username
!!      write(ounit,'(a,i2,a,i3,2a,i2,a,i2,a,i4)')&
!!        ' Starting at ',nhours,' h ',nminutes,' min ',&
!!        ' on ',day,'.',month,'.',year
!!      write(ounit,'(2a)')' in directory ',directory
!!
      write(ounit,*)'-------------------------------------------------------------'
      return
      end
