!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine for all atom-centered symmetry functions

!! called by:
!! - readinput.f90
!!
      subroutine readsymfunctionhextoff(keyword,ndim,&
        maxnum_funcvalues_local,symcount_local,function_type_local,symelement_local,&
        funccutoff_local,eta_local,zeta_local,rshift_local,lambda_local)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2,i3,ndim
      integer icount                                             ! internal
      integer symcount_local(ndim)                              ! in/out
      integer maxnum_funcvalues_local                            ! in
      integer function_type_temp                                 ! internal
      integer function_type_local(maxnum_funcvalues_local,ndim) ! out
      integer symelement_local(maxnum_funcvalues_local,2,ndim)  ! out
!!
      real*8 funccutoff_local(maxnum_funcvalues_local,ndim)     ! out
      real*8 eta_local(maxnum_funcvalues_local,ndim)            ! out
      real*8 zeta_local(maxnum_funcvalues_local,ndim)           ! out
      real*8 rshift_local(maxnum_funcvalues_local,ndim)         ! out
      real*8 lambda_local(maxnum_funcvalues_local,ndim)         ! out
      real*8 funccutoff_temp                                     ! internal
      real*8 lambda_temp                                         ! internal
      real*8 rshift_temp                                         ! internal
      real*8 eta_temp                                            ! internal
      real*8 zeta_temp                                           ! internal
!!
      character*2 elementtemp1,elementtemp2,elementtemp3         ! internal
      integer ztemp1,ztemp2,ztemp3                               ! internal
      character*40 dummy                                         ! internal
      character*40 keyword                                       ! in
!!
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,elementtemp3,function_type_temp
!!        write(*,*) 'function type', function_type_temp
        call checkelement(elementtemp1)
        call checkelement(elementtemp2)
        call checkelement(elementtemp3)
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        call nuccharge(elementtemp3,ztemp3)
        icount=ndim
!!
!! type 1
        if(function_type_temp.eq.1)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            elementtemp1,elementtemp2,elementtemp3,function_type_temp,eta_temp,funccutoff_temp
          symcount_local(icount)=symcount_local(icount)+1
          function_type_local(symcount_local(icount),icount)=function_type_temp
          funccutoff_local(symcount_local(icount),icount)   =funccutoff_temp
          eta_local(symcount_local(icount),icount)=eta_temp
          
!! type 2 
        elseif(function_type_temp.eq.2)then ! just bond length
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            elementtemp1,elementtemp2,elementtemp3,function_type_temp,eta_temp,funccutoff_temp
          symcount_local(icount)=symcount_local(icount)+1
          function_type_local(symcount_local(icount),icount)=function_type_temp
          funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
          eta_local(symcount_local(icount),icount)=eta_temp
!! type 3 
        elseif(function_type_temp.eq.3)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            elementtemp1,elementtemp2,elementtemp3,function_type_temp,eta_temp,funccutoff_temp
          symcount_local(icount)=symcount_local(icount)+1
          function_type_local(symcount_local(icount),icount)=function_type_temp
          funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
          eta_local(symcount_local(icount),icount)=eta_temp
!!
        else
          write(ounit,*)'Error: Unknown global_pairsymfunction_short type in input.nn ',&
           function_type_temp !'
          stop
        endif
!!
!!      write(*,*) 'eta hextoff', eta_local
      return
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments hextoff'
      write(ounit,*)dummy, function_type_temp
      stop
!!
      end
!!
