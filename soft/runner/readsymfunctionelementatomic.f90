!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine for all atom-centered symmetry functions

!! called by:
!! - readinput.f90
!!
      subroutine readsymfunctionelementatomic(keyword,&
        maxnum_funcvalues_local,symcount_local,function_type_local,symelement_local,&
        funccutoff_local,eta_local,zeta_local,rshift_local,lambda_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2                                              ! internal
      integer ztemp                                              ! internal
      integer symcount_local(nelem)                              ! in/out
      integer maxnum_funcvalues_local                            ! in
      integer function_type_temp                                 ! internal
      integer function_type_local(maxnum_funcvalues_local,nelem) ! out
      integer symelement_local(maxnum_funcvalues_local,2,nelem)  ! out
!!
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)     ! out
      real*8 eta_local(maxnum_funcvalues_local,nelem)            ! out
      real*8 zeta_local(maxnum_funcvalues_local,nelem)           ! out
      real*8 rshift_local(maxnum_funcvalues_local,nelem)         ! out
      real*8 lambda_local(maxnum_funcvalues_local,nelem)         ! out
      real*8 eta_temp                                            ! internal
      real*8 zeta_temp                                           ! internal
      real*8 rshift_temp                                         ! internal
      real*8 lambda_temp                                         ! internal
      real*8 funccutoff_temp                                     ! internal
!!
      character*2 elementtemp                                    ! internal
      character*40 dummy                                         ! internal
      character*40 keyword                                       ! in
!!
      backspace(nnunit)
      read(nnunit,*,ERR=99)dummy,elementtemp
      call checkelement(elementtemp)
      call nuccharge(elementtemp,ztemp)
      backspace(nnunit)
      read(nnunit,*,ERR=99)dummy,elementtemp,function_type_temp
!! type 1 radial function
      if(function_type_temp.eq.1)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,&
          function_type_temp,funccutoff_temp
        do i1=1,nelem
          symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
          function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
          funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i1)
        enddo ! i1
!! type 2 radial function
      elseif(function_type_temp.eq.2)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,&
          function_type_temp,eta_temp,rshift_temp,funccutoff_temp
        do i1=1,nelem
          symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
          function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
          eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
          rshift_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=rshift_temp
          funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i1)
        enddo ! i1
!! type 3 angular function
      elseif(function_type_temp.eq.3)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,&
          function_type_temp,eta_temp,lambda_temp,&
          zeta_temp,funccutoff_temp
        do i1=1,nelem
          symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
          function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
          eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
          lambda_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=lambda_temp
          zeta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=zeta_temp
          funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i1)
          symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=nucelem(i1)
        enddo
        do i1=1,nelem
          if(nelem.gt.1)then
            do i2=1,i1-1
              symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
              function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
              eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
              lambda_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=lambda_temp
              zeta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=zeta_temp
              funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
              symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i2)
              symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=nucelem(i1)
            enddo ! i2
          endif
        enddo ! i1
!! type 4 radial function
      elseif(function_type_temp.eq.4)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,&
          function_type_temp,eta_temp,funccutoff_temp
        do i1=1,nelem
          symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
          function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
          eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
          funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i1)
        enddo ! i1
!! type 5 just cartesian coordinate
      elseif(function_type_temp.eq.5)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,function_type_temp,eta_temp
        symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
        function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
        eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
!! type 6 just R_ij
      elseif(function_type_temp.eq.6)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,&
          function_type_temp,funccutoff_temp
        symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
        function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
!! type 8 angular function
      elseif(function_type_temp.eq.8)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,&
          function_type_temp,eta_temp,rshift_temp,&
          funccutoff_temp
        do i1=1,nelem
          symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
          function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
          eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
          rshift_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=rshift_temp
          funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i1)
          symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=nucelem(i1)
        enddo
        do i1=1,nelem
          if(nelem.gt.1)then
            do i2=1,i1-1
              symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
              function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
              eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
              rshift_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=rshift_temp
              funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
              symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i2)
              symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=nucelem(i1)
            enddo ! i2
          endif
        enddo ! i1
!! type 9 angular function
      elseif(function_type_temp.eq.9)then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,&
          function_type_temp,eta_temp,lambda_temp,&
          zeta_temp,funccutoff_temp
        do i1=1,nelem
          symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
          function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
          eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
          lambda_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=lambda_temp
          zeta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=zeta_temp
          funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i1)
          symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=nucelem(i1)
        enddo
        do i1=1,nelem
          if(nelem.gt.1)then
            do i2=1,i1-1
              symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
              function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=function_type_temp
              eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=eta_temp
              lambda_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=lambda_temp
              zeta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=zeta_temp
              funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))=funccutoff_temp
              symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=nucelem(i2)
              symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=nucelem(i1)
            enddo ! i2
          endif
        enddo ! i1
      else
        write(ounit,*)'Error: Unknown element_symfunction_short type in input.nn ',&
         function_type_temp !'
        stop
      endif ! function_type_temp
!!
      return
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop
!!
      end
!!
