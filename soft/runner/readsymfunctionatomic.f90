!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine for all atom-centered symmetry functions

!! called by:
!! - readinput.f90
!!
      subroutine readsymfunctionatomic(keyword,&
        maxnum_funcvalues_local,symcount_local,function_type_local,symelement_local,&
        funccutoff_local,eta_local,zeta_local,rshift_local,lambda_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer ztemp                                              ! internal
      integer itemp                                              ! internal
      integer ztemp1                                             ! internal
      integer ztemp2                                             ! internal
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
!!
      character*2 elementtemp                                    ! internal
      character*2 elementtemp1                                   ! internal
      character*2 elementtemp2                                   ! internal
      character*40 dummy                                         ! internal
      character*40 keyword                                       ! in
!!
      backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp
        call checkelement(elementtemp)
        call nuccharge(elementtemp,ztemp)
        backspace(nnunit)
        symcount_local(elementindex(ztemp))=symcount_local(elementindex(ztemp))+1
        read(nnunit,*,ERR=99)dummy,elementtemp,function_type_temp
!! type 1 radial function
        if(function_type_temp.eq.1)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_temp,elementtemp1
          call checkelement(elementtemp1)
          call nuccharge(elementtemp1,ztemp1)
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            elementtemp1,&
            funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=ztemp1
!! type 2 radial function
        elseif(function_type_temp.eq.2)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_temp,elementtemp1
          call checkelement(elementtemp1)
          call nuccharge(elementtemp1,ztemp1)
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            elementtemp1,&
            eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            rshift_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=ztemp1
!! type 3 angular function
        elseif(function_type_temp.eq.3)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_temp,elementtemp1,elementtemp2
          call checkelement(elementtemp1)
          call checkelement(elementtemp2)
          call nuccharge(elementtemp1,ztemp1)
          call nuccharge(elementtemp2,ztemp2)
!! sort numerically
          if(ztemp1.gt.ztemp2)then
            itemp=ztemp1
            ztemp1=ztemp2
            ztemp2=itemp
          endif
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            elementtemp1,elementtemp2,&
            eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            lambda_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            zeta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=ztemp1
          symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=ztemp2
!! type 4 radial function
        elseif(function_type_temp.eq.4)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_temp,elementtemp1
          call checkelement(elementtemp1)
          call nuccharge(elementtemp1,ztemp1)
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            elementtemp1,&
            eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=ztemp1
!! type 5 just cartesian coordinate
        elseif(function_type_temp.eq.5)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
!! type 6 just R_ij
        elseif(function_type_temp.eq.6)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_temp,elementtemp1
          call checkelement(elementtemp1)
          call nuccharge(elementtemp1,ztemp1)
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            elementtemp1,&
            funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=ztemp1
!! type 8 angular function
        elseif(function_type_temp.eq.8)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_temp,elementtemp1,elementtemp2
          call checkelement(elementtemp1)
          call checkelement(elementtemp2)
          call nuccharge(elementtemp1,ztemp1)
          call nuccharge(elementtemp2,ztemp2)
!! sort numerically
          if(ztemp1.gt.ztemp2)then
            itemp=ztemp1
            ztemp1=ztemp2
            ztemp2=itemp
          endif
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            elementtemp1,elementtemp2,&
            eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            rshift_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=ztemp1
          symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=ztemp2
!! type 9 angular function
        elseif(function_type_temp.eq.9)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_temp,elementtemp1,elementtemp2
          call checkelement(elementtemp1)
          call checkelement(elementtemp2)
          call nuccharge(elementtemp1,ztemp1)
          call nuccharge(elementtemp2,ztemp2)
!! sort numerically
          if(ztemp1.gt.ztemp2)then
            itemp=ztemp1
            ztemp1=ztemp2
            ztemp2=itemp
          endif
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,&
            function_type_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            elementtemp1,elementtemp2,&
            eta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            lambda_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            zeta_local(symcount_local(elementindex(ztemp)),elementindex(ztemp)),&
            funccutoff_local(symcount_local(elementindex(ztemp)),elementindex(ztemp))
          symelement_local(symcount_local(elementindex(ztemp)),1,elementindex(ztemp))=ztemp1
          symelement_local(symcount_local(elementindex(ztemp)),2,elementindex(ztemp))=ztemp2
        else
          write(ounit,*)'Error: Unknown symfunction type in input.nn ',&
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
