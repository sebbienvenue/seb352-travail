!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine for all atom-centered symmetry functions

!! called by:
!! - readinput.f90
!!
      subroutine readsymfunctionglobalatomic(keyword,&
        maxnum_funcvalues_local,symcount_local,function_type_local,symelement_local,&
        funccutoff_local,eta_local,zeta_local,rshift_local,lambda_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3
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
      real*8 funccutoff_temp                                     ! internal
      real*8 rshift_temp                                         ! internal
      real*8 eta_temp                                            ! internal
      real*8 zeta_temp                                           ! internal
      real*8 lambda_temp                                         ! internal
!!
      character*40 dummy                                         ! internal
      character*40 keyword                                       ! in
!!
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_temp
!! type 1 radial function
        if(function_type_temp.eq.1)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            function_type_temp,funccutoff_temp
          do i1=1,nelem
            do i2=1,nelem
              symcount_local(i1)=symcount_local(i1)+1
              function_type_local(symcount_local(i1),i1)=function_type_temp
              funccutoff_local(symcount_local(i1),i1)   =funccutoff_temp
              symelement_local(symcount_local(i1),1,i1) =nucelem(i2)
            enddo ! i2
          enddo ! i1
!! type 2 radial function
        elseif(function_type_temp.eq.2)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            function_type_temp,eta_temp,rshift_temp,funccutoff_temp
          do i1=1,nelem
            do i2=1,nelem
              symcount_local(i1)=symcount_local(i1)+1
              function_type_local(symcount_local(i1),i1)=function_type_temp
              eta_local(symcount_local(i1),i1)          =eta_temp
              rshift_local(symcount_local(i1),i1)       =rshift_temp
              funccutoff_local(symcount_local(i1),i1)   =funccutoff_temp
              symelement_local(symcount_local(i1),1,i1) =nucelem(i2)
            enddo ! i2
          enddo ! i1
!! type 3 angular function
        elseif(function_type_temp.eq.3)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            function_type_temp,eta_temp,lambda_temp,&
            zeta_temp,funccutoff_temp
          do i3=1,nelem
            do i1=1,nelem
              symcount_local(i3)=symcount_local(i3)+1
              function_type_local(symcount_local(i3),i3)=function_type_temp
              eta_local(symcount_local(i3),i3)          =eta_temp
              lambda_local(symcount_local(i3),i3)       =lambda_temp
              zeta_local(symcount_local(i3),i3)         =zeta_temp
              funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
              symelement_local(symcount_local(i3),1,i3) =nucelem(i1)
              symelement_local(symcount_local(i3),2,i3) =nucelem(i1)
            enddo
            do i1=1,nelem
              if(nelem.gt.1)then
                do i2=1,i1-1
                  symcount_local(i3)=symcount_local(i3)+1
                  function_type_local(symcount_local(i3),i3)=function_type_temp
                  eta_local(symcount_local(i3),i3)          =eta_temp
                  lambda_local(symcount_local(i3),i3)       =lambda_temp
                  zeta_local(symcount_local(i3),i3)         =zeta_temp
                  funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
                  symelement_local(symcount_local(i3),1,i3) =nucelem(i2)
                  symelement_local(symcount_local(i3),2,i3) =nucelem(i1)
                enddo ! i2
              endif
            enddo ! i1
          enddo ! i3
!! type 4 radial function
        elseif(function_type_temp.eq.4)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            function_type_temp,eta_temp,funccutoff_temp
          do i3=1,nelem
            do i1=1,nelem
              symcount_local(i3)=symcount_local(i3)+1
              function_type_local(symcount_local(i3),i3)=function_type_temp
              eta_local(symcount_local(i3),i3)          =eta_temp
              funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
              symelement_local(symcount_local(i3),1,i3) =nucelem(i1)
            enddo ! i1
        enddo ! i3
!! type 5 just cartesian coordinate
        elseif(function_type_temp.eq.5)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,function_type_temp,eta_temp
          do i3=1,nelem
            symcount_local(i3)=symcount_local(i3)+1
            function_type_local(symcount_local(i3),i3)=function_type_temp
            eta_local(symcount_local(i3),i3)          =eta_temp
            symelement_local(symcount_local(i3),1,i3) =nucelem(i3)
          enddo ! i3
!! type 6 just R_ij
        elseif(function_type_temp.eq.6)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            function_type_temp,funccutoff_temp
          do i3=1,nelem
            symcount_local(i3)=symcount_local(i3)+1
            function_type_local(symcount_local(i3),i3)=function_type_temp
            funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
            symelement_local(symcount_local(i3),1,i3) =nucelem(i3)
          enddo ! i3
!! type 8 angular function
        elseif(function_type_temp.eq.8)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            function_type_temp,eta_temp,rshift_temp,&
            funccutoff_temp
          do i3=1,nelem
            do i1=1,nelem
              symcount_local(i3)=symcount_local(i3)+1
              function_type_local(symcount_local(i3),i3)=function_type_temp
              eta_local(symcount_local(i3),i3)          =eta_temp
              rshift_local(symcount_local(i3),i3)       =rshift_temp
              funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
              symelement_local(symcount_local(i3),1,i3) =nucelem(i1)
              symelement_local(symcount_local(i3),2,i3) =nucelem(i1)
            enddo
            do i1=1,nelem
              if(nelem.gt.1)then
                do i2=1,i1-1
                  symcount_local(i3)=symcount_local(i3)+1
                  function_type_local(symcount_local(i3),i3)=function_type_temp
                  eta_local(symcount_local(i3),i3)          =eta_temp
                  rshift_local(symcount_local(i3),i3)       =rshift_temp
                  funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
                  symelement_local(symcount_local(i3),1,i3) =nucelem(i2)
                  symelement_local(symcount_local(i3),2,i3) =nucelem(i1)
                enddo ! i2
              endif
            enddo ! i1
          enddo ! i3
!! type 9 angular function
        elseif(function_type_temp.eq.9)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,&
            function_type_temp,eta_temp,lambda_temp,&
            zeta_temp,funccutoff_temp
          do i3=1,nelem
            do i1=1,nelem
              symcount_local(i3)=symcount_local(i3)+1
              function_type_local(symcount_local(i3),i3)=function_type_temp
              eta_local(symcount_local(i3),i3)          =eta_temp
              lambda_local(symcount_local(i3),i3)       =lambda_temp
              zeta_local(symcount_local(i3),i3)         =zeta_temp
              funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
              symelement_local(symcount_local(i3),1,i3) =nucelem(i1)
              symelement_local(symcount_local(i3),2,i3) =nucelem(i1)
            enddo
            do i1=1,nelem
              if(nelem.gt.1)then
                do i2=1,i1-1
                  symcount_local(i3)=symcount_local(i3)+1
                  function_type_local(symcount_local(i3),i3)=function_type_temp
                  eta_local(symcount_local(i3),i3)          =eta_temp
                  lambda_local(symcount_local(i3),i3)       =lambda_temp
                  zeta_local(symcount_local(i3),i3)         =zeta_temp
                  funccutoff_local(symcount_local(i3),i3)   =funccutoff_temp
                  symelement_local(symcount_local(i3),1,i3) =nucelem(i2)
                  symelement_local(symcount_local(i3),2,i3) =nucelem(i1)
                enddo ! i2
              endif
            enddo ! i1
          enddo ! i3
        else
          write(ounit,*)'Error: Unknown global_symfunction_short type in input.nn ',&
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
