!####################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine for all atom-centered symmetry functions

!! called by:
!! - readinput.f90
!!
      subroutine readsymfunctionelementpair(keyword,&
        maxnum_funcvalues_local,symcount_local,function_type_local,symelement_local,&
        funccutoff_local,eta_local,zeta_local,rshift_local,lambda_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3
      integer icount                                             ! internal
      integer ztemp1                                             ! internal
      integer ztemp2                                             ! internal
      integer symcount_local(npairs)                              ! in/out
      integer maxnum_funcvalues_local                            ! in
      integer function_type_temp                                 ! internal
      integer function_type_local(maxnum_funcvalues_local,npairs) ! out
      integer symelement_local(maxnum_funcvalues_local,2,npairs)  ! out
!!
      real*8 funccutoff_local(maxnum_funcvalues_local,npairs)     ! out
      real*8 eta_local(maxnum_funcvalues_local,npairs)            ! out
      real*8 zeta_local(maxnum_funcvalues_local,npairs)           ! out
      real*8 rshift_local(maxnum_funcvalues_local,npairs)         ! out
      real*8 lambda_local(maxnum_funcvalues_local,npairs)         ! out
      real*8 funccutoff_temp                                     ! internal
      real*8 lambda_temp                                         ! internal
      real*8 rshift_temp                                         ! internal
      real*8 eta_temp                                            ! internal
      real*8 zeta_temp                                           ! internal
!!
      character*2 elementtemp1                                   ! internal
      character*2 elementtemp2                                   ! internal
      character*40 dummy                                         ! internal
      character*40 keyword                                       ! in
!!
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,function_type_temp
        call checkelement(elementtemp1)
        call checkelement(elementtemp2)
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
!! type 1 
        if(function_type_temp.eq.1)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,&
            function_type_temp,funccutoff_temp
          icount=0
          do i1=1,nelem
            do i2=i1,nelem
              icount=icount+1
              if((ztemp1.eq.nucelem(i1)).and.(ztemp2.eq.nucelem(i2)))then
                do i3=1,nelem ! loop over neighbor elements of pair
                  symcount_local(icount)=symcount_local(icount)+1
                  function_type_local(symcount_local(icount),icount) =function_type_temp
                  funccutoff_local(symcount_local(icount),icount)    =funccutoff_temp
                  symelement_local(symcount_local(icount),1,icount)  =nucelem(i3)  ! neighbor element
                enddo ! i3
              elseif((ztemp1.eq.nucelem(i2)).and.(ztemp2.eq.nucelem(i1)))then
                do i3=1,nelem ! loop over neighbor elements of pair
                  symcount_local(icount)=symcount_local(icount)+1
                  function_type_local(symcount_local(icount),icount) =function_type_temp
                  funccutoff_local(symcount_local(icount),icount)    =funccutoff_temp
                  symelement_local(symcount_local(icount),1,icount)  =nucelem(i3)  ! neighbor element
                enddo ! i3
              endif
            enddo ! i2
          enddo ! i1
!! type 2 
        elseif(function_type_temp.eq.2)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,&
            function_type_temp,funccutoff_temp
!! in this special case this is the same as pairsymfunction_short
          icount=0
          do i1=1,nelem
            do i2=i1,nelem
              icount=icount+1
              if((ztemp1.eq.nucelem(i1)).and.(ztemp2.eq.nucelem(i2)))then
                symcount_local(icount)=symcount_local(icount)+1
                function_type_local(symcount_local(icount),icount)=function_type_temp
                funccutoff_local(symcount_local(icount),icount)   =funccutoff_temp
              elseif((ztemp1.eq.nucelem(i2)).and.(ztemp2.eq.nucelem(i1)))then
                symcount_local(icount)=symcount_local(icount)+1
                function_type_local(symcount_local(icount),icount)=function_type_temp
                funccutoff_local(symcount_local(icount),icount)   =funccutoff_temp
              endif
            enddo ! i2
          enddo ! i1
!! type 3 
         elseif(function_type_temp.eq.3)then
           backspace(nnunit)
           read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,&
             function_type_temp,eta_temp,rshift_temp,funccutoff_temp
           icount=0
           do i1=1,nelem
             do i2=i1,nelem
               icount=icount+1
               if((ztemp1.eq.nucelem(i1)).and.(ztemp2.eq.nucelem(i2)))then
                  do i3=1,nelem ! loop over neighbor elements of pair
                   symcount_local(icount)=symcount_local(icount)+1
                   function_type_local(symcount_local(icount),icount)=function_type_temp
                   funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
                   symelement_local(symcount_local(icount),1,icount) =nucelem(i3)  ! neighbor element
                   eta_local(symcount_local(icount),icount)=eta_temp
                   rshift_local(symcount_local(icount),icount)=rshift_temp
                 enddo ! i3
               elseif((ztemp1.eq.nucelem(i2)).and.(ztemp2.eq.nucelem(i1)))then
                 do i3=1,nelem ! loop over neighbor elements of pair
                   symcount_local(icount)=symcount_local(icount)+1
                   function_type_local(symcount_local(icount),icount)=function_type_temp
                   funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
                   symelement_local(symcount_local(icount),1,icount) =nucelem(i3)  ! neighbor element
                   eta_local(symcount_local(icount),icount)=eta_temp
                   rshift_local(symcount_local(icount),icount)=rshift_temp
                 enddo ! i3
               endif
             enddo ! i2
           enddo ! i1
!! type 4 
        elseif(function_type_temp.eq.4)then
           backspace(nnunit)
           read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,function_type_temp,eta_temp,lambda_temp,zeta_temp,funccutoff_temp
           icount=0
           do i1=1,nelem
             do i2=i1,nelem
               icount=icount+1
               if((ztemp1.eq.nucelem(i1)).and.(ztemp2.eq.nucelem(i2)))then
                  do i3=1,nelem ! loop over neighbor elements of pair
                   symcount_local(icount)=symcount_local(icount)+1
                   function_type_local(symcount_local(icount),icount)=function_type_temp
                   funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
                   symelement_local(symcount_local(icount),1,icount) =nucelem(i3)  ! neighbor element
                   eta_local(symcount_local(icount),icount)=eta_temp
                   lambda_local(symcount_local(icount),icount)=lambda_temp
                   zeta_local(symcount_local(icount),icount)=zeta_temp
                 enddo ! i3
               elseif((ztemp1.eq.nucelem(i2)).and.(ztemp2.eq.nucelem(i1)))then
                 do i3=1,nelem ! loop over neighbor elements of pair
                   symcount_local(icount)=symcount_local(icount)+1
                   function_type_local(symcount_local(icount),icount)=function_type_temp
                   funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
                   symelement_local(symcount_local(icount),1,icount) =nucelem(i3)  ! neighbor element
                   eta_local(symcount_local(icount),icount)=eta_temp
                   lambda_local(symcount_local(icount),icount)=lambda_temp
                   zeta_local(symcount_local(icount),icount)=zeta_temp
                 enddo ! i3
               endif
             enddo ! i2
           enddo ! i1
!! type 5 
        elseif(function_type_temp.eq.5)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,&
          function_type_temp,funccutoff_temp
          icount=0
          do i1=1,nelem
            do i2=i1,nelem
              icount=icount+1
              if((ztemp1.eq.nucelem(i1)).and.(ztemp2.eq.nucelem(i2)))then
                symcount_local(icount)=symcount_local(icount)+1
                function_type_local(symcount_local(icount),icount)=function_type_temp
                funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
              elseif((ztemp1.eq.nucelem(i2)).and.(ztemp2.eq.nucelem(i1)))then
                symcount_local(icount)=symcount_local(icount)+1
                function_type_local(symcount_local(icount),icount)=function_type_temp
                funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
              endif
            enddo ! i2
          enddo ! i1
!! type 6 
        elseif(function_type_temp.eq.6)then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,&
          function_type_temp,eta_temp,rshift_temp,funccutoff_temp
          icount=0
          do i1=1,nelem
            do i2=i1,nelem
              icount=icount+1
              if((ztemp1.eq.nucelem(i1)).and.(ztemp2.eq.nucelem(i2)))then
                symcount_local(icount)=symcount_local(icount)+1
                function_type_local(symcount_local(icount),icount)=function_type_temp
                eta_local(symcount_local(icount),icount)=eta_temp
                rshift_local(symcount_local(icount),icount)=rshift_temp
                funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
              elseif((ztemp1.eq.nucelem(i2)).and.(ztemp2.eq.nucelem(i1)))then
                symcount_local(icount)=symcount_local(icount)+1
                function_type_local(symcount_local(icount),icount)=function_type_temp
                eta_local(symcount_local(icount),icount)=eta_temp
                rshift_local(symcount_local(icount),icount)=rshift_temp
                funccutoff_local(symcount_local(icount),icount)=funccutoff_temp
              endif
            enddo ! i2
          enddo ! i1
!!
        else
          write(ounit,*)'Error: Unknown element_pairsymfunction_short type in input.nn ',&
           function_type_temp !'
          stop
        endif
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
