!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - readinput.f90
!!
      subroutine sortpairsymfunctions(npairs,&
        maxnum_funcvaluesp,num_funcvaluesp,&
        function_typep,sympelement,&
        etap,zetap,rshiftp,lambdap,funccutoffp,&
        elempair)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2
      integer npairs                                    ! in
      integer maxnum_funcvaluesp                        ! in
      integer num_funcvaluesp(npairs)                   ! in
      integer function_type_temp                        ! internal
      integer function_typep(maxnum_funcvaluesp,npairs) ! in/out
      integer sympelement(maxnum_funcvaluesp,2,npairs)  ! in/out
      integer sympelement_temp(2)                       ! internal
      integer counter(6,npairs)                         ! internal
      integer elempair(npairs,2)                        ! in 
!!
      real*8 funccutoffp(maxnum_funcvaluesp,npairs)     ! in/out
      real*8 etap(maxnum_funcvaluesp,npairs)               ! in/out
      real*8 zetap(maxnum_funcvaluesp,npairs)              ! in/out
      real*8 lambdap(maxnum_funcvaluesp,npairs)            ! in/out
      real*8 rshiftp(maxnum_funcvaluesp,npairs)            ! in/out
      real*8 funccutoff_temp                            ! internal
      real*8 eta_temp                                   ! internal
      real*8 zeta_temp                                  ! internal
      real*8 lambda_temp                                ! internal
      real*8 rshift_temp                                ! internal
      real*8 thres                                      ! internal
!!
!!
      thres=0.0001d0
!!
!! check if there are unkown symmetry function types
      do i1=1,npairs
        do i2=1,num_funcvaluesp(i1)
          if(function_typep(i2,i1).gt.6)then
            write(ounit,*)'ERROR in sortpairsymfunctions'
            write(ounit,*)'unknown function_type ',function_typep(i2,i1)
            write(ounit,*)'redimension also counter!'
            stop
          endif
        enddo
      enddo
!!
!! first sort according to function_type
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 10       continue
          do i2=1,num_funcvaluesp(i1)-1
            if(function_typep(i2,i1).gt.function_typep(i2+1,i1))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 10
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! count symmetry functions of each type
!!
      counter(:,:)=0
      do i1=1,npairs
        do i2=1,num_funcvaluesp(i1)
          counter(function_typep(i2,i1),i1)=counter(function_typep(i2,i1),i1)+1
        enddo ! i2
      enddo ! i1
!!
!! sort according to cutoff
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 11       continue
          do i2=1,num_funcvaluesp(i1)-1
            if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
              .and.(funccutoffp(i2,i1).gt.funccutoffp(i2+1,i1)))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 11
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to eta
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 12       continue
          do i2=1,num_funcvaluesp(i1)-1
            if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
              .and.(funccutoffp(i2,i1).eq.funccutoffp(i2+1,i1))&
              .and.(etap(i2,i1).gt.etap(i2+1,i1)))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 12
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to zeta
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 13       continue
          do i2=1,num_funcvaluesp(i1)-1
            if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
              .and.(funccutoffp(i2,i1).eq.funccutoffp(i2+1,i1))&
              .and.(etap(i2,i1).eq.etap(i2+1,i1))&
              .and.(zetap(i2,i1).gt.zetap(i2+1,i1)))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 13
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to lambda 
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 14       continue
          do i2=1,num_funcvaluesp(i1)-1
            if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
              .and.(funccutoffp(i2,i1).eq.funccutoffp(i2+1,i1))&
              .and.(etap(i2,i1).eq.etap(i2+1,i1))&
              .and.(zetap(i2,i1).eq.zetap(i2+1,i1))&
              .and.(lambdap(i2,i1).gt.lambdap(i2+1,i1)))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 14
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to rshift 
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 15       continue
          do i2=1,num_funcvaluesp(i1)-1
            if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
              .and.(funccutoffp(i2,i1).eq.funccutoffp(i2+1,i1))&
              .and.(etap(i2,i1).eq.etap(i2+1,i1))&
              .and.(zetap(i2,i1).eq.zetap(i2+1,i1))&
              .and.(lambdap(i2,i1).eq.lambdap(i2+1,i1))&
              .and.(rshiftp(i2,i1).gt.rshiftp(i2+1,i1)))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 15
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to sympelement(:,1,:) 
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 16       continue
          do i2=1,num_funcvaluesp(i1)-1
            if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
              .and.(funccutoffp(i2,i1).eq.funccutoffp(i2+1,i1))&
              .and.(etap(i2,i1).eq.etap(i2+1,i1))&
              .and.(zetap(i2,i1).eq.zetap(i2+1,i1))&
              .and.(lambdap(i2,i1).eq.lambdap(i2+1,i1))&
              .and.(rshiftp(i2,i1).eq.rshiftp(i2+1,i1))&
              .and.(sympelement(i2,1,i1).gt.sympelement(i2+1,1,i1)))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 16
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to sympelement(:,2,:)
!!
      do i1=1,npairs
        if(num_funcvaluesp(i1).gt.1)then
 17       continue
          do i2=1,num_funcvaluesp(i1)-1
            if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
              .and.(funccutoffp(i2,i1).eq.funccutoffp(i2+1,i1))&
              .and.(etap(i2,i1).eq.etap(i2+1,i1))&
              .and.(zetap(i2,i1).eq.zetap(i2+1,i1))&
              .and.(lambdap(i2,i1).eq.lambdap(i2+1,i1))&
              .and.(rshiftp(i2,i1).eq.rshiftp(i2+1,i1))&
              .and.(sympelement(i2,1,i1).eq.sympelement(i2+1,1,i1))&
              .and.(sympelement(i2,2,i1).gt.sympelement(i2+1,2,i1)))then
              function_type_temp    =function_typep(i2,i1)
              function_typep(i2,i1)  =function_typep(i2+1,i1)
              function_typep(i2+1,i1)=function_type_temp
              sympelement_temp(1)    =sympelement(i2,1,i1)
              sympelement(i2,1,i1)   =sympelement(i2+1,1,i1)
              sympelement(i2+1,1,i1) =sympelement_temp(1)
              sympelement_temp(2)    =sympelement(i2,2,i1)
              sympelement(i2,2,i1)   =sympelement(i2+1,2,i1)
              sympelement(i2+1,2,i1) =sympelement_temp(2)
              eta_temp    =etap(i2,i1)
              etap(i2,i1)  =etap(i2+1,i1)
              etap(i2+1,i1)=eta_temp
              zeta_temp    =zetap(i2,i1)
              zetap(i2,i1)  =zetap(i2+1,i1)
              zetap(i2+1,i1)=zeta_temp
              lambda_temp    =lambdap(i2,i1)
              lambdap(i2,i1)  =lambdap(i2+1,i1)
              lambdap(i2+1,i1)=lambda_temp
              rshift_temp    =rshiftp(i2,i1)
              rshiftp(i2,i1)  =rshiftp(i2+1,i1)
              rshiftp(i2+1,i1)=rshift_temp
              funccutoff_temp    =funccutoffp(i2,i1)
              funccutoffp(i2,i1)  =funccutoffp(i2+1,i1)
              funccutoffp(i2+1,i1)=funccutoff_temp
              goto 17
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! check for double symmetry functions
      do i1=1,npairs
        do i2=1,num_funcvaluesp(i1)-1
          if((function_typep(i2,i1).eq.function_typep(i2+1,i1))&
            .and.(abs(funccutoffp(i2,i1)-funccutoffp(i2+1,i1)).lt.thres)&
            .and.(abs(etap(i2,i1)-etap(i2+1,i1)).lt.thres)&
            .and.(abs(zetap(i2,i1)-zetap(i2+1,i1)).lt.thres)&
            .and.(abs(lambdap(i2,i1)-lambdap(i2+1,i1)).lt.thres)&
            .and.(abs(rshiftp(i2,i1)-rshiftp(i2+1,i1)).lt.thres)&
            .and.(sympelement(i2,1,i1).eq.sympelement(i2+1,1,i1))&
            .and.(sympelement(i2,2,i1).eq.sympelement(i2+1,2,i1)))then
            write(ounit,*)'ERROR: a pair symmetry function is specified twice!'
            write(ounit,*)'pair is ',elempair(i1,1),elempair(i1,2)
            write(ounit,*)'function_type ',function_typep(i2,i1)
            write(ounit,'(a14,f10.3)')'funccutoff    ',funccutoffp(i2,i1)
            write(ounit,'(a14,f10.3)')'eta           ',etap(i2,i1)
            write(ounit,'(a14,f10.3)')'zeta          ',zetap(i2,i1)
            write(ounit,'(a14,f10.3)')'lambda        ',lambdap(i2,i1)
            write(ounit,'(a14,f10.3)')'rshift        ',rshiftp(i2,i1)
            write(ounit,*)'sympelement(1) ',sympelement(i2,1,i1)
            write(ounit,*)'sympelement(2) ',sympelement(i2,2,i1)
            stop
          endif
        enddo ! i2
      enddo ! i1
!!
      return
      end
