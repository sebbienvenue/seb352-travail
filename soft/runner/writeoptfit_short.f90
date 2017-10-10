!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine writeoptfit_short(optepoch,&
        tounit,optrmse_short,optrmse_short_test,optrmse_force_s,&
        optrmse_force_s_test)
!!
      use fileunits
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer optepoch                 ! in
!!
      real*8 tounit                    ! in
      real*8 optrmse_short             ! in
      real*8 optrmse_short_test        ! in
      real*8 optrmse_force_s           ! in
      real*8 optrmse_force_s_test      ! in
!!
      write(ounit,*)'============================================================='
        if(optepoch.eq.0)then
          write(ounit,'(a,i5)')' No improvement of the short range fit has been obtained'
        else
          write(ounit,'(a,i5)')' Best short range fit has been obtained in epoch ',optepoch
          if((luseforces.and.(.not.lfinalforce))&
            .or.(lfinalforce.and.(optepoch.eq.nepochs)))then
            write(ounit,'(2a)')'                    --- E_short: --- ',&
                           '         --- F_short: ---'
            write(ounit,'(2a)')'                   train         test ',&
                           '       train         test'
            write(ounit,'(a11,x,4f13.6)')' OPTSHORT  ',&
              optrmse_short*tounit,optrmse_short_test*tounit,&
              optrmse_force_s*tounit,optrmse_force_s_test*tounit
          else ! no forces:
            write(ounit,'(a)')'                    --- E_short: --- '
            write(ounit,'(a)')'                   train         test '
            write(ounit,'(a11,x,2f13.6)')' OPTSHORT  ',&
              optrmse_short*tounit,optrmse_short_test*tounit
          endif ! luseforces
        endif ! optepoch.eq.0
        write(ounit,*)'-------------------------------------------------------------'
!! write maximum errors
        if(fitting_unit.eq.1)then
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (train set): ',&
            maxerror_eshort_train*tounit,' eV/atom (structure ',imaxerror_eshort_train,')' 
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (test set) : ',&
            maxerror_eshort_test*tounit, ' eV/atom (structure ',imaxerror_eshort_test,')' 
        elseif(fitting_unit.eq.2)then
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (train set): ',&
            maxerror_eshort_train*tounit,' Ha/atom (structure ',imaxerror_eshort_train,')' 
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (test set) : ',&
            maxerror_eshort_test*tounit, ' Ha/atom (structure ',imaxerror_eshort_test,')' 
        endif
!!
      return
      end
