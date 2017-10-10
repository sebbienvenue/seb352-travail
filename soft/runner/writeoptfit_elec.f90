!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine writeoptfit_elec(optepoche,tounit,&
        optrmse_ewald,optrmse_ewald_test,&
        optrmse_charge,optrmse_charge_test,optrmse_totalcharge,&
        optrmse_totalcharge_test)
!!
      use fileunits
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer optepoche                ! in
!!
      real*8 tounit                    ! in
      real*8 optrmse_ewald             ! in
      real*8 optrmse_ewald_test        ! in
      real*8 optrmse_charge            ! in
      real*8 optrmse_charge_test       ! in
      real*8 optrmse_totalcharge       ! in
      real*8 optrmse_totalcharge_test  ! in
!!
      write(ounit,*)'============================================================='
      if(optepoche.eq.0)then
        write(ounit,'(a,i5)')' No improvement of the charge fit has been obtained'
      else
        write(ounit,'(a,i5)')' Best charge fit has been obtained in epoch ',optepoche
        write(ounit,'(3a)')'                    --- E_elec: --- ',&
                     '          --- charges: ---',&
                     '      --- total charges: ---'
        write(ounit,'(3a)')'                   train         test ',&
                     '       train         test',&
                     '        train         test'
        write(ounit,'(a11,x,6f13.6)')' OPTCHARGE ',&
        optrmse_ewald*tounit,optrmse_ewald_test*tounit,&
        optrmse_charge,optrmse_charge_test,&
        optrmse_totalcharge,optrmse_totalcharge_test
      endif ! optepoche.eq.0
      write(ounit,*)'-------------------------------------------------------------'
!! write maximum errors
      if(fitting_unit.eq.1)then
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (train set): ',&
          maxerror_elec_train*tounit,' eV/atom (structure ',imaxerror_elec_train,')' 
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (test set) : ',&
          maxerror_elec_test*tounit, ' eV/atom (structure ',imaxerror_elec_test,')' 
      elseif(fitting_unit.eq.2)then
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (train set): ',&
          maxerror_elec_train*tounit,' Ha/atom (structure ',imaxerror_elec_train,')' 
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (test set) : ',&
          maxerror_elec_test*tounit, ' Ha/atom (structure ',imaxerror_elec_test,')' 
      endif
!!
      return
      end
