!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine writeoptfit(optepoch,optepoche,&
        tounit,optrmse_short,optrmse_short_test,optrmse_force_s,&
        optrmse_force_s_test,optrmse_ewald,optrmse_ewald_test,&
        optrmse_charge,optrmse_charge_test,optrmse_totalcharge,&
        optrmse_totalcharge_test,optrmse_force_e,optrmse_force_e_test,&
        optrmse_etot,optrmse_etot_test)
!!
      use fileunits
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer optepoch                 ! in
      integer optepoche                ! in
!!
      real*8 tounit                    ! in
      real*8 optrmse_short             ! in
      real*8 optrmse_short_test        ! in
      real*8 optrmse_force_s           ! in
      real*8 optrmse_force_s_test      ! in
      real*8 optrmse_ewald             ! in
      real*8 optrmse_ewald_test        ! in
      real*8 optrmse_charge            ! in
      real*8 optrmse_charge_test       ! in
      real*8 optrmse_totalcharge       ! in
      real*8 optrmse_totalcharge_test  ! in
      real*8 optrmse_force_e           ! in
      real*8 optrmse_force_e_test      ! in
      real*8 optrmse_etot              ! in
      real*8 optrmse_etot_test         ! in
!!
      write(ounit,*)'============================================================='
      if(lshort)then
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
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (train set): ',maxerroreshorttrain*tounit,' eV/atom (structure ',imaxerror_eshort_train,')' 
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (test set) : ',maxerroreshorttest*tounit, ' eV/atom (structure ',imaxerror_eshort_test,')' 
        elseif(fitting_unit.eq.2)then
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (train set): ',maxerroreshorttrain*tounit,' Ha/atom (structure ',imaxerror_eshort_train,')' 
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eshort error in last epoch (test set) : ',maxerroreshorttest*tounit, ' Ha/atom (structure ',imaxerror_eshort_test,')' 
        endif
      endif ! lshort
      if(lshort.and.lelec)then
        write(ounit,*)'-------------------------------------------------------------'
      endif
      if(lelec)then
        if(optepoche.eq.0)then
          write(ounit,'(a,i5)')' No improvement of the charge fit has been obtained'
        else
          write(ounit,'(a,i5)')' Best charge fit has been obtained in epoch ',optepoche
          if(luseforces)then
            write(ounit,'(4a)')'                    --- E_elec: --- ',&
                         '          --- charges: ---',&
                         '      --- total charges: ---',&
                         '        --- F_elec: ---'
            write(ounit,'(4a)')'                   train         test ',&
                         '       train         test',&
                         '        train         test',&
                         '        train         test'
            write(ounit,'(a11,x,8f13.6)')' OPTCHARGE ',&
            optrmse_ewald*tounit,optrmse_ewald_test*tounit,&
            optrmse_charge,optrmse_charge_test,&
            optrmse_totalcharge,optrmse_totalcharge_test,&
            optrmse_force_e*tounit,optrmse_force_e_test*tounit
          else ! no forces:
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
          endif ! luseforces
        endif ! optepoche.eq.0
        write(ounit,*)'-------------------------------------------------------------'
!! write maximum errors
        if(fitting_unit.eq.1)then
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (train set): ',maxerrorewaldtrain*tounit,' eV/atom (structure ',imaxerror_elec_train,')' 
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (test set) : ',maxerrorewaldtest*tounit, ' eV/atom (structure ',imaxerror_elec_test,')' 
        elseif(fitting_unit.eq.2)then
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (train set): ',maxerrorewaldtrain*tounit,' Ha/atom (structure ',imaxerror_elec_train,')' 
          write(ounit,'(a,f14.6,a21,i8,a2)')' max Eelec error in last epoch (test set) : ',maxerrorewaldtest*tounit, ' Ha/atom (structure ',imaxerror_elec_test,')' 
        endif
      endif ! lelec
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,'(a)')'                      --- E_tot: --- '
      write(ounit,'(a)')'                   train         test '
      write(ounit,'(a11,x,2f13.6)')' OPTETOT   ',&
        optrmse_etot*tounit,optrmse_etot_test*tounit
      write(ounit,*)'-------------------------------------------------------------'
!! write maximum errors
      if(fitting_unit.eq.1)then
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Etot error in last epoch (train set): ',maxerroretottrain*tounit,' eV/atom (structure ',imaxerror_etot_train,')' 
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Etot error in last epoch (test set) : ',maxerroretottest*tounit, ' eV/atom (structure ',imaxerror_etot_test,')' 
      elseif(fitting_unit.eq.2)then
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Etot error in last epoch (train set): ',maxerroretottrain*tounit,' Ha/atom (structure ',imaxerror_etot_train,')' 
        write(ounit,'(a,f14.6,a21,i8,a2)')' max Etot error in last epoch (test set) : ',maxerroretottest*tounit, ' Ha/atom (structure ',imaxerror_etot_test,')' 
      endif
!!
      return
      end
