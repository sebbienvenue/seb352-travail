!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!! - geterror.f90 
!! - geterrorpair.f90
!!
      subroutine getrmse(ntrain,ncharges,nenergies,netot,nforces,nforcese,&
        nforcest,&
        rmse_short,rmse_charge,rmse_totalcharge,rmse_ewald,rmse_etot,&
        rmse_force_s,rmse_force_e,rmse_force_t,&
        mad_short,mad_charge,mad_totalcharge,mad_ewald,mad_etot,&
        mad_force_s,mad_force_e,mad_force_t)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer ntrain,ncharges,nenergies,nforces,nforcese,nforcest,netot  ! in
!!
      real*8 rmse_short                         ! in/out
      real*8 rmse_ewald                         ! in/out
      real*8 rmse_charge                        ! in/out
      real*8 rmse_totalcharge                   ! in/out
      real*8 rmse_etot                          ! in/out
      real*8 rmse_force_s                       ! in/out
      real*8 rmse_force_e                       ! in/out
      real*8 rmse_force_t                       ! in/out
      real*8 mad_short                          ! in/out
      real*8 mad_ewald                          ! in/out
      real*8 mad_charge                         ! in/out
      real*8 mad_totalcharge                    ! in/out
      real*8 mad_etot                           ! in/out
      real*8 mad_force_s                        ! in/out
      real*8 mad_force_e                        ! in/out
      real*8 mad_force_t                        ! in/out
!!
!!
!! some security checks if settings for maxenergy or maxforce are inappropriate:
      if(lshort)then
        if(nenergies.eq.0)then
          write(ounit,*)'ERROR in getrmse: nenergies is zero '
          write(ounit,*)'Probably maxenergy and/or energy_threshold is set too low'
          stop  !'
        endif
      endif
!!
!! short range
      if(lshort)then
        rmse_short          =rmse_short/dble(nenergies)
        rmse_short          =dsqrt(rmse_short)
!        write(ounit,*)'getrmse nforces ',nforces
        rmse_force_s        =rmse_force_s/dble(nforces)
        rmse_force_s        =sqrt(rmse_force_s)
        mad_short           =mad_short/dble(nenergies)
        mad_force_s         =mad_force_s/dble(nforces)
      endif
!! electrostatic
      if(lelec)then
        rmse_charge         =rmse_charge/dble(ncharges)
        rmse_charge         =dsqrt(rmse_charge)
        rmse_totalcharge    =rmse_totalcharge/dble(ntrain)
        rmse_totalcharge    =dsqrt(rmse_totalcharge)
        rmse_ewald          =rmse_ewald/dble(ntrain)
        rmse_ewald          =dsqrt(rmse_ewald)
        rmse_force_e        =rmse_force_e/dble(nforcese)
        rmse_force_e        =sqrt(rmse_force_e)
        mad_charge          =mad_charge/dble(ncharges)
        mad_totalcharge     =mad_totalcharge/dble(ntrain)
        mad_ewald           =mad_ewald/dble(ntrain)
        mad_force_e         =mad_force_e/dble(nforcese)
      endif
!! total
      rmse_etot             =rmse_etot/dble(netot)
      rmse_etot             =dsqrt(rmse_etot)
      rmse_force_t          =rmse_force_t/dble(nforcest)
      rmse_force_t          =sqrt(rmse_force_t)
      mad_etot              =mad_etot/dble(netot)
      mad_force_t           =mad_force_t/dble(nforcest)
!!
      return
      end

