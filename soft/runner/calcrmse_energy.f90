!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - ewaldenergies_para.f90
!! - geterror.f90
!! - geterrorpair.f90
!! - getshortenergies_para.f90
!! - getshortenergies_parapair.f90
!!
      subroutine calcrmse_energy(ndim,npoints,nenergies,&
          ndone,imaxerror_local,&
          rmse_local,mad_local,maxerror_local,&
          maxenergy_local,energy_local,nnenergy_local)
!!
      use fileunits
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer ndim                        ! in
      integer ndone                       ! in
      integer npoints                     ! in
      integer i1                          ! internal
      integer nenergies                   ! in/out
      integer imaxerror_local             ! in/out
!!
      real*8 rmse_local                   ! in/out
      real*8 mad_local                    ! in/out
      real*8 energy_local(ndim)           ! in
      real*8 nnenergy_local(ndim)         ! in
      real*8 maxenergy_local              ! in
      real*8 maxerror_local               ! in/out
!!
!!
      do i1=1,npoints
        if(energy_local(i1).le.maxenergy_local)then
          nenergies = nenergies+1
          rmse_local= rmse_local +(energy_local(i1)-nnenergy_local(i1))**2
          mad_local = mad_local  +abs(energy_local(i1)-nnenergy_local(i1))
          if(abs(energy_local(i1)-nnenergy_local(i1)).gt.maxerror_local)then
            maxerror_local = abs(energy_local(i1)-nnenergy_local(i1)) 
            imaxerror_local= ndone + i1
          endif
        endif
      enddo
!!
      return
      end
