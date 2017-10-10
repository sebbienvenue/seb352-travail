!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine adjustkalman_elec(numberq,kalmanthresholde_temp)
!!
      use fileunits
      use nnflags 
      use globaloptions
      use fittingoptions
      use mpi_mod
!!
      implicit none
!!
      integer numberq                              ! in
!!
      real*8 kalmanthresholde_temp                 ! in/out
      real*8 ztemp                                 ! internal
!!
!! check if a non-zero number of charges has been used for the update
      if((numberq.eq.0).and.(kalmanthresholde.gt.0.0d0).and.(chargernd.gt.0.0d0))then
        kalmanthresholde_temp=kalmanthresholde_temp*0.9d0
        if(mpirank.eq.0)then
          write(ounit,*)'### WARNING ### kalmanthresholde has been adjusted to ',kalmanthresholde_temp
        endif
      elseif((numberq.gt.0).and.lelec.and.(nn_type_elec.eq.1).and.(kalmanthresholde.gt.0.0d0))then
        ztemp=kalmanthresholde_temp
        kalmanthresholde_temp=kalmanthresholde_temp/0.9d0
        kalmanthresholde_temp=min(kalmanthresholde_temp,kalmanthresholde)
        if(abs(ztemp-kalmanthresholde_temp).gt.0.000000001d0)then
          if(mpirank.eq.0)then
            write(ounit,*)'### WARNING ### kalmanthresholde has been adjusted to ',kalmanthresholde_temp
          endif
        endif
      endif
!!
      return
      end
