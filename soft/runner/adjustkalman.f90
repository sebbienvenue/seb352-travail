!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine adjustkalman(numbere,numberf,numberq,&
        kalmanthreshold_temp,kalmanthresholdf_temp,kalmanthresholde_temp)
!!
      use fileunits
      use nnflags 
      use globaloptions
      use fittingoptions
      use mpi_mod
!!
      implicit none
!!
      integer numbere                              ! in
      integer numberf                              ! in
      integer numberq                              ! in
!!
      real*8 kalmanthreshold_temp                  ! in/out
      real*8 kalmanthresholdf_temp                 ! in/out
      real*8 kalmanthresholde_temp                 ! in/out
      real*8 ztemp                                 ! internal
!!
!! check if a non-zero number of structures has been used for the energy update
      if((numbere.eq.0).and.(lshort).and.(kalmanthreshold.gt.0.0d0).and.(energyrnd.gt.0.0d0))then
        kalmanthreshold_temp=kalmanthreshold_temp*0.9d0
        if(mpirank.eq.0)then
          write(ounit,*)'### WARNING ### kalmanthreshold has been adjusted to ',kalmanthreshold_temp
        endif
      elseif((numbere.gt.0).and.(lshort).and.(kalmanthreshold.gt.0.0d0))then
        ztemp=kalmanthreshold_temp
        kalmanthreshold_temp=kalmanthreshold_temp/0.9d0
        kalmanthreshold_temp=min(kalmanthreshold_temp,kalmanthreshold)
        if(abs(ztemp-kalmanthreshold_temp).gt.0.000000001d0)then
          if(mpirank.eq.0)then
            write(ounit,*)'### WARNING ### kalmanthreshold has been adjusted to ',kalmanthreshold_temp
          endif
        endif
      endif
!! check if a non-zero number of structures has been used for the force update
      if((numberf.eq.0).and.(lshort).and.(luseforces).and.(kalmanthresholdf.gt.0.0d0).and.(forcernd.gt.0.0d0))then
        kalmanthresholdf_temp=kalmanthresholdf_temp*0.9d0
        if(mpirank.eq.0)then
          write(ounit,*)'### WARNING ### kalmanthresholdf has been adjusted to ',kalmanthresholdf_temp
        endif
      elseif((numberf.gt.0).and.(lshort).and.(luseforces).and.(kalmanthresholdf.gt.0.0d0))then
        ztemp=kalmanthresholdf_temp
        kalmanthresholdf_temp=kalmanthresholdf_temp/0.9d0
        kalmanthresholdf_temp=min(kalmanthresholdf_temp,kalmanthresholdf)
        if(abs(ztemp-kalmanthresholdf_temp).gt.0.000000001d0)then
          if(mpirank.eq.0)then
            write(ounit,*)'### WARNING ### kalmanthresholdf has been adjusted to ',kalmanthresholdf_temp
          endif
        endif
      endif
!! check if a non-zero number of charges has been used for the update
      if((numberq.eq.0).and.lelec.and.(nn_type_elec.eq.1).and.(kalmanthresholde.gt.0.0d0).and.(chargernd.gt.0.0d0))then
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
