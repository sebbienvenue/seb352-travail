!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - calcfunctions.f90
!! - prediction.f90 
!! - getallelectrostatic.f90
!!
      subroutine getewaldenergy(max_num_neighbors_elec,&
        neighboridx_elec,num_neighbors_elec,num_atoms,zelem,&
        lattice,xyzstruct,atomcharge,elecenergy,&
        dchargedxyz,nnelecforce,ldoforces_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer max_num_neighbors_elec       ! in
      integer num_neighbors_elec(num_atoms)! in
      integer neighboridx_elec(num_atoms,0:max_num_neighbors_elec)
      integer num_atoms                    ! in
      integer zelem(max_num_atoms)         ! in
      integer day                          ! internal
!!
      real*8 lattice(3,3)                  ! in
      real*8 xyzstruct(3,max_num_atoms)    ! in
      real*8 dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)    ! in  
      real*8 atomcharge(max_num_atoms)     ! in
      real*8 elecenergy                    ! out
      real*8 eself                         ! internal
      real*8 ereal                         ! internal
      real*8 erecip                        ! internal
      real*8 timestart                     ! internal
      real*8 timeend                       ! internal
      real*8 nnelecforce(3,max_num_atoms)  ! out 
      real*8 erealforce(3,max_num_atoms)   ! internal 
      real*8 erecipforce(3,max_num_atoms)  ! internal 
      real*8 eselfforce(3,max_num_atoms)   ! internal 
      real*8 ewaldcorr                     ! internal
      real*8 ewaldcorrforce(3,max_num_atoms) ! internal
!!
      logical ldoforces_local                    ! in
!!
!! initializations
      elecenergy      =0.0d0
      eself           =0.0d0
      ereal           =0.0d0
      erecip          =0.0d0
      nnelecforce(:,:)=0.0d0
      erealforce(:,:)  =0.0d0
      erecipforce(:,:) =0.0d0
      eselfforce(:,:)  =0.0d0
      ewaldcorr        =0.0d0
      ewaldcorrforce(:,:)=0.0d0
!!
!! calculate the real space energy part ereal
      if(lfinetime)then
        timestart=0.0d0
        day=0
        call abstime(timestart,day)
      endif
!!
      call ewaldreal(max_num_neighbors_elec,num_neighbors_elec,&
        neighboridx_elec,num_atoms,zelem,&
        atomcharge,ereal,dchargedxyz,&
        lattice,xyzstruct,erealforce,&
        .true.,ldoforces_local)
!!
      if(lfinetime)then
        call abstime(timeend,day)
        write(ounit,'(a,f14.3)')&
        'time calculation Ewald real ',timeend-timestart
      endif !'
!!
!! calculate the reciprocal space energy part erecip
      if(lfinetime)then
        timestart=0.0d0
        day=0
        call abstime(timestart,day)
      endif
!!
      call ewaldrecip(max_num_neighbors_elec,num_atoms,xyzstruct,&
        atomcharge,erecipforce,erecip,lattice,dchargedxyz,ldoforces_local)
      if(lfinetime)then
        call abstime(timeend,day)
        write(ounit,'(a,f14.3)')&
        'time calculation Ewald recip ',timeend-timestart
      endif !'
!!
!! calculate the self-term energy part eself
      if(lfinetime)then
        timestart=0.0d0
        day=0
        call abstime(timestart,day)
      endif
!!
      call ewaldself(max_num_neighbors_elec,neighboridx_elec,&
        num_neighbors_elec,num_atoms,&
        atomcharge,dchargedxyz,eselfforce,eself,&
        ldoforces_local)
!!
      if(lfinetime)then
        call abstime(timeend,day)
        write(ounit,'(a,f14.3)')&
        'time calculation Ewald self ',timeend-timestart
      endif !'
!!
      if(lscreen) then
        call getewaldcorr(max_num_neighbors_elec,neighboridx_elec,num_atoms,zelem,&
          atomcharge,dchargedxyz,lattice,xyzstruct,&
          .true.,ldoforces_local,ewaldcorr,ewaldcorrforce)
!!
        elecenergy = ereal+erecip+eself-ewaldcorr
        nnelecforce(:,:)=erealforce(:,:)+erecipforce(:,:)+eselfforce(:,:)-ewaldcorrforce(:,:)
      else
!! sum the ewald energy contributions
!! elecenergy is the total energy here, not per atom!
        elecenergy = ereal+erecip+eself  ! note sign of eself!
!! sum the ewald force contributions
        nnelecforce(:,:)=erealforce(:,:)+erecipforce(:,:)+eselfforce(:,:)
      endif
!!
      return
      end 
