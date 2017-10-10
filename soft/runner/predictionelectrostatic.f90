!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - predict.f90
!!
      subroutine predictionelectrostatic(&
        num_atoms,zelem,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        lattice,xyzstruct,&
        nntotalcharge,nnatomcharge,&
        chargemin,chargemax,nnelecenergy,&
        nnelecforce,nnstress_elec,sense,lperiodic)
!!
      use globaloptions
      use fileunits
      use predictionoptions
      use nnewald
      use timings
!!
      implicit none
!!
      integer zelem(max_num_atoms)                                      ! in
      integer num_atoms                                                 ! in
!!
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)                ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)                ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)                 ! in
      real*8 nnelecenergy                                               ! out, total electrostatic energy 
      real*8 chargemin(nelem)                                           ! in
      real*8 chargemax(nelem)                                           ! in
      real*8 nntotalcharge                                              ! out 
      real*8 nnatomcharge(max_num_atoms)                                ! out 
      real*8 nnstress_elec(3,3)                                         ! out 
      real*8 sense(nelem,maxnum_funcvalues_elec)                        ! out  
      real*8 nnelecforce(3,max_num_atoms)                               ! out 
      real*8 lattice(3,3)                                               ! in
      real*8 xyzstruct(3,max_num_atoms)                                 ! in
!! 
      logical lperiodic                                                 ! in

!!
!!======================================================================
!! Initialzations
!!======================================================================
!! initializations electrostatic part
      nntotalcharge          = 0.0d0
      nnatomcharge(:)        = 0.0d0
      nnelecenergy           = 0.0d0
      nnelecforce(:,:)       = 0.0d0
      nnstress_elec(:,:)     = 0.0d0
!! initialization of timings
      timeelec               = 0.0d0
!! initialization of sensitivities
      sense(:,:)             = 0.0d0
!!
!!======================================================================
!!======================================================================
!! Start electrostatic part
!!======================================================================
!!======================================================================
      if(lfinetime)then
        dayelec=0
        call abstime(timeelecstart,dayelec)
      endif ! lfinetime
!!
      call predictelec(&
        num_atoms,zelem,lattice,xyzstruct,&
        minvalue_elec,maxvalue_elec,avvalue_elec,chargemin,chargemax,&
        nnelecenergy,nnatomcharge,nntotalcharge,nnelecforce,&
        nnstress_elec,sense,lperiodic)
!!
      if(lfinetime)then
        call abstime(timeelecend,dayelec)
        timeelec=timeelec+timeelecend-timeelecstart
      endif ! lfinetime
!!======================================================================
!!======================================================================
!! End Electrostatic part 
!!======================================================================
!!======================================================================
      return
      end
