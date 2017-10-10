!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!! - main.f90
!!
      subroutine getsymmetryfunctions(iseed,numtrain,numtest,numrej)
!!
      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions
      use mode1options
      use symfunctions
      use nnshort_atomic
      use nnshort_pair
      use nnewald
      use nnham
      use structures
      use timings
!!
      implicit none
!!
      integer npoints                                                      ! internal
      integer ncount                                                       ! internal
      integer iseed                                                        ! in/out
      integer numtrain                                                     ! out 
      integer numtest                                                      ! out 
      integer numrej                                                       ! out 
      integer pointnumber                                                  ! internal
      integer num_atoms_element_list(nblock,nelem)                         ! internal 
      integer i1,i2                                                        ! internal
      integer ndone                                                        ! internal

      real*8 chargemin(nelem)                                              ! internal 
      real*8 chargemax(nelem)                                              ! internal 
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)                   ! internal 
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)                   ! internal 
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)                    ! internal 
      real*8 minvalue_ham(nelem,maxnum_funcvalues_ham)                     ! internal 
      real*8 maxvalue_ham(nelem,maxnum_funcvalues_ham)                     ! internal 
      real*8 avvalue_ham(nelem,maxnum_funcvalues_ham)                      ! internal
      real*8 avvalue_s(nelem,maxnum_funcvalues_s)                          ! internal 
      real*8 minvalue_s(nelem,maxnum_funcvalues_s)                         ! internal
      real*8 maxvalue_s(nelem,maxnum_funcvalues_s)                         ! internal
      real*8 avvalue_hexton(nelem,maxnum_funcvalues_hexton)                ! internal
      real*8 minvalue_hexton(nelem,maxnum_funcvalues_hexton)               ! internal
      real*8 maxvalue_hexton(nelem,maxnum_funcvalues_hexton)               ! internal
      real*8 avvalue_hextoff(1,maxnum_funcvalues_hextoff)                  ! internal
      real*8 minvalue_hextoff(1,maxnum_funcvalues_hextoff)                 ! internal
      real*8 maxvalue_hextoff(1,maxnum_funcvalues_hextoff)                 ! internal
      real*8 avvalue_dens(nelem,maxnum_funcvalues_dens)                    ! internal
      real*8 minvalue_dens(nelem,maxnum_funcvalues_dens)                   ! internal
      real*8 maxvalue_dens(nelem,maxnum_funcvalues_dens)                   ! internal

      real*8 dummy                                                         ! internal
!!
!!===================================================================================
!! initializations 
!!===================================================================================
      pointnumber=0
!! 
!!===================================================================================
!! get NN data for charges if electrostatics shall be removed from reference E and F
!!===================================================================================
      if(lshort.and.lelec.and.(nn_type_elec.eq.1))then
        weights_elec(:,:)       = 0.0d0
        minvalue_elec(:,:)      = 0.0d0
        maxvalue_elec(:,:)      = 0.0d0
        avvalue_elec(:,:)       = 0.0d0
        write(ounit,*)'Reading charge fit data for determination of electrostatic forces'
        write(ounit,*)'-------------------------------------------------------------'
!!'
!!===================================================================================
!! read scalinge data for electrostatics
!!===================================================================================
        call readscale(nelem,3,&
          maxnum_funcvalues_elec,num_funcvalues_elec,&
          minvalue_elec,maxvalue_elec,avvalue_elec,&
          dummy,dummy,chargemin,chargemax)
!!===================================================================================
!! read electrostatic weights
!!===================================================================================
        call readweights(1,nelem,&
          maxnum_weights_elec,num_weights_elec,&
          weights_elec)
      endif ! lelec
!!
!!===================================================================================
!! get NN data for Hamiltonian if NNTB E and F shall be removed from reference data 
!!===================================================================================
      if(lshort.and.lnntb)then
        write(ounit,*)'ERROR in getsymmetryfunctions, implement reading Hamiltonian'
        stop
      endif ! lnntb'
!!
!!===================================================================================
!!===================================================================================
!! calculate the symmetry functions for all structures in blocks of structures
!!===================================================================================
!!===================================================================================
      ncount=totnum_structures
      ndone=0
!! process next group of points
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!      write(*,*) ncount, npoints
      
!!
!!===================================================================================
!! read data file(s) into the arrays of module structures for a block of npoints structures
!!===================================================================================
      call readstructures(npoints,num_atoms_element_list)
      
!!
!!===================================================================================
!! remove atomic energies from total energies if requested
!!===================================================================================
      if(lremoveatomenergies)then
        call removeatoms(nblock,npoints, &
          num_atoms_list,zelem_list,&
          num_atoms_element_list,&
          totalenergy_list,atomenergy_list)
      endif
!!
!!===================================================================================
!! if charges are fixed, overwrite DFT reference charges by charges from input.nn
!! this is needed to subtract the correct electrostatic energies from the total energy
!!===================================================================================
      if(nn_type_elec.eq.3)then
        do i1=1,npoints
          do i2=1,num_atoms_list(i1)
            atomcharge_list(i1,i2)=fixedcharge(elementindex(zelem_list(i1,i2)))
          enddo
        enddo
      endif
      
!!
!!===================================================================================
!! calculate and write symmetry functions for a block of structures
!!===================================================================================
      call calcfunctions(npoints,ndone,&
        iseed,numtrain,numtest,numrej,pointnumber,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        minvalue_ham,maxvalue_ham,avvalue_ham,&
        minvalue_s,maxvalue_s,avvalue_s,&
        minvalue_hexton,maxvalue_hexton,avvalue_hexton,&
        minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
        minvalue_dens,maxvalue_dens,avvalue_dens)
!!
      ndone=ndone+npoints
!!===================================================================================
!!===================================================================================
!! if there are structures left go to next group of structures
      if(ncount.gt.0) goto 10
!!===================================================================================
!!===================================================================================
!!
      return
      end
