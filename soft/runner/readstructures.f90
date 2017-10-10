!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getsymmetryfunctions.f90
!! - getpairsymfunctions.f90
!!
      subroutine readstructures(npoints,num_atoms_element_list)
!!
      use fileunits
      use globaloptions
      use structures
      use nnflags
      use basismod
!!
      implicit none
!!
      integer npoints                                 ! in
      integer num_atoms                               ! internal
      integer zelem(max_num_atoms)                    ! internal
      integer num_atoms_element(nelem)                ! internal 
      integer num_atoms_element_list(nblock,nelem)    ! out
      integer i1                                      ! internal
!!
      real*8 totalcharge                              ! internal
      real*8 totalenergy                              ! internal
      real*8 atomcharge(max_num_atoms)                ! internal
      real*8 atomenergy(max_num_atoms)                ! internal
      real*8 xyzstruct(3,max_num_atoms)               ! internal
      real*8 totalforce(3,max_num_atoms)              ! internal
      real*8 hextoff_value(maxnum_basis,maxnum_basis) ! internal
!!      real*8 hextoff_list(nblock,maxnum_basis,maxnum_basis) ! internal NEEDS TO GO TO A MODULE LATER
!!
      character*2 elementsymbol(max_num_atoms)        ! internal
!!
      logical lperiodic                               ! internal
!!
!! initialization
      num_atoms_element_list(:,:)=0
      originatom_id = 0
      zatom_id = 0
!!
      do i1=1,npoints
        num_atoms = 0
!!
!!      write(*,*) 'Structure', i1, ' read'
        if(lnntb) then
          call readonestructure_ham(num_atoms,&
            zelem,num_atoms_element,lattice_list(1,1,i1),&
            totalcharge,totalenergy,&
            atomcharge,atomenergy,xyzstruct,&
            totalforce,hextoff_value,elementsymbol,lperiodic,i1)
        else
          call readonestructure(num_atoms,&
            zelem,num_atoms_element,lattice_list(1,1,i1),&
            totalcharge,totalenergy,&
            atomcharge,atomenergy,xyzstruct,&
            totalforce,elementsymbol,lperiodic,i1)
        endif
!!
        num_atoms_list(i1)          =num_atoms
        zelem_list(i1,:)            =zelem(:)
        totalcharge_list(i1)        =totalcharge
        totalenergy_list(i1)        =totalenergy
        atomcharge_list(i1,:)       =atomcharge(:)
        atomenergy_list(i1,:)       =atomenergy(:)
        xyzstruct_list(:,:,i1)      =xyzstruct(:,:)
        totalforce_list(:,:,i1)     =totalforce(:,:)
        elementsymbol_list(i1,:)    =elementsymbol(:)
        lperiodic_list(i1)          =lperiodic
        num_atoms_element_list(i1,:)=num_atoms_element(:)
        if(lnntb) then
          hextoff_list(i1,:,:)      =hextoff_value(:,:)
        endif
!!
      enddo ! i1
!!
      
      return
      end
