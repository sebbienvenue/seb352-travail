!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readstructures.f90 
!! - predict.f90
!!
!! It is assumed that all inconsistencies in the input.data file
!! have been detected before by checkonestructure
!!
      subroutine readonestructure(num_atoms,&
      zelem,num_atoms_element,lattice,&
      totalcharge,totalenergy,&
      atomcharge,atomenergy,xyzstruct,&
      totalforce,elementsymbol,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer nlattice                         ! internal 
      integer num_atoms                        ! out 
      integer zelem(max_num_atoms)             ! out 
      integer num_atoms_element(nelem)         ! out
      integer i                                ! internal 
!!
      real*8 lattice(3,3)                      ! out 
      real*8 xyzstruct(3,max_num_atoms)        ! out
      real*8 totalforce(3,max_num_atoms)       ! out
      real*8 atomcharge(max_num_atoms)         ! out
      real*8 atomenergy(max_num_atoms)         ! out
      real*8 totalcharge                       ! out
      real*8 totalenergy                       ! out 
      real*8 threshold                         ! internal 
      real*8 atomicchargesum                   ! internal 
      real*8 atomicenergysum                   ! internal 
!!
      character*20 keyword                     ! internal
      character*2 elementsymbol(max_num_atoms) ! out
!!
      logical lbegin                           ! internal, true if 'begin' flag of structure found
      logical lend                             ! internal, true if 'end' flag of structure found
      logical lcharge                          ! internal, true if 'charge' keyword is found
      logical lenergy                          ! internal, true if 'energy' keyword is found
      logical lperiodic                        ! out, true if 3 lattice vectors are found 
!!
!!=====================================================
!! initializations
!!=====================================================
      nlattice             = 0
      lbegin               = .false.
      lend                 = .false.
      lcharge              = .false.
      lenergy              = .false.
      lperiodic            = .false.
      num_atoms            = 0
      threshold            = 0.000001d0
      atomicenergysum      = 0.0d0
      atomicchargesum      = 0.0d0
      lattice(:,:)         = 0.0d0
      zelem(:)             = 0
      xyzstruct(:,:)       = 0.0d0
      totalforce(:,:)      = 0.0d0
      elementsymbol(:)     = '  '
      num_atoms_element(:) = 0
      totalcharge          = 0.0d0
!!
 10   continue
      read(dataunit,*)keyword
!! 
      if(keyword.eq.'begin') lbegin=.true. 
      if(keyword.eq.'lattice') then
        nlattice=nlattice+1
        backspace(dataunit)
        read(dataunit,*)keyword,(lattice(nlattice,i),i=1,3)
      endif
      if(keyword.eq.'atom') then
        num_atoms=num_atoms+1
        backspace(dataunit)
        read(dataunit,*)keyword,(xyzstruct(i,num_atoms),i=1,3),&
            elementsymbol(num_atoms),atomcharge(num_atoms),&
            atomenergy(num_atoms),(totalforce(i,num_atoms),i=1,3)

        call nuccharge(elementsymbol(num_atoms),zelem(num_atoms))
      endif
      if(keyword.eq.'charge') then
        backspace(dataunit)
        if(lcharge) then
          write(*,*)'Error: 2 total charges given for structure'
          stop
        endif 
        read(dataunit,*) keyword,totalcharge
        lcharge=.true.
      endif
      if(keyword.eq.'energy') then
        backspace(dataunit)
        if(lenergy) then
          write(*,*)'Error: 2 total energies given for structure'
          stop
        endif
        read(dataunit,*) keyword,totalenergy
        lenergy=.true.
      endif
      if((keyword.eq.'end').and.(lbegin)) lend=.true. 
      if(keyword.ne.'end') goto 10

      if(nlattice.eq.3) then
        lperiodic=.true.
      endif
!!
!!=====================================================
!! translate atoms back into cell if necessary
!!=====================================================
      if(lperiodic)then
        call translate(num_atoms,lattice,xyzstruct)
      endif
!!
!!=====================================================
!! count the number of atoms of each element
!!=====================================================
      do i=1,num_atoms
        num_atoms_element(elementindex(zelem(i)))=&
           num_atoms_element(elementindex(zelem(i)))+1
      enddo
!!
      return
      end
