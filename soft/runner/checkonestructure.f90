!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - checkstructures.f90
!!
      subroutine checkonestructure(inumber,lelement)
!!
      use fileunits
      use globaloptions
      use basismod
      use nnflags
!!
      implicit none
!!
      integer nlattice ! number of lattice vectors found
      integer num_atoms ! number of atoms in structure
      integer zelem(max_num_atoms) ! nuclear charges of the atoms
      integer inumber
      integer i ! internal 
!!
      real*8 lattice(3,3) ! lattice vectors if present
      real*8 xyzstruct(3,max_num_atoms)
      real*8 totalforce(3,max_num_atoms)
      real*8 atomcharge(max_num_atoms)
      real*8 atomenergy(max_num_atoms)
      real*8 totalcharge  ! total charge of structure
      real*8 totalenergy  ! total energy of structure
      integer ham_struct(6)  ! hamiltonian sub-block structure
      real*8 threshold          ! internal threshold for consistency checks of energy and charge
      real*8 atomicchargesum    ! temp variable
      real*8 atomicenergysum    ! temp variable
      real*8 volume             ! internal
!!
      character*21 keyword
      character*2 elementsymbol(max_num_atoms)
!!
      logical lbegin     ! true if 'begin' flag of structure found
      logical lend       ! true if 'end' flag of structure found
      logical lcharge    ! true if 'charge' keyword is found
      logical lenergy    ! true if 'energy' keyword is found
      logical lperiodic  ! true if 3 lattice vectors are found 
      logical lelement(102)
!!
!! initializations
      nlattice   = 0
      lbegin     = .false.
      lend       = .false.
      lcharge    = .false.
      lenergy    = .false.
      lperiodic  = .false.
      num_atoms  = 0
      threshold  = 0.005d0
      atomicenergysum = 0.0d0
      atomicchargesum = 0.0d0
      lattice(:,:)    = 0.0d0
      zelem(:) = 0
      ham_struct = 0
!!
 10   continue
      read(dataunit,*)keyword
!! CMH      write(*,*) keyword
!! 
      if(keyword.eq.'begin') lbegin=.true. 
      if(keyword.eq.'lattice') then
        nlattice=nlattice+1
        backspace(dataunit)
        read(dataunit,*,err=90)keyword,(lattice(nlattice,i),i=1,3)
      endif
      if(keyword.eq.'atom') then
        num_atoms=num_atoms+1
        backspace(dataunit)
        read(dataunit,*,err=91)keyword,(xyzstruct(i,num_atoms),i=1,3),&
            elementsymbol(num_atoms),atomcharge(num_atoms),&
            atomenergy(num_atoms),(totalforce(i,num_atoms),i=1,3)
        call nuccharge(elementsymbol(num_atoms),zelem(num_atoms))
        lelement(zelem(num_atoms))=.true. ! element found
      endif
      if(keyword.eq.'charge') then
        backspace(dataunit)
        if(lcharge)  then
          write(ounit,*)'Error: 2 total charges given for structure ',inumber
          stop
        endif
        read(dataunit,*,err=92) keyword,totalcharge
        lcharge=.true.
      endif
      if(keyword.eq.'energy')  then
        backspace(dataunit)
        if(lenergy) then
          write(ounit,*)'Error: 2 total energies given for structure ',inumber
          stop
        endif
        read(dataunit,*,err=93) keyword,totalenergy
        lenergy=.true.
      endif

      if((keyword.eq.'end').and.(lbegin)) lend=.true. 
      if(keyword.ne.'end') goto 10
!!
!! final check of structure
      if(num_atoms.eq.0) then
        write(ounit,*)'Error: No atoms found in structure ',inumber
        stop
      endif
      if(.not.lbegin) then
        write(ounit,*)'Error: No begin flag in structure ',inumber
        stop
      endif
      if(.not.lend) then
        write(ounit,*)'Error: No end flag in structure ',inumber
        stop
      endif
      if((nlattice.ne.0).and.(nlattice.ne.3))then
        write(ounit,*)'Error: Strange number of lattice vectors ',nlattice,inumber
        stop
      endif
      if(nlattice.eq.3) then
        lperiodic=.true.
      endif
      if(.not.lcharge) then
        write(ounit,*)'### WARNING ### no total charge given, assuming 0'
        totalcharge=0.0d0
      endif
      if((.not.lenergy) .and. (.not.lnntb)) then
        write(ounit,*)'Error: no total energy given for structure ',inumber
!!        stop
      endif
!! check if sum of atomic energies is equal to total energy
      if(luseatomenergies)then
      do i=1,num_atoms
        atomicenergysum=atomicenergysum + atomenergy(i)
      enddo
      if(abs(atomicenergysum-totalenergy).gt.threshold) then
        write(ounit,*)'Error: sum of atomic energies differs from total energy',inumber
        write(ounit,*)'Sum of energies:  Total energy:'
        write(ounit,'(2f20.8)')atomicenergysum,totalenergy
        stop
      endif
      endif
!!
!! check if lattice vectors make sense
      if(lperiodic)then
        call getvolume(lattice,volume)
        if(volume.lt.0.0000001d0)then
          write(ounit,*)'ERROR: volume of a periodic structure is very small ',volume
          stop
        endif
      endif
!!
!! check if sum of atomic charges is equal to total charge
      if(lelec.and.(nn_type_elec.eq.1))then
        if(luseatomcharges)then
        do i=1,num_atoms
          atomicchargesum=atomicchargesum + atomcharge(i)
        enddo
        if(abs(atomicchargesum-totalcharge).gt.threshold) then
          write(ounit,*)'### WARNING ### sum of atomic charges differs from total charge',inumber
          write(ounit,'(a,2f20.8)')'Sum of charges:  Total charge: ',atomicchargesum,totalcharge
        endif
        endif
      endif
!!
!! write all read-in information for debugging
!!      if(ldebug)then
      if(.false.)then
        call writeinputdata(num_atoms,&
             lattice,xyzstruct,totalforce,&
             totalcharge,totalenergy,atomcharge,&
             atomenergy,elementsymbol,lperiodic)
      endif
!!
!!
!! CMH check hamiltonian information
!!

!!
      return
!!
 90   continue
      write(ounit,*)'ERROR in reading lattice vectors of structure ',inumber
      stop
      return
!!
 91   continue
      write(ounit,'(a,i10,a,i10)')'ERROR in reading atom ',num_atoms,' in structure ',inumber
      stop
      return
!!
 92   continue
      write(ounit,*)'ERROR in reading total charge of structure ',inumber
      stop
      return
!!
 93   continue
      write(ounit,*)'ERROR in reading total energy of structure ',inumber
      stop
      return
!!
 94   continue
      write(ounit,*)'ERROR in reading hamiltonian of structure ',inumber
      stop
      return
!!

      end
