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
      subroutine readonestructure_ham(num_atoms,&
        zelem,num_atoms_element,lattice,&
        totalcharge,totalenergy,&
        atomcharge,atomenergy,xyzstruct,&
        totalforce,hextoff_value,elementsymbol,lperiodic,structnum)
!!
      use fileunits
      use globaloptions
      use nnflags
      use basismod
!!
      implicit none
!!
      integer nlattice                         ! internal 
      integer num_atoms                        ! out 
      integer zelem(max_num_atoms)             ! out 
      integer num_atoms_element(nelem)         ! out
      integer structnum                        ! out
      integer i,i1,i2                          ! internal 
      integer zid                              ! internal
      integer zatom_id_temp                    ! internal
      integer matrixelem_count                 ! internal
      integer num_subblock                     ! internal
      integer at1                              ! internal
      integer at2                              ! internal
      integer n1                               ! internal
      integer n2                               ! internal
      integer l1                               ! internal
      integer l2                               ! internal
      integer m1                               ! internal
      integer m2                               ! internal
      integer hnum1                            ! internal
      integer hnum2                            ! internal
      integer hcount                           ! internal
      integer count_ham                        ! internal
      integer count_hoverlap                   ! internal
      integer count_hexton                     ! internal
      integer count_hextoff                    ! internal
      integer count_hdens                      ! internal
      integer hextoff_line                     ! internal
      integer hextoff_column                   ! internal
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
      real*8 zcoord                            ! internal
      real*8 hextoff_value(maxnum_basis,maxnum_basis) ! out 
!!
      character*20 keyword                     ! internal
      character*2 elementsymbol(max_num_atoms) ! out
      character*10 matrixelem                  ! internal
      character*6 dummy1, dummy2               ! internal
!!
      logical lbegin                           ! internal, true if 'begin' flag of structure found
      logical lend                             ! internal, true if 'end' flag of structure found
      logical lcharge                          ! internal, true if 'charge' keyword is found
      logical lenergy                          ! internal, true if 'energy' keyword is found
      logical lperiodic                        ! out, true if 3 lattice vectors are found 
      logical originflag                       ! internal
      logical zflag                            ! internal
      logical found_hextoff(maxnum_basis,maxnum_basis) ! internal
!!
!! initializations
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
      originflag           = .false.
      zflag                = .false.
      zid                  = 0
      matrixelem_count     = 0
      count_ham            = 0
      count_hoverlap       = 0
      count_hexton         = 0
      count_hextoff        = 0
      count_hdens          = 0
      hcount               = 0
      if(nntb_flag(3)) then ! Hextoff case
       found_hextoff(:,:)=.false.
      endif
!!
 10   continue
      read(dataunit,*)keyword
!!
!! read in all information from input.data for this structure 
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
      if(keyword.eq.'Ham')then     
        count_ham=count_ham+1
        write(ounit,*)'ERROR: reading of ham not yet implemented in readonestructure'
        stop
      endif
      if(keyword.eq.'Hoverlap')then
        count_hoverlap=count_hoverlap+1
        write(ounit,*)'ERROR: reading of overlap not yet implemented in readonestructure'
        stop
      endif
      if(keyword.eq.'Hexton')then  
        count_hexton=count_hexton+1
        write(ounit,*)'ERROR: reading of hexton not yet implemented in readonestructure'
        stop
      endif
      if(keyword.eq.'Hextoff')then 
        count_hextoff=count_hextoff+1
!! check if the following 2 lines are really needed
        backspace(dataunit)
        read(dataunit,*)keyword,hextoff_line,hextoff_column
        backspace(dataunit)
        read(dataunit,*)keyword,hextoff_line,hextoff_column,hextoff_value(hextoff_line,hextoff_column)
        found_hextoff(hextoff_line,hextoff_column)=.true.
      endif
      if(keyword.eq.'Hdens')then   
        count_hdens=count_hdens+1
        !!write(ounit,*)'ERROR: reading of hdens not yet implemented in readonestructure'
        !!stop
      endif
      if((keyword.eq.'end').and.(lbegin)) lend=.true. 
      if(keyword.ne.'end') goto 10

      if(nlattice.eq.3) then
        lperiodic=.true.
      endif
!!
!! translate atoms back into cell if necessary
      if(lperiodic)then
        call translate(num_atoms,lattice,xyzstruct)
      endif
!!
!! count the number of atoms of each element
      do i=1,num_atoms
        num_atoms_element(elementindex(zelem(i)))=&
           num_atoms_element(elementindex(zelem(i)))+1
      enddo
!!
!!    Lets look at the structure information!!! CMH
!!      write(*,*) 'num atoms', num_atoms
!!      write(*,*) 'coords'
!!      write(*,*) xyzstruct(:,1),elementsymbol(1)
!!      write(*,*) xyzstruct(:,2),elementsymbol(2)
!!      write(*,*) xyzstruct(:,3),elementsymbol(3)
      if(lnntb.and.(mode.eq.1))then
        if((num_atoms.ne.2).and.(((nntb_flag(0)).or.(nntb_flag(1)).or.(nntb_flag(2)))))then
          write(ounit,*) 'Error: For Overlap and Hexton in mode 1 there must be 2 atoms in the system',structnum
          stop
        elseif((num_atoms.ne.3).and.(nntb_flag(3)))then ! Hextoff
          write(ounit,*) 'Error: For Hextoff in mode 1 there must be only 3 atoms in the system',structnum
          stop
        elseif((num_atoms.le.2).and.(nntb_flag(4)))then
          write(ounit,*) 'Error: For NNTB density there must be 2 atoms or more in the system',structnum
          stop
        endif
      endif
!!
!! check standard orientation for hext
!! Now lets check the structure for origin and z axis alignment - Not sure if we will need this for the Density part
      if(lnntb.and.(nntb_flag(3)).and.(mode.lt.3))then
        do i=1,num_atoms
          if(((abs(xyzstruct(1,i))).le.(0.0E-10))&
            .and.(abs((xyzstruct(2,i))).le.(0.0E-10))&
            .and.(abs((xyzstruct(3,i))).le.(0.0E-10))) then
            originflag = .true.
!! define what element is the atom at the origin
            if(structnum.eq.1)then
              originatom_id = zelem(i)
              if((zelem(1) /= hextoff_training_triplet(1)).or.&
                 (zelem(1).lt.zelem(2))) then
                write(ounit,*) 'Error: atom type 1 is lower than atom type 2, and error in structure', structnum
                stop
              endif
!! check if element of the atom at the origin is always the same
            else
              if(zelem(i).ne.originatom_id)then
                write(ounit,*) zelem(i), originatom_id
                write(ounit,*) 'Error: atom type at origin inconsistent in stucture ', structnum
                stop
              endif
            endif
!!           write(*,*) originatom_id
!!           write(*,*) 'atom at origin is ', i
!! check if atom is at the positive z-axis
           elseif((abs(xyzstruct(1,i)).le.0.0E-14)&
             .and.(abs(xyzstruct(2,i)).le.0.0E-14)&
             .and.(xyzstruct(3,i).gt.0.5)) then
             if(zflag) then
!! check in case of 2 atoms on positive z axis which is closer to the origin
               if (xyzstruct(3,i).lt.zcoord) then
                 zcoord = xyzstruct(3,i)
                 zid = i
                 zatom_id_temp = zelem(i)
                 !write(*,*) zatom_id_temp
               endif
             else
               zflag = .true.
               zcoord = xyzstruct(3,i)
               zid = i
               zatom_id_temp = zelem(i)
               !!write(*,*) zatom_id_temp
             endif
           endif
        enddo
!! check for all structures if the element of the second atom in pair is the same like in first structure
        if (structnum.gt.1)then
           if(zatom_id_temp.ne.zatom_id)then
             write(*,*) 'Error: Z axis id not consistent in structure ', structnum
             stop
           endif
        else
          zatom_id = zatom_id_temp
        endif
!!      write(*,*) 'z axis atom ', zid
!!
!! calculate the number of subblock matrix elements to be expected
!! and check if we have the correct number of lines in input.data
!! check if all required matrix elements have been found exactly once
        if(nntb_flag(3)) then ! Hextoff case
          num_subblock = num_basis(elementindex(zelem(1)))*num_basis(elementindex(zelem(2)))
          if(num_subblock.ne.count_hextoff)then
            write(ounit,*)'ERROR: number of hextoff matrix elements incorrect in structure ',structnum
            stop
          endif
!! check if all matrix elements have been found
          do i1=1,num_basis(elementindex(zelem(1)))
            do i2=1,num_basis(elementindex(zelem(2)))
              if(found_hextoff(i1,i2).eqv..false.)then
                write(ounit,*)'ERROR: Hextoff has not been found in input.data ',i1,i2
                write(ounit,*)'structure number is ',structnum
                stop
              endif
            enddo
          enddo
        endif






      endif ! lnntb

      return
      end
