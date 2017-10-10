!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: - determine max_num_pairs
!!          - check if pair symmetry functions are given in input.nn in case of nn_type_short 2
!!
!! called by:
!! - initialization.f90
!! 
      subroutine paircount()
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer count_struct                     ! internal
      integer function_type_temp               ! internal
      integer i1                               ! internal
      integer nlattice                         ! internal
      integer num_atoms                        ! internal
      integer nn_type_short_local                    ! internal
      integer num_pairs                        ! in
      integer num_triplets                     ! in
      integer zelem(max_num_atoms)             ! internal 

      real*8  funccutoff_local                 ! internal 
      real*8  lattice(3,3)                     ! internal
      real*8  maxcutoff_local                  ! internal 
      real*8  xyzstruct(3,max_num_atoms)       ! internal
!!
      character*2  elementsymbol(max_num_atoms)! internal 
      character*2  elementtemp1                ! internal
      character*2  elementtemp2                ! internal
      character*40 keyword                     ! internal
      character*40 dummy
      character*7  dummy1

      logical lperiodic                        ! internal
!!
!! initializations:

      max_num_pairs    = 0
      max_num_triplets = 0
      maxcutoff_local  = 0.0d0
      funccutoff_local = 0.0d0
      nlattice         = 0
      num_atoms        = 0
      nn_type_short_local    = 1  ! use the same default as in readinput.f90 here

!! read input.nn file  to get the maxcutoff_local
      open(nnunit,file='input.nn')
      rewind(nnunit)
90    read(nnunit,*,END=80) keyword

      if(keyword.eq.'nn_type_short')then
        backspace(nnunit)
        read(nnunit,*)dummy1,nn_type_short_local
      endif

      goto 90
80    continue
      close(nnunit)

      if((nn_type_short_local.eq.2).or.(lnntb))then
        open(nnunit,file='input.nn')
        rewind(nnunit)

70      read(nnunit,*,END=60) keyword

          if(keyword.eq.'symfunction_hextoff')then
            backspace(nnunit)
            read(nnunit,*)dummy,dummy,dummy,dummy,function_type_temp
            if(function_type_temp.eq.1)then
              backspace(nnunit)
              read(nnunit,*)dummy,dummy,dummy,dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.2)then
              backspace(nnunit)
              read(nnunit,*)dummy,dummy,dummy,dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.3)then
              backspace(nnunit)
              read(nnunit,*)dummy,dummy,dummy,dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.4)then
              backspace(nnunit)
              read(nnunit,*)dummy,dummy,dummy,function_type_temp,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.5)then
              backspace(nnunit)
              read(nnunit,*)dummy,dummy,dummy,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.6)then
              backspace(nnunit)
              read(nnunit,*)dummy,dummy,dummy,function_type_temp,dummy,dummy,funccutoff_local
            else
              write(ounit,*)'Error: unknown symfunction in paircount'
              stop
            endif

          elseif(keyword.eq.'global_symfunction_hexton')then
            backspace(nnunit)
            read(nnunit,*)dummy,function_type_temp
            if(function_type_temp.eq.1)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.2)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.3)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.4)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.5)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.6)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,funccutoff_local
            else
              write(ounit,*)'Error: unknown symfunction in paircount'
              stop
            endif
          elseif(keyword.eq.'global_symfunction_s')then
            backspace(nnunit)
            read(nnunit,*)dummy,function_type_temp
            if(function_type_temp.eq.1)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.2)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.3)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.4)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.5)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.6)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,funccutoff_local
            else
              write(ounit,*)'Error: unknown symfunction in paircount'
              stop
            endif
          elseif(keyword.eq.'global_symfunction_dens')then
            backspace(nnunit)
            read(nnunit,*)dummy,function_type_temp
            if(function_type_temp.eq.1)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.2)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.3)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.4)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.5)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.6)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,funccutoff_local
            else
              write(ounit,*)'Error: unknown symfunction in paircount'
              stop
            endif







          elseif(keyword.eq.'global_pairsymfunction_short')then
            backspace(nnunit)
            read(nnunit,*)dummy,function_type_temp
            if(function_type_temp.eq.1)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.2)then
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.3)then   
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.4)then   
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.5)then    
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.6)then   
              backspace(nnunit)
              read(nnunit,*)dummy,function_type_temp,dummy,dummy,funccutoff_local
            else
              write(ounit,*)'Error: unknown symfunction in paircount'
              stop
            endif
!!
          elseif(keyword.eq.'element_pairsymfunction_short')then
            backspace(nnunit)
            read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp
            if(function_type_temp.eq.1)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.2)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.3)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.4)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.5)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.6)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,dummy,funccutoff_local
            else
              write(ounit,*)'Error: unknown symfunction in paircount'
              stop
           endif
!!
          elseif(keyword.eq.'pairsymfunction_short')then       
            backspace(nnunit)
            read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp
            if(function_type_temp.eq.1)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.2)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,funccutoff_local
            elseif(function_type_temp.eq.3)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.4)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,dummy,dummy,dummy,funccutoff_local
            elseif(function_type_temp.eq.5)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,funccutoff_local
            elseif(function_type_temp.eq.6)then
              backspace(nnunit)
              read(nnunit,*)dummy,elementtemp1,elementtemp2,function_type_temp,dummy,dummy,funccutoff_local
            else
              write(ounit,*)'Error: unknown symfunction in paircount'
              stop
           endif
          endif
!!
          maxcutoff_local=max(maxcutoff_local,funccutoff_local)
          goto 70
60        continue
        close(nnunit)

!! check if maxcutoff_local has been found properly
        if(maxcutoff_local.eq.0.0d0)then
          write(ounit,*)'Error: maxcutoff_local is zero in paircount'
          write(ounit,*)'Did you forget to specify symmetry functions???'
          stop !'
        endif
      endif ! --> nn_type_short_local = 2


!! Determine max_num_pairs:
      if((nn_type_short_local.eq.2).or.(lnntb))then
        count_struct=0
        open(dataunit,file='input.data',form='formatted')
          backspace(dataunit)
 10       continue
          read(dataunit,*,END=30) keyword
          if(keyword.eq.'begin') then
            count_struct= count_struct+1
            nlattice  = 0
            num_atoms = 0
            lperiodic =.false.
          endif
          if(keyword.eq.'lattice') then
            nlattice=nlattice+1
            backspace(dataunit)
            read(dataunit,*)keyword,(lattice(nlattice,i1),i1=1,3)
          endif
          if(keyword.eq.'atom') then
            backspace(dataunit)
            num_atoms=num_atoms+1
            read(dataunit,*)keyword,(xyzstruct(i1,num_atoms),i1=1,3),elementsymbol(num_atoms)
            call nuccharge(elementsymbol(num_atoms),zelem(num_atoms))
          endif
          if(keyword.eq.'end') then
            if(nlattice.eq.3)then
              lperiodic=.true.
            endif 
            if(lperiodic)then
               call translate(num_atoms,lattice,xyzstruct)
            endif
!! determine num_pairs
            call getnumpairs(num_atoms,num_pairs,zelem,&
              maxcutoff_local,lattice,xyzstruct,lperiodic)
!!          
            max_num_pairs=max(max_num_pairs,num_pairs)
            call getnumtriplets(num_atoms,num_triplets,zelem,&
              maxcutoff_local,lattice,xyzstruct,lperiodic)
            max_num_triplets=max(max_num_triplets,num_triplets)
          endif
          goto 10
 30       continue
        close(dataunit)
      endif ! nn_type_short_local.eq.2     



      if((max_num_pairs.eq.0).and.((nn_type_short_local.eq.2).or.(lnntb)))then
        write(ounit,*)'Error: max_num_pairs is 0 in paircount.f90 ',max_num_pairs
        stop
      endif
!!
      return
      end
