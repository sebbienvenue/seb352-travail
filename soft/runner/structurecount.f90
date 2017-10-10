!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - initialization.f90
!!
      subroutine structurecount()
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer n_begin
      integer n_end
      integer n_atom
      integer n_pairs                 ! internal
      integer i1
      integer idummy                  ! internal
!!
      real*8 edummy                   ! internal
!!
      character*20 keyword
!!
!! initialization
      n_begin       = 0
      n_end         = 0
      n_atom        = 0
      max_num_atoms = 0
!!
      if((mode.eq.1).or.(mode.eq.3))then
        open(dataunit,file='input.data',form='formatted')
 10       continue
          read(dataunit,*,END=30) keyword
          if(keyword.eq.'begin') n_begin=n_begin+1
          if(keyword.eq.'atom') then
            n_atom=n_atom+1
            max_num_atoms=max(max_num_atoms,n_atom)
          endif
          if(keyword.eq.'end') then
            n_end  = n_end+1
            n_atom = 0
          endif
          goto 10
 30     continue
        close(dataunit)
!!
        if(n_begin.ne.n_end) then
          write(ounit,*)'ERROR: in input.data'
          write(ounit,*)'n_begin,n_end = ',n_begin,n_end
          stop
        endif
        totnum_structures = n_begin
!!
      elseif(mode.eq.2)then

        if(lshort)then
!! check function.data file
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit) !'
 40       continue
          if(nn_type_short.eq.1)then
            read(symunit,*,END=41) n_atom 
            max_num_atoms=max(max_num_atoms,n_atom)
            do i1=1,n_atom
              read(symunit,*)idummy
            enddo
          elseif(nn_type_short.eq.2)then
            read(symunit,*,END=41) n_atom,n_pairs
            max_num_atoms=max(max_num_atoms,n_atom)
            do i1=1,n_pairs 
              read(symunit,*)idummy
            enddo
          else
            write(ounit,*)'ERROR: unknown nn_type_short in structurecount'
            stop
          endif
          read(symunit,*)edummy 
          goto 40
 41       continue
          close(symunit)
!! check testing.data file
          open(tymunit,file='testing.data',form='formatted',status='old')
          rewind(tymunit)
 50       continue
          if(nn_type_short.eq.1)then
            read(tymunit,*,END=51) n_atom 
            max_num_atoms=max(max_num_atoms,n_atom)
            do i1=1,n_atom
              read(tymunit,*)idummy
            enddo
          else
            read(tymunit,*,END=51) n_atom,n_pairs
            max_num_atoms=max(max_num_atoms,n_atom)
            do i1=1,n_pairs 
              read(tymunit,*)idummy
            enddo
          endif
          read(tymunit,*)edummy 
          goto 50
 51       continue
          close(tymunit)
!!
        elseif(lelec.and.(nn_type_elec.eq.1))then ! get max_num_atoms from here only if no short range NN is used
          open(symeunit,file='functione.data',form='formatted',status='old')
          rewind(symeunit) !'
 60       continue
          read(symeunit,*,END=61) n_atom 
          max_num_atoms=max(max_num_atoms,n_atom)
          do i1=1,n_atom
            read(symeunit,*)idummy
          enddo
          read(symeunit,*)edummy 
          goto 60
 61       continue
          close(symeunit)
!! check testinge.data file
          open(tymeunit,file='testinge.data',form='formatted',status='old')
          rewind(tymeunit) !'
 70       continue
          read(tymeunit,*,END=71) n_atom 
          max_num_atoms=max(max_num_atoms,n_atom)
          do i1=1,n_atom
            read(tymeunit,*)idummy
          enddo
          read(tymeunit,*)edummy 
          goto 70
 71       continue
          close(tymeunit)
        endif ! lshort
        if(lnntb) then
          if(nntb_flag(3))then
            max_num_atoms = 3
          endif
        endif
!!
      else
        write(ounit,*)'Error: illegal runner_mode detected in structurecount ',mode
        stop
      endif
!!
!!
      return
      end
