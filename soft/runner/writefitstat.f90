!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fittingpair.f90

      subroutine writefitstat(ntrain,fitstat,fitstatf,fitstatq)
!!
      use nnflags
      use globaloptions
      use fileunits
!!
      implicit none
!!
      integer i1,i2,i3          ! internal
      integer ntrain            ! in
      integer num_atoms         ! internal
      integer num_pairs         ! internal
      integer idummy
      integer fitstat(ntrain)                  ! in
      integer fitstatf(3,max_num_atoms,ntrain) ! in
      integer fitstatq(max_num_atoms,ntrain)   ! in
!!
      real*8 edummy
!!
      write(ounit,*)'============================================================='
      write(ounit,*)'Fitting statistics:'
      if(lshort)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Short range energies used:'
        write(ounit,*)'           Point   Usage'
        do i1=1,ntrain
          write(ounit,'(a,2i8)')' NNstatE ',i1,fitstat(i1)
        enddo
        if(luseforces)then
!! we use a quick and dirty way to get the number of atoms of each structure from function.data here
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit)
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,*)'Short range forces used:' !'
          write(ounit,*)'           Point    Atom      fx      fy      fz'
          do i1=1,ntrain
            if(nn_type_short.eq.1)then
              read(symunit,*)num_atoms
              do i2=1,num_atoms
                read(symunit,*)idummy
              enddo
            elseif(nn_type_short.eq.2)then
              read(symunit,*)num_atoms,num_pairs
              do i2=1,num_pairs
                read(symunit,*)idummy
              enddo ! i2
            endif
            read(symunit,*)edummy
            do i2=1,num_atoms
              write(ounit,'(a,5i8)')' NNstatF ',i1,i2,(fitstatf(i3,i2,i1),i3=1,3)
            enddo
          enddo ! i1
          close(symunit)
        endif ! luseforces
      endif ! lshort
      if(lelec.and.(nn_type_elec.eq.1))then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Atomic charges used:'
        write(ounit,*)'           Point    Atom   Usage'
!! we use a quick and dirty way to get the number of atoms of each structure from function.data here
        open(symeunit,file='functione.data',form='formatted',status='old')
        rewind(symeunit) !'
        do i1=1,ntrain
          read(symeunit,*)num_atoms
          do i2=1,num_atoms
            read(symeunit,*)idummy
            write(ounit,'(a,5i8)')' NNstatQ ',i1,i2,fitstatq(i2,i1)
          enddo
          read(symeunit,*)edummy
        enddo ! i1
        close(symeunit)
      endif ! lelec
!!
      return
      end

