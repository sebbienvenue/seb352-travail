!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:

      subroutine writefitstat_short(ntrain,fitstat,fitstatf)
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
!!
      real*8 edummy
!!
      write(ounit,*)'============================================================='
      write(ounit,*)'Fitting statistics:'
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
!!
      return
      end

