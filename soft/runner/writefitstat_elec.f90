!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:

      subroutine writefitstat_elec(ntrain,fitstatq)
!!
      use globaloptions
      use fileunits
!!
      implicit none
!!
      integer i1,i2             ! internal
      integer ntrain            ! in
      integer num_atoms         ! internal
      integer idummy
      integer fitstatq(max_num_atoms,ntrain)   ! in
!!
      real*8 edummy
!!
      write(ounit,*)'============================================================='
      write(ounit,*)'Fitting statistics:'
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
!!
      return
      end

