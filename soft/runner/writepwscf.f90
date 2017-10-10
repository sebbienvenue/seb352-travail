!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - prediction.f90
!!
      subroutine writepwscf(num_atoms,&
           lattice,xyzstruct,elementsymbol,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer num_atoms                ! in
      integer i1,i2                    ! internal
!! 
      real*8 xyzstruct(3,max_num_atoms) ! in
      real*8 lattice(3,3)               ! in
!!
      character*2 elementsymbol(max_num_atoms) ! in
!!
      logical lperiodic                ! in
!!
!!
      open(pwunit,file='pwscf.in',form='formatted',status='replace')
!!
      write(pwunit,'(a8)')' &system'
      write(pwunit,'(a9)')' ibrav= 0'
      write(pwunit,'(a6,x,i5,x,a2)')' nat= ',num_atoms
      write(pwunit,'(a7,x,i5,x,a2)')' ntyp= ',nelem
!!      write(pwunit,'(a)')'    ecutwfc =30.0,nosym=.true.'
!!      write(pwunit,'(a)')'    occupations='smearing''
!!      write(pwunit,'(a)')'    degauss=0.0073501d0'
!!      write(pwunit,'(a)')'    smearing='fermi-dirac''
      write(pwunit,'(a2)')' /'
!!
      write(pwunit,'(a)')'ATOMIC_SPECIES'
      do i1=1,nelem
        write(pwunit,'(a2)')elementsymbol(i1)
      enddo
!!
      write(pwunit,*)'CELL_PARAMETERS'
      if(lperiodic)then
        do i1=1,3
          write(pwunit,'(3f14.8)')(lattice(i1,i2),i2=1,3) 
        enddo
      else
        write(pwunit,'(3f14.8)') 30.0d0,0.0d0,0.0d0
        write(pwunit,'(3f14.8)') 0.0d0,30.0d0,0.0d0
        write(pwunit,'(3f14.8)') 0.0d0,0.0d0,30.0d0
      endif ! lperiodic
!!
      write(pwunit,*)'ATOMIC_POSITIONS (bohr)'
      do i1=1,num_atoms
        write(pwunit,'(a2,x,3f14.8)')elementsymbol(i1),(xyzstruct(i2,i1),i2=1,3)
      enddo
      close(pwunit)
!!
      return
      end
