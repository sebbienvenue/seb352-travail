!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - checkonestructure.f90
!!
      subroutine writeinputdata(num_atoms,&
             lattice,xyzstruct,totalforce,&
             totalcharge,totalenergy,atomcharge,&
             atomenergy,elementsymbol,lperiodic)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer num_atoms
      integer i,j
!!
      real*8 lattice(3,3)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 totalforce(3,max_num_atoms)
      real*8 totalenergy
      real*8 totalcharge
      real*8 atomenergy(max_num_atoms)
      real*8 atomcharge(max_num_atoms)
!!
      character*2 elementsymbol(max_num_atoms)
!!
      logical lperiodic
!!
      if(ldebug)then
!!
      write(ounit,*)'-------------------------------------------------------------'
      if(lperiodic) then
        write(ounit,*)'Found periodic structure with atoms ',num_atoms
      else
        write(ounit,*)'Found nonperiodic structure with atoms ',num_atoms
      endif
      write(ounit,*)'begin '
      if(lperiodic)then
        do i=1,3
          write(ounit,'(a8,x,3f14.8)')' lattice',(lattice(i,j),j=1,3)
        enddo
      endif
      do i=1,num_atoms
        write(ounit,'(a5,x,3f14.8,x,a2,5f14.8)')' atom',&
              (xyzstruct(j,i),j=1,3),&
              elementsymbol(i),atomcharge(i),atomenergy(i),(totalforce(j,i),j=1,3)
      enddo
      write(ounit,'(a,f20.8)')' energy ',totalenergy
      write(ounit,'(a,f20.8)')' charge ',totalcharge
      write(ounit,*)'end '
!!
      endif ! ldebug'
!!
      return
      end
