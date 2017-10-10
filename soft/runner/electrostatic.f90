!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - calcfunctions.f90
!!            - calcpairfunctions.f90
!!            - getallelectrostatic.f90
!!            - prediction.f90
!!            - predictionpair.f90
!!
      subroutine electrostatic(num_atoms,&
                    atomcharge,xyzstruct,elecenergy)
!!
      use globaloptions
      use fileunits
!!
      implicit none
!!
      integer num_atoms
      integer i1,i2
!!
      real*8 xyzstruct(3,max_num_atoms)               ! in
      real*8 atomcharge(max_num_atoms)                ! in
      real*8 elecenergy                               ! out
      real*8 distance                                 ! internal
      real*8 fscreen                                  ! internal
      real*8 fscreenderiv                             ! internal
!!
!!
!! For this subroutine the use of Hartree energy and Bohr length units is mandatory!
!! 1 Hartree = \frac{\hbar^2}{m_e a_0^2} = \frac{e^2}{4\pi \epsilon_0 a_0}
!!
      elecenergy     =0.0d0
!!
      if(num_atoms.gt.1) then
      do i1=1,num_atoms
        do i2=i1+1,num_atoms
          distance=(xyzstruct(1,i1)-xyzstruct(1,i2))**2 + &
                   (xyzstruct(2,i1)-xyzstruct(2,i2))**2 + &
                   (xyzstruct(3,i1)-xyzstruct(3,i2))**2
          distance=dsqrt(distance)
!!
          if(lscreen) then
!! Question JB: Why is drijdxyz put in as 0.0d0?
            call getscreenfunctionforelectrostatics(&
                   distance,fscreen,&
                   fscreenderiv,0.0d0)
!!
             elecenergy=elecenergy + atomcharge(i1)*atomcharge(i2)/ &
                           distance*fscreen
          else
            elecenergy=elecenergy + atomcharge(i1)*atomcharge(i2)/distance
          endif
        enddo ! i2
      enddo ! i1 
      endif
!!
      return
!!
      end
