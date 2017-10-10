!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - calcfunctions.f90
!! - calcpairfunctions.f90
!!
      subroutine splitcoulombforces(&
      num_atoms,atomcharge,xyzstruct,elecforce)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2,i3,i4
      integer num_atoms
!!      
      real*8 atomcharge(max_num_atoms)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 elecforce(3,max_num_atoms) ! out
!!
      real*8 distance
      real*8 drijdxyz(max_num_atoms,3)   ! internal
      real*8 invrij2 
      real*8 fscreen                          ! internal
      real*8 fscreenderiv                     ! internal
!!
!!
      elecforce(:,:)=0.0d0
!!
      if(num_atoms.gt.1) then
      do i3=1,num_atoms
      do i4=1,3
!!
      do i1=1,num_atoms-1
        do i2=i1+1,num_atoms
          distance=(xyzstruct(1,i1)-xyzstruct(1,i2))**2 + &
                   (xyzstruct(2,i1)-xyzstruct(2,i2))**2 + &
                   (xyzstruct(3,i1)-xyzstruct(3,i2))**2
          invrij2 =1.d0/distance
          distance=dsqrt(distance)
!! calculation of \frac{\partial r_{ij}}{\partial \alpha}
          drijdxyz(i3,i4)=0.0d0 ! initialization
          if(i1.eq.i3)then
            drijdxyz(i3,i4)=(xyzstruct(i4,i1)-xyzstruct(i4,i2))/distance
          elseif(i2.eq.i3)then
            drijdxyz(i3,i4)=-1.d0*(xyzstruct(i4,i1)-xyzstruct(i4,i2))/distance
          endif
!!
!! final force calculation
!!------------------------
          if(lscreen) then
            call getscreenfunctionforelectrostatics(&
                   distance,fscreen,&
                   fscreenderiv,drijdxyz(i3,i4))
!!
            elecforce(i4,i3)=elecforce(i4,i3) &
             - atomcharge(i1)*atomcharge(i2)* &
                (-invrij2*drijdxyz(i3,i4)*fscreen + & 
                  fscreenderiv/distance)
          else
            elecforce(i4,i3)=elecforce(i4,i3) &
               +invrij2*atomcharge(i1)*atomcharge(i2)*drijdxyz(i3,i4)
          endif
!!
        enddo ! i2
      enddo ! i1
!!
      enddo ! i4
      enddo ! i3
      endif 
!!
      return
      end
