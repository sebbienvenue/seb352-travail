!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - calcfunctions.f90
!!            - calcpairfunctions.f90
!!
!! this subroutine is for the non-periodic case
!!
      subroutine getcoulombforcesone(max_num_neighbors_elec,&
         num_atoms,dchargedxyz,&
         nnatomcharge,xyzstruct,nnelecforce)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer num_atoms                                                 ! in
      integer i1,i2,i3,i4                                               ! internal
      integer max_num_neighbors_elec                                    ! in
!!
      real*8 nnatomcharge(max_num_atoms)                                ! in
      real*8 xyzstruct(3,max_num_atoms)                                 ! in
      real*8 nnelecforce(3,max_num_atoms)                               ! out
      real*8 dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)      ! in
      real*8 distance                                                   ! internal
      real*8 invrij2                                                    ! internal
      real*8 drijdxyz(max_num_atoms,3)                                  ! internal
      real*8 fscreen                          ! internal
      real*8 fscreenderiv                     ! internal
!!
!!
      nnelecforce(:,:)=0.0d0 
!!
      if(num_atoms.gt.1) then
!!      do i3=1,num_atoms
      do i3=1,num_atoms ! FIXME
      do i4=1,3
!!
!! caution: do we need double counting for the forces???
      do i1=1,num_atoms
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
            nnelecforce(i4,i3)=nnelecforce(i4,i3) &
             - invrij2*((dchargedxyz(i1,i3,i4)*nnatomcharge(i2) &
                        + nnatomcharge(i1)*dchargedxyz(i2,i3,i4))*&
                        distance*fscreen - nnatomcharge(i1)* &
                        nnatomcharge(i2)*drijdxyz(i3,i4)* &
                        fscreen+fscreenderiv*nnatomcharge(i1)* &
                        nnatomcharge(i2)*distance)
          else
            nnelecforce(i4,i3)=nnelecforce(i4,i3) &
               - invrij2*((dchargedxyz(i1,i3,i4)*nnatomcharge(i2) &
                         + nnatomcharge(i1)*dchargedxyz(i2,i3,i4))*distance&
               -nnatomcharge(i1)*nnatomcharge(i2)*drijdxyz(i3,i4))
          endif
!!
          enddo ! i2
        enddo ! i1
      enddo ! i4
      enddo ! i3
!! 
      endif ! num_atoms.gt.1
!!
      return
      end