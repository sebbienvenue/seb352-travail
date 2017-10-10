!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getewaldenergy_para.f90
!!
      subroutine ewaldself_para(max_num_neighbors_elec,neighboridx_elec,&
                 num_neighbors_elec,natoms,atomindex,&
                 atomcharge,&
                 dchargedxyz,eselfforce,eself,ldoforces)
!!
      use globaloptions
!!
      implicit none
!!
      integer max_num_neighbors_elec                      ! in
      integer num_neighbors_elec(natoms)                  ! in 
      integer neighboridx_elec(natoms,0:max_num_neighbors_elec) ! in
      integer natoms                                      ! in
      integer atomindex(natoms)                           ! in
!!
      integer i1,i2,i3                                    ! internal
!!
      real*8 atomcharge(max_num_atoms)                    ! in
      real*8 eselfforce(3,max_num_atoms)                  ! out
      real*8 eself                                        ! out
      real*8 sqrtpiinv
      parameter(sqrtpiinv=0.564189583d0)
      real*8 dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)    ! in  
!!
      logical ldoforces                                   ! in
!!
      eself=0.0d0
!!
      do i1=1,natoms
        eself=eself+atomcharge(atomindex(i1))*atomcharge(atomindex(i1)) 
      enddo
!!
      if(ldoforces)then
        do i2=1,3
         do i3=1,natoms
          do i1=0,num_neighbors_elec(i3)
            eselfforce(i2,neighboridx_elec(i2,i1))=eselfforce(i2,neighboridx_elec(i2,i1)) &
              -2.d0*atomcharge(atomindex(i3))*dchargedxyz(atomindex(i3),i1,i2)
          enddo ! i1
         enddo ! i2
        enddo ! i3
      endif
!!
      eself=-1.d0*eself*ewaldalpha*sqrtpiinv ! note sign!
      eselfforce(:,:)=-1.d0*eselfforce(:,:)*ewaldalpha*sqrtpiinv ! note sign
!!
      return
      end
