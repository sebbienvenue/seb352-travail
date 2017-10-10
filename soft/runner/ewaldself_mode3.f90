!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: Calculate self-energy term of standard Ewald summation, parallel case
!!
!! called by:
!! - getewald.f90
!!
      subroutine ewaldself_mode3(max_num_neighbors_elec,&
        num_atoms,zelem,invneighboridx_elec,&
        natoms,atomindex,num_funcvalues_elec,&
        nnatomcharge,dchargedsfunc,dsfuncdxyze,&
        eselfforce,eself,ldoforces_local)
!!
      use globaloptions
      use fileunits
      use nnflags
!!
      implicit none
!!
      integer max_num_neighbors_elec                      ! in
      integer invneighboridx_elec(natoms,max_num_atoms)   ! in
      integer natoms                                      ! in
      integer num_atoms                                   ! in
      integer atomindex(natoms)                           ! in
      integer num_funcvalues_elec(nelem)                  ! in
      integer zelem(max_num_atoms)                        ! in
!!
      integer i1,i2,i3,i4                                 ! internal
!!
      real*8 nnatomcharge(max_num_atoms)                  ! in
      real*8 eselfforce(3,max_num_atoms)                  ! out
      real*8 eself                                        ! out
      real*8 sqrtpiinv
      parameter(sqrtpiinv=0.564189583d0)
      real*8 dsfuncdxyze(maxnum_funcvalues_elec,natoms,0:max_num_neighbors_elec,3) ! in
      real*8 dchargedsfunc(natoms,maxnum_funcvalues_elec) ! in
      real*8 tempsum                                      ! internal
!!
      logical ldoforces_local                             ! in
!!
      eself=0.0d0
      eselfforce(:,:)=0.0d0
!!
      do i1=1,natoms
        eself=eself+nnatomcharge(atomindex(i1))*nnatomcharge(atomindex(i1)) 
      enddo
!!
      if(ldoforces_local)then
        if((nn_type_elec.ne.3).and.(nn_type_elec.ne.4))then
          do i2=1,3
            do i3=1,num_atoms
              do i1=1,natoms
                tempsum=0.0d0
                do i4=1,num_funcvalues_elec(elementindex(zelem(atomindex(i1))))
                  tempsum=tempsum+dchargedsfunc(i1,i4)*dsfuncdxyze(i4,i1,invneighboridx_elec(i1,i3),i2)
                enddo ! i4
                eselfforce(i2,i3)=eselfforce(i2,i3)& 
                  -2.d0*nnatomcharge(atomindex(i1))*tempsum
              enddo ! i1
            enddo ! i2
          enddo ! i3
        endif
      endif
!!
      eself=-1.d0*eself*ewaldalpha*sqrtpiinv ! note sign!
      eselfforce(:,:)=-1.d0*eselfforce(:,:)*ewaldalpha*sqrtpiinv ! note sign
!!
      return
      end
