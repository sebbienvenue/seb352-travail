!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate a part of the electrostatic energy and forces for non-periodic case for a block of atoms

!! called by: 
!!
      subroutine coulomb_para(natoms,num_atoms,atomindex,zelem,&
        num_funcvalues_elec,max_num_neighbors_elec,invneighboridx_elec,&
        nnatomcharge,xyzstruct,dsfuncdxyze,dchargedsfunc,&
        maxcutoff_local,nnelecenergy,nnelecforce,ldoforces)
!!
      use fileunits
      use globaloptions
      use nnflags
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5                                            ! internal
      integer natoms                                                    ! in
      integer atomindex(natoms)                                         ! in
      integer num_atoms                                                 ! in
      integer num_funcvalues_elec(nelem)                                ! in
      integer max_num_neighbors_elec                                    ! in
      integer invneighboridx_elec(natoms,max_num_atoms)                 ! in
      integer zelem(max_num_atoms)                                      ! in
!!
      real*8 nnatomcharge(max_num_atoms)                                ! in
      real*8 xyzstruct(3,max_num_atoms)                                 ! in
      real*8 nnelecenergy                                               ! in/out
      real*8 nnelecforce(3,max_num_atoms)                               ! in/out
      real*8 fscreen                                                    ! internal
      real*8 fscreenderiv                                               ! internal
      real*8 dsfuncdxyze(maxnum_funcvalues_elec,natoms,0:max_num_neighbors_elec,3) ! in
      real*8 dchargedsfunc(natoms,maxnum_funcvalues_elec)               ! in
      real*8 rij                                                        ! internal
      real*8 distance                                                   ! internal
      real*8 invrij                                                     ! internal
      real*8 drijdxyz(max_num_atoms,3)                                  ! internal
      real*8 tempsum                                                    ! internal
      real*8 maxcutoff_local                                            ! in
!!
      logical ldoforces                                                 ! in
!!
!!======================================================================
!! calculate Coulomb energy  
!!======================================================================
      do i3=1,natoms
        do i4=1,num_atoms
          if(i4.gt.atomindex(i3))then
            rij=(xyzstruct(1,i4)-xyzstruct(1,atomindex(i3)))**2 + &
                (xyzstruct(2,i4)-xyzstruct(2,atomindex(i3)))**2 + &
                (xyzstruct(3,i4)-xyzstruct(3,atomindex(i3)))**2
            rij=dsqrt(rij)
            invrij =1.d0/rij
            if(lscreen) then
!! Question JB: Why is drijdxyz put in as 0.0d0?
              call getscreenfunctionforelectrostatics(&
                rij,fscreen,fscreenderiv,0.0d0)
              nnelecenergy=nnelecenergy &
                + nnatomcharge(atomindex(i3))*nnatomcharge(i4)*invrij*fscreen
            else
              nnelecenergy=nnelecenergy &
                + nnatomcharge(atomindex(i3))*nnatomcharge(i4)*invrij
!              write(ounit,'(a,2i5,f14.8)')'adding ',i4,atomindex(i3),nnatomcharge(atomindex(i3))*nnatomcharge(i4)*invrij
            endif
          endif ! (i4.ne.atomindex(i3))then
        enddo ! i4
      enddo ! i3
!!
!!======================================================================
!! loop over force components i1 of all atoms i2 
!!======================================================================
      if(ldoforces)then
        do i1=1,3
          do i2=1,num_atoms
!!======================================================================
!! loop over all natoms atoms i3 of this block 
!!======================================================================
            do i3=1,natoms
!!======================================================================
!! loop over all neighboring atoms i4  
!!======================================================================
!! TODO: check: must this be over num_neighbors instead(?) Maybe not because of infinite range of Coulomb 
              do i4=1,num_atoms
                if(i4.ne.atomindex(i3))then
                  rij=(xyzstruct(1,i4)-xyzstruct(1,atomindex(i3)))**2 + &
                      (xyzstruct(2,i4)-xyzstruct(2,atomindex(i3)))**2 + &
                      (xyzstruct(3,i4)-xyzstruct(3,atomindex(i3)))**2
                  rij=dsqrt(rij)
                  invrij =1.d0/rij
!!======================================================================
!! calculation of \frac{\partial r_{ij}}{\partial \alpha}
!!======================================================================
!! TODO CHECK: Is the sign correct here?:
                  drijdxyz(i2,i1)=0.0d0 ! initialization
                  if(i2.eq.i4)then
                    drijdxyz(i2,i1)=(xyzstruct(i1,i4)-xyzstruct(i1,atomindex(i3)))*invrij
                  elseif(i2.eq.atomindex(i3))then
                    drijdxyz(i2,i1)=-1.d0*(xyzstruct(i1,i4)-xyzstruct(i1,atomindex(i3)))*invrij
                  endif
!!
!!======================================================================
!! calculate dchargedxyz implicitly as
!! \sum_i5 \frac{\partial Q(i3)}{\partial G(i3)_{i5}}*\frac{\partial G(i3)_{i5}}{\partial alpha(i2,i1)} 
!! Caution: dsfuncdxyze is only defined for neighbors, therefore check distance
!!======================================================================
                  tempsum=0.0d0
                  if((nn_type_elec.ne.3).and.(nn_type_elec.ne.4))then ! don't do this if we use fixed charges
                    distance=(xyzstruct(1,i2)-xyzstruct(1,atomindex(i3)))**2 + &
                             (xyzstruct(2,i2)-xyzstruct(2,atomindex(i3)))**2 + &
                             (xyzstruct(3,i2)-xyzstruct(3,atomindex(i3)))**2
                    if(distance.le.maxcutoff_local**2)then 
                      do i5=1,num_funcvalues_elec(elementindex(zelem(atomindex(i3))))
                       tempsum=tempsum+dchargedsfunc(i3,i5)*dsfuncdxyze(i5,i3,invneighboridx_elec(i3,i2),i1)
                      enddo
                    endif
                  endif
!!
                  if(lscreen) then
                    call getscreenfunctionforelectrostatics(rij,fscreen,fscreenderiv,drijdxyz(i2,i1))
!!
                    nnelecforce(i1,i2)=nnelecforce(i1,i2) &
                      +nnatomcharge(i4)*invrij*fscreen &
                      *(0.5d0*nnatomcharge(atomindex(i3))*drijdxyz(i2,i1)*invrij-tempsum) &
                      -nnatomcharge(atomindex(i3))*nnatomcharge(i4)*0.5d0*invrij*fscreenderiv
                  else
                    nnelecforce(i1,i2)=nnelecforce(i1,i2) &
                      +nnatomcharge(i4)*invrij &
                      *(0.5d0*nnatomcharge(atomindex(i3))*drijdxyz(i2,i1)*invrij-tempsum)
                  endif
!                  write(ounit,*)'nnelecforce ',i1,i2,i3,i4,nnelecforce(i1,i2)
                endif ! (i4.ne.atomindex(i3))
              enddo ! i4
            enddo ! i3
          enddo ! i2
        enddo ! i1
      endif ! ldoforces 
!!
      return
      end
