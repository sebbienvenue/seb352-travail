!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getewaldenergy_para.f90 
!!
      subroutine  ewaldrecip_para(max_num_neighbors_elec,&
            num_atoms,xyzstruct,&
            atomcharge,erecipforce,erecip,lattice,dchargedxyz,&
            ldoforces)
!!
      use mpi_mod
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i,j
      integer i0,numk,totnumk
      integer k_start,k_end
      integer i1,i2,i3
      integer i4,i5
      integer max_num_neighbors_elec
!!
      integer num_atoms 
!!
      real*8 atomcharge(max_num_atoms)
      real*8 erecip
      real*8 determinant                           ! internal
      real*8 lattice(3,3)
      real*8 reclattice(3,3)
      real*8 kvec(3)
      real*8 ksquared
      real*8 kri,krj
      real*8 xyzstruct(3,max_num_atoms)
      real*8 temp1
      real*8 twopi
      real*8 factor
      real*8 sumsin,sumcos
      real*8 sqsumsin,sqsumcos
      real*8 sumsinj,sumcosj
      real*8 erecipforce(3,max_num_atoms)         ! out
      real*8 dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)    ! in  
      real*8 sumcosdev
      real*8 sumsindev
!!
      logical ldoforces
!!
!! initializations
      erecip=0.0d0
      twopi=6.283185307d0
!!
!! calculate the reciprocal lattice vectors
      call getreclat(lattice,reclattice,determinant)
!!
!! generate the reciprocal lattice vectors
!! later we should implement here a more clever way that automatically adapts to the
!! lengths of the reciprocal lattice vectors to really have a sphere in reciprocal space
!!
!! determine the total number of loops over kmax
      totnumk=(2*ewaldkmax+1)*(2*ewaldkmax+1)*(2*ewaldkmax+1)
!!      write(ounit,*)'number of reciprocal points ',totnumk
!! determine the number of loops for this process
      if((mpirank+1).le.(mod(totnumk,mpisize)))then
        numk=int(totnumk/mpisize)+1
      else
        numk=int(totnumk/mpisize)
      endif
!!      write(ounit,*)'process ',mpirank,' does ',numk,' reciprocal points'
!!' 
!! calculate the first and last point to be calculated by this process
      if((mpirank).le.(mod(totnumk,mpisize)))then
        k_start=int(totnumk/mpisize)*mpirank+1
      else
        k_start=int(totnumk/mpisize)*mpirank+1
      endif
      if(mpirank.gt.0)k_start=k_start+1
      k_end=k_start+numk-1
!!      write(ounit,'(a8,i4,a18,i6,a4,i6,a2,i6,a8)')&
!!        'process ',mpirank,' calculates point ',k_start,' to ',k_end,&
!!        '(',numk,' points)'
!!
      do i0=k_start,k_end
!!
        i3 = mod(i0,(2*ewaldkmax+1))
        if(i3.eq.0) i3=2*ewaldkmax+1
        i3 = i3-ewaldkmax-1
!!
        i2 = mod(int((i0-1)/(2*ewaldkmax+1))+1,(2*ewaldkmax+1))
        if(i2.eq.0) i2=2*ewaldkmax+1
        i2 = i2-ewaldkmax-1
!!
        i1 = int((i0-1)/((2*ewaldkmax+1)*(2*ewaldkmax+1)))+1
        if(i1.eq.0) i1=2*ewaldkmax+1
        i1 = i1-ewaldkmax-1
!!
!!        write(ounit,*)i0,i1,i2,i3
!!
!!      do i1=-ewaldkmax,ewaldkmax
!!        do i2=-ewaldkmax,ewaldkmax
!!          do i3=-ewaldkmax,ewaldkmax
            kvec(1)=dble(i1)*reclattice(1,1)&
                   +dble(i2)*reclattice(2,1)&
                   +dble(i3)*reclattice(3,1)
            kvec(2)=dble(i1)*reclattice(1,2)&
                   +dble(i2)*reclattice(2,2)&
                   +dble(i3)*reclattice(3,2)
            kvec(3)=dble(i1)*reclattice(1,3)&
                   +dble(i2)*reclattice(2,3)&
                   +dble(i3)*reclattice(3,3)
            ksquared=kvec(1)**2 + kvec(2)**2 + kvec(3)**2 ! in Bohr-2
            factor=(1.d0/ksquared)*&
                   dexp(-1.d0*ksquared/(4.d0*(ewaldalpha)**2))
!!
      if((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0))then
      else
      sumcos=0.0d0
      sumsin=0.0d0
!! loop over all atoms
      do i=1,num_atoms
        kri=kvec(1)*xyzstruct(1,i)+kvec(2)*xyzstruct(2,i)+kvec(3)*xyzstruct(3,i)
        sumcos=sumcos+atomcharge(i)*dcos(kri)
        sumsin=sumsin+atomcharge(i)*dsin(kri)
      enddo !i
      sqsumcos=sumcos*sumcos
      sqsumsin=sumsin*sumsin
      erecip=erecip+factor*(sqsumcos+sqsumsin)
!!
      if(ldoforces)then
        sumsinj=0.0d0
        sumcosj=0.0d0
!! These sums are the same for all derivatives
        do j=1,num_atoms
          krj=kvec(1)*xyzstruct(1,j)+kvec(2)*xyzstruct(2,j)+kvec(3)*xyzstruct(3,j)
          sumcosj=sumcosj+atomcharge(j)*dcos(krj)
          sumsinj=sumsinj+atomcharge(j)*dsin(krj)
        enddo
!!
!! These sums are different for all derivatives
        do i5=1,num_atoms
          do i4=1,3
            sumcosdev=0.0d0
            sumsindev=0.0d0
            do j=1,num_atoms
              krj=kvec(1)*xyzstruct(1,j)+kvec(2)*xyzstruct(2,j)+kvec(3)*xyzstruct(3,j)
!!
              if(i5.eq.j)then
              if(i4.eq.1)then
                temp1=kvec(1)
              elseif(i4.eq.2)then
                temp1=kvec(2)
              elseif(i4.eq.3)then
                temp1=kvec(3)
              endif
              else
                temp1=0.0d0
              endif
!!
              sumcosdev=sumcosdev+dcos(krj)*dchargedxyz(j,i5,i4)&
                       -atomcharge(j)*dsin(krj)*temp1
              sumsindev=sumsindev+dsin(krj)*dchargedxyz(j,i5,i4)&
                       +atomcharge(j)*dcos(krj)*temp1
            enddo
!!
            erecipforce(i4,i5)=erecipforce(i4,i5)&
              -factor*(2.d0*sumcosj*sumcosdev +2.0d0*sumsinj*sumsindev)
!!
          enddo ! i4
        enddo ! i5
!!
      endif ! ldoforces
      endif
!!
!!          enddo ! i3
!!        enddo ! i2
!!      enddo ! i1
!!
      enddo ! i0
!!
      erecip=erecip*twopi/determinant
      erecipforce(:,:)=erecipforce(:,:)*twopi/determinant
!!
!!
      return
      end
