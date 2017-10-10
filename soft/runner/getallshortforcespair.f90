!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by: 
!! - getshortenergies_parapair.f90
!! - optimize_short_combinedpair.f90
!!
         subroutine getallshortforcespair(ndim,npoints,&
            num_atoms_local,zelem_local,zelemp_local,&
            symfunctionp_local,nnshortforce_local,&
            lattice_local,xyzstruct_local,minvalue_short_pair,maxvalue_short_pair,&
            lperiodic_local)
!!
      use fileunits
      use globaloptions
      use symfunctions
      use nnshort_pair
!!
      implicit none
!!
      integer ndim                                                      ! in
      integer npoints                                                   ! in
      integer num_atoms                                                 ! internal
      integer num_pairs                                                 ! internal
      integer num_atoms_local(ndim)                                    ! in
      integer zelem(max_num_atoms)                                      ! internal
      integer zelemp(2,max_num_pairs)                                   ! internal
      integer zelemp_local(2,ndim,max_num_pairs)                       ! in
      integer zelem_local(ndim,max_num_atoms)
      integer pairs_charge(2,listdim)                                   ! internal
      integer i1                                                        ! internal
      integer ndummy                                                    ! internal
!!
      real*8 dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3)
      real*8 nnshortforce(3,max_num_atoms)                                ! internal 
      real*8 nnshortforce_local(3,max_num_atoms,ndim)                     ! out 

      real*8 depairdsfunc(max_num_pairs,maxnum_funcvalues_short_pair) 
      real*8 symfunctionpdummy(maxnum_funcvalues_short_pair,max_num_pairs)                     ! in
      real*8 symfunctionp_local(maxnum_funcvalues_short_pair,max_num_pairs,ndim)
      real*8 lattice_local(3,3,ndim)                                   ! in
      real*8 xyzstruct_local(3,max_num_atoms,ndim)                     ! in
      real*8 shortforce(3,max_num_pairs)                               ! internal

      real*8 strs_pair(3,3,maxnum_funcvalues_short_pair,max_num_pairs)          

      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)     ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)     ! in
!!
      logical lperiodic_local(ndim)                                    ! in
      logical lperiodic                                                 ! internal
      logical ldummy

!!
      shortforce(:,:)=0.0d0 ! just dummy here
      nnshortforce_local(:,:,:)=0.0d0
      zelemp(:,:)  =0
      zelem(:)     =0
      ndummy = 0
!!
      do i1=1,npoints
        num_atoms               = num_atoms_local(i1)
        zelemp(:,:)             = zelemp_local(:,i1,:)
        zelem(:)                = zelem_local(i1,:)
        lperiodic               = lperiodic_local(i1)
        symfunctionpdummy(:,:)  = symfunctionp_local(:,:,i1)
!!
!! in principle dsfuncdxyz could be precalculated, but needs too much storage on disk
!!
!! get dsfuncdxyz:
!! Caution: symmetry functions cannot be used because here they are not scaled
        call calconefunction_pair(1,&
           ndummy,num_atoms,zelem,&
           num_pairs,pairs_charge,&
           lattice_local(1,1,i1),xyzstruct_local(1,1,i1),symfunctionpdummy,&
           dsfuncdxyz_pair,strs_pair,&
           lperiodic,.true.,ldummy)
!!
!! scale dsfuncdxyz
        if(lscalesym)then
           call scaledsfuncpair(num_pairs,&
             minvalue_short_pair,maxvalue_short_pair,&
             zelemp,dsfuncdxyz_pair,strs_pair)
        endif
!!
        call getshortforcespair(&
          num_atoms,num_pairs,zelemp,&
          symfunctionp_local(1,1,i1),&
          dsfuncdxyz_pair,depairdsfunc,nnshortforce)
!!        
        nnshortforce_local(:,:,i1)=nnshortforce(:,:)
!!
      enddo ! i1
!!
!!
      return
      end
