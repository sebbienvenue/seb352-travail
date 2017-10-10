!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: add empirical vdW term based on S. Grimme, J. Comp. Chem. 27, 1787 (2006).

!! FIXME: at the moment vdw energies and forces are simply cut off at cutoffvdw
!! this is not very elegant and should be improved to avoid jumps in energies and forces
!! Probably the best way is to use our symmetry cutoff function
!!
!! called by:
!! - predict.f90
!!
      subroutine addvdw(num_atoms,zelem,&
        xyzstruct,lattice,energyvdw,forcevdw,&
        lperiodic)
!!
      use fileunits
      use globaloptions
      use predictionoptions
!!
      implicit none
!!
      integer i1,i2                                                  ! internal
      integer num_atoms                                              ! in
      integer num_pairs_local                                        ! internal
      integer zelem(max_num_atoms)                                   ! in
      integer lsta(2,max_num_atoms)                                  ! internal
      integer pairs_charge(2,listdim)                                ! internal
      integer pair_lsta(listdim,2)                                   ! internal
      integer pair_lstc(listdim)                                     ! internal
      integer pairs_atom(listdim,2)                                  ! internal   
      integer pair_lste(listdim)                                     ! internal 
!!
      real*8 lstb(listdim,4)                                         ! internal
      real*8 xyzstruct(3,max_num_atoms)                              ! in
      real*8 lattice(3,3)                                            ! in
      real*8 pairs_rr(listdim)                                       ! internal
      real*8 pairs_xyz(listdim,2,3)                                  ! internal
      real*8 pair_lstb(listdim,5)                                    ! internal
      real*8 s6                                                      ! internal
      real*8 c6                                                      ! internal
      real*8 d                                                       ! internal
      real*8 Rr                                                      ! internal
      real*8 energyvdw                                               ! out
      real*8 forcevdw(3,max_num_atoms)                               ! out
      real*8 fdamp                                                   ! internal
      real*8 dfdampdxi                                               ! internal
      real*8 dfdampdyi                                               ! internal
      real*8 dfdampdzi                                               ! internal
      real*8 dfdampdxj                                               ! internal
      real*8 dfdampdyj                                               ! internal
      real*8 dfdampdzj                                               ! internal
      real*8 drijdxi,drijdyi,drijdzi                                 ! internal
      real*8 drijdxj,drijdyj,drijdzj                                 ! internal
      real*8 deltaxj,deltayj,deltazj                                 ! internal
      real*8 temp1                                                   ! internal
      real*8 temp2                                                   ! internal
!!
      logical lperiodic                                              ! in
      logical ldebug_local                                           ! internal
!!
      energyvdw    =0.0d0
      forcevdw(:,:)=0.0d0
!!     
!! FIXME: check units of Grimme parameters, we need Ha and Bohr 
!!
      ldebug_local=.true.
!!
!! get the unique pairs considering periodic boundary conditions
      call neighborpair_para(num_atoms,num_pairs_local,zelem,&
        lsta,lstb,cutoffvdw,lattice,xyzstruct,&
        pairs_rr,pairs_atom,pairs_xyz,pairs_charge,&
        pair_lsta,pair_lstb,pair_lstc,pair_lste,lperiodic)
!!
!! calculate vdw energy for each pair
      do i1=1,num_pairs_local
        s6=vdw_param(pairindex(pairs_charge(1,i1),pairs_charge(2,i1)),1)
        c6=vdw_param(pairindex(pairs_charge(1,i1),pairs_charge(2,i1)),2)
        d =vdw_param(pairindex(pairs_charge(1,i1),pairs_charge(2,i1)),3)
        Rr=vdw_param(pairindex(pairs_charge(1,i1),pairs_charge(2,i1)),4)
!!
        temp1=-dexp(-d*pairs_rr(i1)/Rr -1.0d0)
        fdamp=1.d0/(1.d0+temp1)
        energyvdw=energyvdw -s6*c6*fdamp/(pairs_rr(i1))**6

        if(ldoforces)then
          deltaxj=-1.d0*(pairs_xyz(i1,1,1)-pairs_xyz(i1,2,1))
          deltayj=-1.d0*(pairs_xyz(i1,1,2)-pairs_xyz(i1,2,2))
          deltazj=-1.d0*(pairs_xyz(i1,1,3)-pairs_xyz(i1,2,3))
          drijdxi=-deltaxj/pairs_rr(i1)
          drijdyi=-deltayj/pairs_rr(i1)
          drijdzi=-deltazj/pairs_rr(i1)
          drijdxj=-1.d0*drijdxi
          drijdyj=-1.d0*drijdyi
          drijdzj=-1.d0*drijdzi
!!
          temp2= fdamp**2 * temp1 * d/Rr  
!!
          dfdampdxi= temp2 * drijdxi
          dfdampdyi= temp2 * drijdyi
          dfdampdzi= temp2 * drijdzi
          dfdampdxj= temp2 * drijdxj
          dfdampdyj= temp2 * drijdyj
          dfdampdzj= temp2 * drijdzj
!!
          forcevdw(1,pairs_atom(i1,1))=forcevdw(1,pairs_atom(i1,1))&
            - 6.d0 * s6 * c6 * dfdampdxi * drijdxi / (pairs_rr(i1))**7 
          forcevdw(2,pairs_atom(i1,1))=forcevdw(2,pairs_atom(i1,1))&
            - 6.d0 * s6 * c6 * dfdampdyi * drijdyi / (pairs_rr(i1))**7
          forcevdw(3,pairs_atom(i1,1))=forcevdw(3,pairs_atom(i1,1))&
            - 6.d0 * s6 * c6 * dfdampdzi * drijdzi / (pairs_rr(i1))**7
          forcevdw(1,pairs_atom(i1,2))=forcevdw(1,pairs_atom(i1,2))&
            - 6.d0 * s6 * c6 * dfdampdxj * drijdxj / (pairs_rr(i1))**7
          forcevdw(2,pairs_atom(i1,2))=forcevdw(2,pairs_atom(i1,2))&
            - 6.d0 * s6 * c6 * dfdampdyj * drijdyj / (pairs_rr(i1))**7
          forcevdw(3,pairs_atom(i1,2))=forcevdw(3,pairs_atom(i1,2))&
            - 6.d0 * s6 * c6 * dfdampdzj * drijdzj / (pairs_rr(i1))**7
        endif
!!
        if(ldostress)then
          write(ounit,*)'ERROR: stress is not implemented in addvdw'
          stop
        endif
      enddo
!!
!! output for debugging
      if(ldebug_local)then
        write(ounit,'(a,f14.6)')' addvdw: vdw energy ',energyvdw
        do i1=1,num_atoms
          write(ounit,'(a,i5,3f14.8)')' addvdw: vdw force ',i1,(forcevdw(i2,i1),i2=1,3)
        enddo
      endif
!!
      return
      end
