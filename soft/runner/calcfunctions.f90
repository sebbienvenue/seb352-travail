!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by: -getsymmetryfunctions.f90
!!
!! FIXME: implement effective neighbor lists for pair case to reduce memory
!!
      subroutine calcfunctions(npoints,ndone,&
        iseed,numtrain,numtest,numrej,pointnumber,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        minvalue_ham,maxvalue_ham,avvalue_ham,&
        minvalue_s,maxvalue_s,avvalue_s,&
        minvalue_hexton,maxvalue_hexton,avvalue_hexton,&
        minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
        minvalue_dens,maxvalue_dens,avvalue_dens)
!!
      use fileunits
      use nnflags
      use globaloptions 
      use mode1options
      use symfunctions
      use nnshort_atomic
      use nnshort_pair
      use nnewald
      use nnham
      use basismod
      use structures
!!
      implicit none
!!
      integer npoints                                                     ! in
      integer num_atoms                                                   ! internal
      integer zelem(max_num_atoms)                                        ! internal
      integer iseed                                                       ! in/out
      integer numtrain                                                    ! in/out
      integer numtest                                                     ! in/out
      integer numrej                                                      ! in/out
      integer ndone                                                       ! in
      integer pointnumber                                                 ! in/out
      integer i1,i2,i3                                                    ! internal
      integer istruct                                                     ! internal
      integer num_pairs                                                   ! internal
      integer num_triplets                                                ! internal
      integer, allocatable :: pairs_charge(:,:)                           ! internal
      integer, allocatable :: pairs_charge_list(:,:,:)                    ! internal
!!      integer, allocatable :: triplet_charge(:,:)                         ! internal
!!      integer, allocatable :: triplet_charge_list(:,:,:)                  ! internal 
!! neighbor arrays for short range atomic case
      integer, allocatable :: lsta_shortatomic(:,:)                       ! internal, numbers of neighbors
      integer, allocatable :: lstc_shortatomic(:)                         ! internal, identification of atom
      integer, allocatable :: lste_shortatomic(:)                         ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_short_atomic(:)               ! internal
      integer max_num_neighbors_short_atomic                              ! internal
      integer, allocatable :: neighboridx_short_atomic(:,:)               ! internal
      integer, allocatable :: invneighboridx_short_atomic(:,:)            ! internal
!! neighbor arrays for electrostatic case
      integer, allocatable :: lsta_elec(:,:)                              ! internal, numbers of neighbors
      integer, allocatable :: lstc_elec(:)                                ! internal, identification of atom
      integer, allocatable :: lste_elec(:)                                ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_elec(:)                       ! internal
      integer max_num_neighbors_elec                                      ! internal
      integer, allocatable :: neighboridx_elec(:,:)                       ! internal
      integer, allocatable :: invneighboridx_elec(:,:)                    ! internal
      integer, allocatable :: atomindex_dummy(:)                          ! internal
!! neighbor arrays for Hamiltonian case
      integer, allocatable :: lsta_ham(:,:)                              ! internal, numbers of neighbors
      integer, allocatable :: lstc_ham(:)                                ! internal, identification of atom
      integer, allocatable :: lste_ham(:)                                ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_ham(:)                       ! internal
      integer max_num_neighbors_ham                                      ! internal
      integer, allocatable :: neighboridx_ham(:,:)                       ! internal
      integer, allocatable :: invneighboridx_ham(:,:)                    ! internal
!!      integer, allocatable :: atomindex_dummy(:)                          ! internal
!! neighbor arrays for Hamiltonian overlap case
      integer, allocatable :: lsta_s(:,:)                              ! internal, numbers of neighbors
      integer, allocatable :: lstc_s(:)                                ! internal, identification of atom
      integer, allocatable :: lste_s(:)                                ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_s(:)                       ! internal
      integer max_num_neighbors_s                                      ! internal
      integer, allocatable :: neighboridx_s(:,:)                       ! internal
      integer, allocatable :: invneighboridx_s(:,:)                    ! internal
!!      integer, allocatable :: atomindex_dummy(:)                          ! internal
!! neighbor arrays for Hamiltonian Hexton case
      integer, allocatable :: lsta_hexton(:,:)                              ! internal, numbers of neighbors
      integer, allocatable :: lstc_hexton(:)                                ! internal, identification of atom
      integer, allocatable :: lste_hexton(:)                                ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_hexton(:)                       ! internal
      integer max_num_neighbors_hexton                                      ! internal
      integer, allocatable :: neighboridx_hexton(:,:)                       ! internal
      integer, allocatable :: invneighboridx_hexton(:,:)                    ! internal
!!      integer, allocatable :: atomindex_dummy(:)                          ! internal
!! neighbor arrays for Hamiltonian Hextoff case
      integer, allocatable :: lsta_hextoff(:,:)                              ! internal, numbers of neighbors
      integer, allocatable :: lstc_hextoff(:)                                ! internal, identification of atom
      integer, allocatable :: lste_hextoff(:)                                ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_hextoff(:)                       ! internal
      integer max_num_neighbors_hextoff                                      ! internal
      integer, allocatable :: neighboridx_hextoff(:,:)                       ! internal
      integer, allocatable :: invneighboridx_hextoff(:,:)                    ! internal
!!      integer, allocatable :: atomindex_dummy(:)                          ! internal
!! neighbor arrays for Hamiltonian density case
      integer, allocatable :: lsta_dens(:,:)                              ! internal, numbers of neighbors
      integer, allocatable :: lstc_dens(:)                                ! internal, identification of atom
      integer, allocatable :: lste_dens(:)                                ! internal, nuclear charge of atom
      integer, allocatable :: num_neighbors_dens(:)                       ! internal
      integer max_num_neighbors_dens                                      ! internal
      integer, allocatable :: neighboridx_dens(:,:)                       ! internal
      integer, allocatable :: invneighboridx_dens(:,:)                    ! internal
!!      integer, allocatable :: atomindex_dummy(:)                          ! internal

!!
      real*8 lattice(3,3)                                                 ! internal
      real*8 totalcharge                                                  ! internal
      real*8 totalenergy                                                  ! internal
      real*8 elecenergy                                                   ! internal
      real*8 nntbenergy                                                   ! internal
      real*8 nntb_s_energy                                                ! internal
      real*8 nntb_hexton_energy                                           ! internal
      real*8 nntb_hextoff_energy                                          ! internal
      real*8 nntb_dens_energy                                             ! internal
      real*8 atomcharge(max_num_atoms)                                    ! internal
      real*8 atomenergy(max_num_atoms)                                    ! internal
      real*8 totalforce(3,max_num_atoms)                                  ! internal
      real*8 elecforce(3,max_num_atoms)                                   ! internal
      real*8 nntbforce(3,max_num_atoms)                                   ! internal
      real*8 nntb_s_force(3,max_num_atoms)                                ! internal
      real*8 nntb_hexton_force(3,max_num_atoms)                           ! internal
      real*8 nntb_hextoff_force(3,max_num_atoms)                          ! internal
      real*8 nntb_dens_force(3,max_num_atoms)                             ! internal
      real*8 ran0                                                         ! internal
      real*8 random                                                       ! internal
      real*8 symfunction(maxnum_funcvalues_short_atomic,max_num_atoms)    ! internal
      real*8 symfunctionp(maxnum_funcvalues_short_pair,max_num_pairs)     ! internal
      real*8 symfunctione(maxnum_funcvalues_elec,max_num_atoms)           ! internal
      real*8 symfunctionework(maxnum_funcvalues_elec,max_num_atoms)       ! internal
      real*8 symfunction_ham(maxnum_funcvalues_ham,max_num_atoms)         ! internal
      real*8 symfunction_s(maxnum_funcvalues_s,max_num_atoms)             ! internal
      real*8 symfunction_hexton(maxnum_funcvalues_hexton,max_num_atoms)   ! internal
      real*8 symfunction_hextoff(maxnum_funcvalues_hextoff)               ! internal
      real*8 symfunction_dens(maxnum_funcvalues_dens,max_num_atoms)       ! internal
      real*8 strs(3,3,maxnum_funcvalues_short_atomic,max_num_atoms)       ! internal
      real*8 strs_pair(3,3,maxnum_funcvalues_short_pair,max_num_pairs)    ! internal
      real*8 strse(3,3,maxnum_funcvalues_elec,max_num_atoms)              ! internal
      real*8 strs_ham(3,3,maxnum_funcvalues_ham,max_num_atoms)            ! internal
      real*8 strs_s(3,3,maxnum_funcvalues_s,max_num_atoms)               ! internal
      real*8 strs_hexton(3,3,maxnum_funcvalues_hexton,max_num_atoms)      ! internal
      real*8 strs_hextoff(3,3,maxnum_funcvalues_hextoff,ntriplets)        ! internal
      real*8 strs_dens(3,3,maxnum_funcvalues_dens,max_num_atoms)          ! internal
      real*8, allocatable :: dsfuncdxyz(:,:,:,:)                          ! internal
      real*8, allocatable :: dsfuncdxyz_pair(:,:,:,:)                     ! internal
      real*8, allocatable :: dsfuncdxyze(:,:,:,:)                         ! internal
      real*8, allocatable :: dchargedxyz(:,:,:)                           ! internal
      real*8, allocatable :: dsfuncdxyz_ham(:,:,:,:)                      ! internal
      real*8, allocatable :: dsfuncdxyz_s(:,:,:,:)                        ! internal
      real*8, allocatable :: dsfuncdxyz_hexton(:,:,:,:)                   ! internal
      real*8, allocatable :: dsfuncdxyz_hextoff(:,:,:,:)                  ! internal
      real*8, allocatable :: dsfuncdxyz_dens(:,:,:,:)                     ! internal
      real*8 forcesum(3)                                                  ! internal
      real*8 absforcesum                                                  ! internal
      real*8 highfcomponent                                               ! internal

      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)                  ! in 
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)                  ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)                   ! in
      real*8 minvalue_ham(nelem,maxnum_funcvalues_ham)                    ! in 
      real*8 maxvalue_ham(nelem,maxnum_funcvalues_ham)                    ! in
      real*8 minvalue_s(nelem,maxnum_funcvalues_s)                        ! in
      real*8 maxvalue_s(nelem,maxnum_funcvalues_s)                        ! in
      real*8 minvalue_hexton(nelem,maxnum_funcvalues_hexton)              ! in
      real*8 maxvalue_hexton(nelem,maxnum_funcvalues_hexton)              ! in
      real*8 minvalue_hextoff(1,maxnum_funcvalues_hextoff)            ! in
      real*8 maxvalue_hextoff(1,maxnum_funcvalues_hextoff)            ! in
      real*8 minvalue_dens(nelem,maxnum_funcvalues_dens)                  ! in
      real*8 maxvalue_dens(nelem,maxnum_funcvalues_dens)                  ! in

      real*8 avvalue_ham(nelem,maxnum_funcvalues_ham)                     ! in
      real*8 avvalue_s(nelem,maxnum_funcvalues_s)                         ! in
      real*8 avvalue_hexton(nelem,maxnum_funcvalues_hexton)               ! in
      real*8 avvalue_hextoff(1,maxnum_funcvalues_hextoff)             ! in
      real*8 avvalue_dens(nelem,maxnum_funcvalues_dens)                   ! in

      real*8 nnatomcharge(max_num_atoms)                                  ! internal 
      real*8, allocatable :: lstb_shortatomic(:,:)                        ! xyz and r_ij 
      real*8, allocatable :: lstb_elec(:,:)                               ! xyz and r_ij 
      real*8, allocatable :: lstb_ham(:, :)                               ! xyz and r_ij
      real*8, allocatable :: lstb_s(:,:)                                  ! xyz and r_ij
      real*8, allocatable :: lstb_hexton(:,:)                             ! xyz and r_ij
      real*8, allocatable :: lstb_hextoff(:,:)                            ! xyz and r_ij
      real*8, allocatable :: lstb_dens(:,:)                               ! xyz and r_ij

!!
      character*2 elementsymbol(max_num_atoms)                            ! internal
!!
      logical lperiodic                                                   ! internal
      logical ldoforces                                                   ! internal here!!!
      logical lrmin(npoints)                                              ! internal
      logical lforceok                                                    ! internal
!!
!!=============================================================================
!! initializations
!!=============================================================================
      elecforce(:,:)        = 0.0d0
      elecforce_list(:,:,:) = 0.0d0
      elecenergy            = 0.0d0
      nntbenergy            = 0.0d0
      nntbforce(:,:)        = 0.0d0
      elecenergy_list(:)    = 0.0d0
      lrmin(:)              = .true.
      ldoforces             = luseforces
      symfunction(:,:)      = 0.0d0
      symfunctionp(:,:)     = 0.0d0
      symfunctione(:,:)     = 0.0d0
      symfunction_ham(:,:)  = 0.0d0
      symfunction_s(:,:)    = 0.0d0
      symfunction_hextoff(:)       = 0.0d0
      symfunction_hexton(:,:)      = 0.0d0
      symfunction_dens(:,:)        = 0.0d0
      lforceok              = .true.

!! caution: pairs_charge arrays must not be allocated for atomic case, because max_num_pairs is not available => segmentation fault
      if((nn_type_short.eq.2).or.(lnntb))then
        allocate(pairs_charge_list(2,listdim,nblock)) 
        allocate(pairs_charge(2,listdim))
      endif
!! caution: triplet_charge arrays must not be allocated for atomic case, because max_num_triplets is not available => segmentation fault
      if((lnntb).and.(nntb_flag(3)))then             
!!        allocate(triplet_charge_list(3,listdim,nblock))
!!        allocate(triplet_charge(3,listdim))
      endif
!!
!!=============================================================================
!! prepare input for calculating the symmetry functions of one structure
!!=============================================================================
      do i1=1,npoints
        num_atoms        = num_atoms_list(i1)
        zelem(:)         = zelem_list(i1,:) 
        lattice(:,:)     = lattice_list(:,:,i1)
        totalcharge      = totalcharge_list(i1)
        totalenergy      = totalenergy_list(i1)
        atomcharge(:)    = atomcharge_list(i1,:)
        atomenergy(:)    = atomenergy_list(i1,:)
        totalforce(:,:)  = totalforce_list(:,:,i1)
        elementsymbol(:) = elementsymbol_list(i1,:)
        lperiodic        = lperiodic_list(i1)
        istruct          = ndone+i1

!!
!!=============================================================================
!! allocate arrays for neighbor lists
!!=============================================================================
        if(lshort.and.(nn_type_short.eq.1))then
          allocate(lsta_shortatomic(2,max_num_atoms))
          allocate(lstc_shortatomic(listdim))
          allocate(lste_shortatomic(listdim))
          allocate(lstb_shortatomic(listdim,4))
          allocate(num_neighbors_short_atomic(num_atoms))
        endif
        if(lelec.and.(nn_type_elec.eq.1))then
          allocate(lsta_elec(2,max_num_atoms))
          allocate(lstc_elec(listdim))
          allocate(lste_elec(listdim))
          allocate(lstb_elec(listdim,4))
          allocate(num_neighbors_elec(num_atoms))
        endif
        if(lnntb)then
          if(nntb_flag(0))then
            allocate(lsta_ham(2,max_num_atoms))
            allocate(lstc_ham(listdim))
            allocate(lste_ham(listdim))
            allocate(lstb_ham(listdim,4))
            allocate(num_neighbors_ham(num_atoms))
          endif
          if(nntb_flag(1))then
            allocate(lsta_s(2,max_num_atoms))
            allocate(lstc_s(listdim))
            allocate(lste_s(listdim))
            allocate(lstb_s(listdim,4))
            allocate(num_neighbors_s(num_atoms))
          endif
          if(nntb_flag(2))then
            allocate(lsta_hexton(2,max_num_atoms))
            allocate(lstc_hexton(listdim))
            allocate(lste_hexton(listdim))
            allocate(lstb_hexton(listdim,4))
            allocate(num_neighbors_hexton(num_atoms))
          endif
          if(nntb_flag(3))then
            allocate(lsta_hextoff(2,max_num_atoms))
            allocate(lstc_hextoff(listdim))
            allocate(lste_hextoff(listdim))
            allocate(lstb_hextoff(listdim,4))
            allocate(num_neighbors_hextoff(num_atoms))
          endif
          if(nntb_flag(4))then
            allocate(lsta_dens(2,max_num_atoms))
            allocate(lstc_dens(listdim))
            allocate(lste_dens(listdim))
            allocate(lstb_dens(listdim,4))
            allocate(num_neighbors_dens(num_atoms))
          endif

        endif
!!
!!=============================================================================
!! determine neighbor lists of all atoms (lsta,lstb,lstc,lste), max_num_neighbors and num_neighbors
!!=============================================================================
        if(lshort.and.(nn_type_short.eq.1))then
          call getneighborsatomic(&
            num_atoms,num_neighbors_short_atomic,zelem,&
            max_num_neighbors_short_atomic,&
            lsta_shortatomic,lstc_shortatomic,lste_shortatomic,&
            maxcutoff_short_atomic,lattice,xyzstruct_list(1,1,i1),&
            lstb_shortatomic,lperiodic)
        endif
        if(lelec.and.(nn_type_elec.eq.1))then
          call getneighborsatomic(&
            num_atoms,num_neighbors_elec,zelem,&
            max_num_neighbors_elec,&
            lsta_elec,lstc_elec,lste_elec,&
            maxcutoff_elec,lattice,xyzstruct_list(1,1,i1),&
            lstb_elec,lperiodic)
        endif
        if(lnntb)then
          if(nntb_flag(0))then
            call getneighborsatomic(&
             num_atoms,num_neighbors_ham,zelem,&
             max_num_neighbors_ham,&
             lsta_ham,lstc_ham,lste_ham,&
             maxcutoff_ham,lattice,xyzstruct_list(1,1,i1),&
             lstb_ham,lperiodic)
          endif
          if(nntb_flag(1))then
            call getneighborsatomic(&
             num_atoms,num_neighbors_s,zelem,&
             max_num_neighbors_s,&
             lsta_s,lstc_s,lste_s,&
             maxcutoff_s,lattice,xyzstruct_list(1,1,i1),&
             lstb_s,lperiodic)
          endif
          if(nntb_flag(2))then
            call getneighborsatomic(&
             num_atoms,num_neighbors_hexton,zelem,&
             max_num_neighbors_hexton,&
             lsta_hexton,lstc_hexton,lste_hexton,&
             maxcutoff_hexton,lattice,xyzstruct_list(1,1,i1),&
             lstb_hexton,lperiodic)
          endif
          if(nntb_flag(3))then
            call getneighborsatomic(&
             num_atoms,num_neighbors_hextoff,zelem,&
             max_num_neighbors_hextoff,&
             lsta_hextoff,lstc_hextoff,lste_hextoff,&
             maxcutoff_hextoff,lattice,xyzstruct_list(1,1,i1),&
             lstb_hextoff,lperiodic)
          endif
          if(nntb_flag(4))then
            call getneighborsatomic(&
             num_atoms,num_neighbors_dens,zelem,&
             max_num_neighbors_dens,&
             lsta_dens,lstc_dens,lste_dens,&
             maxcutoff_dens,lattice,xyzstruct_list(1,1,i1),&
             lstb_dens,lperiodic)
          endif

        endif

!!
!!=============================================================================
!! allocate further arrays and determine neighboridx arrays 
!!=============================================================================
        if(lshort.and.(nn_type_short.eq.1))then
          allocate(dsfuncdxyz(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_short_atomic,3))
          allocate(neighboridx_short_atomic(num_atoms,0:max_num_neighbors_short_atomic))  
          allocate(invneighboridx_short_atomic(num_atoms,max_num_atoms))  
          call getneighboridxatomic(num_atoms,listdim,&
            max_num_atoms,max_num_neighbors_short_atomic,&
            lsta_shortatomic,lstc_shortatomic,neighboridx_short_atomic,&
            invneighboridx_short_atomic)
        elseif(lshort.and.(nn_type_short.eq.2))then
          allocate(dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3)) 
        endif
        if(lelec.and.(nn_type_elec.eq.1))then
          allocate(dsfuncdxyze(maxnum_funcvalues_elec,max_num_atoms,0:max_num_neighbors_elec,3))
          allocate(neighboridx_elec(num_atoms,0:max_num_neighbors_elec))  
          allocate(invneighboridx_elec(num_atoms,max_num_atoms))  
          call getneighboridxatomic(num_atoms,listdim,&
            max_num_atoms,max_num_neighbors_elec,&
            lsta_elec,lstc_elec,neighboridx_elec,invneighboridx_elec)
        endif
!!
        if(lnntb)then
!!          if(maxnum_layers_ham.ge.1)then
!!            allocate(dsfuncdxyz_ham(maxnum_funcvalues_ham,max_num_atoms,0:max_num_neighbors_ham,3))
!!            allocate(neighboridx_ham(num_atoms,0:max_num_neighbors_ham))
!!            allocate(invneighboridx_ham(num_atoms,max_num_atoms))
!!            call getneighboridxatomic(num_atoms,listdim,&
!!              max_num_atoms,max_num_neighbors_ham,&
!!              lsta_ham,lstc_ham,neighboridx_ham,invneighboridx_ham)
!!          endif
!!           if(nntb_flag(1))then
!!            allocate(dsfuncdxyz_s(maxnum_funcvalues_s,max_num_atoms,0:max_num_neighbors_s,3))
!!            allocate(neighboridx_s(num_atoms,0:max_num_neighbors_s))
!!            allocate(invneighboridx_s(num_atoms,max_num_atoms))
!!            call getneighboridxatomic(num_atoms,listdim,&
!!              max_num_atoms,max_num_neighbors_s,&
!!              lsta_s,lstc_s,neighboridx_s,invneighboridx_s)
!!          endif
!!           if(nntb_flag(2))then
!!            allocate(dsfuncdxyz_hexton(maxnum_funcvalues_hexton,max_num_atoms,0:max_num_neighbors_hexton,3))
!!            allocate(neighboridx_hexton(num_atoms,0:max_num_neighbors_hexton))
!!            allocate(invneighboridx_hexton(num_atoms,max_num_atoms))
!!            call getneighboridxatomic(num_atoms,listdim,&
!!              max_num_atoms,max_num_neighbors_hexton,&
!!              lsta_hexton,lstc_hexton,neighboridx_hexton,invneighboridx_hexton)
!!          endif
           if(nntb_flag(3))then
            allocate(dsfuncdxyz_hextoff(maxnum_funcvalues_hextoff,max_num_atoms,0:max_num_triplets,3))
            allocate(neighboridx_hextoff(num_atoms,0:max_num_neighbors_hextoff))
            allocate(invneighboridx_hextoff(num_atoms,max_num_atoms))
            call getneighboridxatomic(num_atoms,listdim,&
              max_num_atoms,max_num_neighbors_hextoff,&
              lsta_hextoff,lstc_hextoff,neighboridx_hextoff,invneighboridx_hextoff)
!!            write(*,*) 'getneighbouridxatomic',num_atoms,listdim
!!            write(*,*) invneighboridx_hextoff(1,:)
!!            write(*,*) invneighboridx_hextoff(2,:)
!!            write(*,*) invneighboridx_hextoff(3,:)
!!            write(*,*) neighboridx_hextoff(1,:)
!!            write(*,*) neighboridx_hextoff(2,:)
!!            write(*,*) neighboridx_hextoff(3,:)           
            
          endif
!!           if(nntb_flag(4))then
!!            allocate(dsfuncdxyz_dens(maxnum_funcvalues_dens,max_num_atoms,0:max_num_neighbors_dens,3))
!!            allocate(neighboridx_dens(num_atoms,0:max_num_neighbors_dens))
!!            allocate(invneighboridx_dens(num_atoms,max_num_atoms))
!!            call getneighboridxatomic(num_atoms,listdim,&
!!              max_num_atoms,max_num_neighbors_dens,&
!!              lsta_dens,lstc_dens,neighboridx_dens,invneighboridx_dens)
!!          endif
!!
        endif
!!
        allocate(atomindex_dummy(num_atoms))
        do i2=1,num_atoms
          atomindex_dummy(i2)=i2
        enddo
!!=============================================================================
!! calculate the short range symmetry functions for one structure here
!!=============================================================================
        if(lshort.and.(nn_type_short.eq.1))then
          call calconefunction_atomic(cutoff_type,max_num_neighbors_short_atomic,&
            max_num_atoms,1,num_atoms,atomindex_dummy,num_atoms,elementindex,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic, &
            nelem,zelem,listdim,&
            lsta_shortatomic,lstc_shortatomic,lste_shortatomic,&
            invneighboridx_short_atomic,&
            function_type_short_atomic,symelement_short_atomic,&
            xyzstruct_list(1,1,i1),symfunction,rmin,&
            funccutoff_short_atomic,eta_short_atomic,rshift_short_atomic,&
            lambda_short_atomic,zeta_short_atomic,dsfuncdxyz,strs,lstb_shortatomic,&
            lperiodic,ldoforces,ldostress,lrmin(i1))
        elseif(lshort.and.(nn_type_short.eq.2))then
          call calconefunction_pair(1,&
            istruct,num_atoms,zelem,&
            num_pairs,pairs_charge,&
            lattice,xyzstruct_list(1,1,i1),symfunctionp,&
            dsfuncdxyz_pair,strs_pair,&
            lperiodic,ldoforces,lrmin(i1))
!!
            pairs_charge_list(:,:,i1)=pairs_charge(:,:)
        endif
!!
!!=============================================================================
!! calculate the electrostatic symmetry functions for one structure here 
!!=============================================================================
        if(lelec.and.(nn_type_elec.eq.1))then
          call calconefunction_atomic(cutoff_type,max_num_neighbors_elec,&
            max_num_atoms,1,num_atoms,atomindex_dummy,num_atoms,elementindex,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            nelem,zelem,listdim,&
            lsta_elec,lstc_elec,lste_elec,invneighboridx_elec,&
            function_type_elec,symelement_elec,&
            xyzstruct_list(1,1,i1),symfunctione,rmin,&
            funccutoff_elec,eta_elec,rshift_elec,lambda_elec,zeta_elec,dsfuncdxyze,strse,lstb_elec,&
            lperiodic,ldoforces,ldostress,lrmin(i1))
        endif

!!=============================================================================
!! calculate the hamiltonian symmetry functions for one structure here
!!=============================================================================
        if(lnntb)then
          if(nntb_flag(0))then
          call calconefunction_pair(1,&
            istruct,num_atoms,zelem,&
            num_pairs,pairs_charge,&
            lattice,xyzstruct_list(1,1,i1),symfunctionp,&
            dsfuncdxyz_pair,strs_pair,&
            lperiodic,ldoforces,lrmin(i1))
!!
            pairs_charge_list(:,:,i1)=pairs_charge(:,:)
          endif
          if(nntb_flag(1))then
          call calconefunction_pair(1,&
            istruct,num_atoms,zelem,&
            num_pairs,pairs_charge,&
            lattice,xyzstruct_list(1,1,i1),symfunctionp,&
            dsfuncdxyz_pair,strs_pair,&
            lperiodic,ldoforces,lrmin(i1))
!!
            pairs_charge_list(:,:,i1)=pairs_charge(:,:)
          endif
          if(nntb_flag(2))then
          call calconefunction_pair(1,&
            istruct,num_atoms,zelem,&
            num_pairs,pairs_charge,&
            lattice,xyzstruct_list(1,1,i1),symfunctionp,&
            dsfuncdxyz_pair,strs_pair,&
            lperiodic,ldoforces,lrmin(i1))
!!
            pairs_charge_list(:,:,i1)=pairs_charge(:,:)
          endif
          if(nntb_flag(3))then
          call calconefunction_hextoff_mode1(1,&
            istruct,num_atoms,zelem,&
            xyzstruct_list(1,1,i1),symfunction_hextoff)
!!
!!            write(*,*) 'Check output'
!!            write(*,*) symfunction_hextoff
          endif
          if(nntb_flag(4))then
          call calconefunction_pair(1,&
            istruct,num_atoms,zelem,&
            num_pairs,pairs_charge,&
            lattice,xyzstruct_list(1,1,i1),symfunctionp,&
            dsfuncdxyz_pair,strs_pair,&
            lperiodic,ldoforces,lrmin(i1))
!!
            pairs_charge_list(:,:,i1)=pairs_charge(:,:)          
          endif

        endif
!!
!!=============================================================================
!! deallocate arrays 
!!=============================================================================
        deallocate(atomindex_dummy)
        if(lshort)then
          if(nn_type_short.eq.1)then
            deallocate(lsta_shortatomic)
            deallocate(lstc_shortatomic)
            deallocate(lste_shortatomic)
            deallocate(lstb_shortatomic)
            deallocate(invneighboridx_short_atomic)  
          endif
          if(nn_type_short.eq.1)then
            deallocate(dsfuncdxyz)
          elseif(nn_type_short.eq.2)then
            deallocate(dsfuncdxyz_pair)
          endif
        endif !lshort
        if(lelec.and.(nn_type_elec.eq.1))then
          deallocate(lsta_elec)
          deallocate(lstc_elec)
          deallocate(lste_elec)
          deallocate(lstb_elec)
          deallocate(invneighboridx_elec)  
        endif
        if(lnntb)then
          if(nntb_flag(0))then
            deallocate(lsta_ham)
            deallocate(lstc_ham)
            deallocate(lste_ham)
            deallocate(lstb_ham)
            deallocate(num_neighbors_ham)
          endif
          if(nntb_flag(1))then
            deallocate(lsta_s)
            deallocate(lstc_s)
            deallocate(lste_s)
            deallocate(lstb_s)
            deallocate(num_neighbors_s)
          endif
          if(nntb_flag(2))then
            deallocate(lsta_hexton)
            deallocate(lstc_hexton)
            deallocate(lste_hexton)
            deallocate(lstb_hexton)
            deallocate(num_neighbors_hexton)
          endif
          if(nntb_flag(3))then
            deallocate(lsta_hextoff)
            deallocate(lstc_hextoff)
            deallocate(lste_hextoff)
            deallocate(lstb_hextoff)
            deallocate(num_neighbors_hextoff)
            deallocate(dsfuncdxyz_hextoff)
            deallocate(neighboridx_hextoff)
            deallocate(invneighboridx_hextoff)
          endif
          if(nntb_flag(4))then
            deallocate(lsta_dens)
            deallocate(lstc_dens)
            deallocate(lste_dens)
            deallocate(lstb_dens)
            deallocate(num_neighbors_dens)
          endif
        endif
!!
!!=============================================================================
!! first initialize shortforce_list as totalforce_list 
!!=============================================================================
        shortforce_list(:,:,i1) =totalforce_list(:,:,i1)
!!
!!=============================================================================
!! calculate the electrostatic contribution to the total energy and forces
!!=============================================================================
        if(lelec)then
          allocate(dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)) 
          dchargedxyz(:,:,:)=0.0d0
!!=============================================================================
!! fixed charge case: 
!!=============================================================================
          if(nn_type_elec.eq.3)then 
            if(lperiodic) then
              call getewaldenergy(max_num_neighbors_elec,&
                neighboridx_elec,num_neighbors_elec,num_atoms,zelem,&
                lattice,xyzstruct_list(1,1,i1),atomcharge,elecenergy,&
                dchargedxyz,elecforce,.false.)
            else ! not periodic
!! calculate electrostatic energy for non-periodic system and fixed charges
              call electrostatic(num_atoms,&
                atomcharge,xyzstruct_list(1,1,i1),elecenergy)
!! calculate electrostatic force for non-periodic system and fixed charges
              if(ldoforces)then
                call splitcoulombforces(&
                  num_atoms,atomcharge,xyzstruct_list(1,1,i1),elecforce)
              endif ! ldoforces
            endif ! lperiodic
!!=============================================================================
!! environment-dependent charges:
!! We need to get the charges and forces from a prepared NN fit here!
!!=============================================================================
          else
!!
!!=============================================================================
!! Step 1: prepare symmetry functions
!!=============================================================================
!! for the scaling we have to make a working copy of symfuncione -> symfunctionework
!! (symfunctione itself must be written to file unscaled!)
!! for dsfuncdxyze this is not necessary, because it is not written and can be modified
            symfunctionework(:,:)=symfunctione(:,:)
!!=============================================================================
!! scale the symmetry functions for the charge prediction
!!=============================================================================
            call scalesymone(nelem,&
              maxnum_funcvalues_elec,num_funcvalues_elec,num_atoms,&
              zelem,symfunctionework,&
              minvalue_elec,maxvalue_elec,avvalue_elec,&
              scmin_elec,scmax_elec)
!!=============================================================================
!! we also need to scale the derivative terms dsfuncdxyze and strse 
!!=============================================================================
            call scaledsfunc(max_num_neighbors_elec,&
              maxnum_funcvalues_elec,num_funcvalues_elec,&
              nelem,num_atoms,minvalue_elec,maxvalue_elec,&
              scmin_elec,scmax_elec,&
              zelem,dsfuncdxyze,strse)
!!=============================================================================
!! predict the atomic charges 'nnatomcharge' for this structure
!!=============================================================================
            nnatomcharge(:)=0.0d0
            call calconecharge(num_atoms,&
              zelem,symfunctionework,nnatomcharge)
!!=============================================================================
!! Step 2: calculate dchargedxyz array for forces
!!=============================================================================
            if(lperiodic)then
              call getdchargedxyz(max_num_neighbors_elec,&
                num_neighbors_elec,neighboridx_elec,num_atoms,zelem,&
                dsfuncdxyze,dchargedxyz,symfunctionework)
            else ! not periodic
              call getcoulombdchargedxyz(max_num_neighbors_elec,&
                num_neighbors_elec,num_atoms,zelem,&
                dsfuncdxyze,symfunctionework,dchargedxyz)
            endif ! lperiodic
!!
!! FIXME: Dirty workaround in case we are intending to do the first charge fit without short range part
!! (if we want to fit short range part then NN electrostatic E and F should be removed, but otherwise not)
            if(.not.lshort)then
              nnatomcharge(:)=atomcharge(:)
            endif
!!
!!=============================================================================
!! Step 3: calculate the NN electrostatic energy and forces
!!=============================================================================
            if(lperiodic) then
              call getewaldenergy(max_num_neighbors_elec,&
                neighboridx_elec,num_neighbors_elec,num_atoms,zelem,&
                lattice,xyzstruct_list(1,1,i1),nnatomcharge,elecenergy,&
                dchargedxyz,elecforce,.true.)
            else ! not periodic
!!=============================================================================
!! calculate ewald energy for non-periodic system
!!=============================================================================
              call electrostatic(num_atoms,&
                nnatomcharge,xyzstruct_list(1,1,i1),elecenergy)
!!=============================================================================
!! calculate ewald force for non-periodic system 
!!=============================================================================
              if(ldoforces)then
                call getcoulombforcesone(max_num_neighbors_elec,&
                  num_atoms,dchargedxyz,&
                  nnatomcharge,xyzstruct_list(1,1,i1),elecforce)
              endif ! ldoforces
            endif ! lperiodic
          endif ! nn_type_elec.eq.3 
!!
          elecforce_list(:,:,i1)  =elecforce(:,:)
          shortforce_list(:,:,i1) =shortforce_list(:,:,i1)-elecforce_list(:,:,i1)
          deallocate(dchargedxyz) 
!!
        else ! no electrostatics
          elecenergy=0.0d0
        endif ! lelec
!!




        if(lnntb)then

!!          nntbforce_list(:,:,i1)  =nntbforce(:,:)
!!          shortforce_list(:,:,i1) =shortforce_list(:,:,i1)-nntbforce_list(:,:,i1)
        else ! no nntb 
          nntbenergy=0.0d0
        endif ! lelec








        if(lshort)then
          if(nn_type_short.eq.1)then
            symfunction_short_atomic_list(:,:,i1) = symfunction(:,:)
            deallocate(num_neighbors_short_atomic)
            deallocate(neighboridx_short_atomic)
          elseif(nn_type_short.eq.2)then
            symfunction_short_pair_list(:,:,i1) = symfunctionp(:,:)
          endif
        endif ! lshort
        if(lelec.and.(nn_type_elec.eq.1))then
          symfunction_elec_list(:,:,i1)= symfunctione(:,:)
          elecenergy_list(i1)          = elecenergy
          shortenergy_list(i1)         = totalenergy_list(i1)-elecenergy_list(i1)
          deallocate(num_neighbors_elec)
          deallocate(neighboridx_elec)
          deallocate(dsfuncdxyze)
        else
          shortenergy_list(i1)     = totalenergy_list(i1)
        endif
        if(lnntb) then
          if(nntb_flag(3)) then
            symfunction_hextoff_list(:,i1) = symfunction_hextoff
          endif
        endif
      enddo ! i1 over npoints
!!
!!=============================================================================
!!=============================================================================
!! write points of this block of structures to files
!!=============================================================================
!!=============================================================================
      do i1=1,npoints
        pointnumber=pointnumber+1

!!
!!=============================================================================
!! check if the sum force vector is close to zero
!!=============================================================================
        if(lcheckinputforces)then
          forcesum(:)=0.0d0
          do i2=1,num_atoms_list(i1)
            forcesum(1)=forcesum(1)+totalforce_list(1,i2,i1)
            forcesum(2)=forcesum(2)+totalforce_list(2,i2,i1)
            forcesum(3)=forcesum(3)+totalforce_list(3,i2,i1)
          enddo 
          absforcesum=dsqrt(forcesum(1)**2 + forcesum(2)**2 + forcesum(3)**2)
!!          absforcesum=absforcesum/dble(num_atoms) ! should we do this to be independent of system size???
          if(absforcesum.gt.inputforcethreshold)then
            write(ounit,'(a,i10,f14.6,4x,3f14.6)')&
              'WARNING: net force vector is large for structure ',&
              pointnumber,absforcesum,forcesum(1),forcesum(2),forcesum(3)
          endif
        endif

!!=============================================================================
!! decide if the point should be used (no high energies/forces and no too close atoms)
!!=============================================================================
        if(lfitfthres)then
          lforceok=.true.     ! by default use structure
          do i2=1,num_atoms_list(i1)
            do i3=1,3
              if(abs(totalforce_list(i3,i2,i1)).gt.fitfthres)then
                lforceok=.false. ! one too large force component found
                highfcomponent=totalforce_list(i3,i2,i1)
              endif
            enddo
          enddo
        endif

        if(((lfitethres.and.(totalenergy_list(i1).lt.(num_atoms_list(i1)*fitethres)))&
           .or.(.not.lfitethres)).and.lrmin(i1).and.&
           ((lfitfthres.and.lforceok).or.(.not.lfitfthres)))then
!!
!!=============================================================================
!! get random number for splitting in training and test set
!!=============================================================================
          random=ran0(iseed)
!!
!!=============================================================================
!! check if write format statements are sufficient for the number of short range symmetry functions
!!=============================================================================
          if(lshort)then
            if(nn_type_short.eq.1)then
              do i2=1,nelem
                if(num_funcvalues_short_atomic(i2).gt.2500)then
                  write(ounit,*)'Error: only 500 funcvalues possible'
                  stop
                endif
              enddo
            elseif(nn_type_short.eq.2)then
              do i2=1,npairs
                if(num_funcvalues_short_pair(i2).gt.2500)then
                  write(ounit,*)'Error: only 500 funcvalues possible'
                  stop
                endif
              enddo
            endif
          endif ! lshort
          if(lelec.and.(nn_type_elec.eq.1))then
            do i2=1,nelem
              if(num_funcvalues_elec(i2).gt.2500)then
                write(ounit,*)'Error: only 500 funcvaluese possible'
                stop
              endif
            enddo
          endif
!!
!!=============================================================================
!! decide if point is for training or test set
!!=============================================================================
          if(random.gt.splitthres)then ! this is a training point
!!
!!=============================================================================
!! write function.data
!!=============================================================================
            if(lshort)then
              if(nn_type_short.eq.1)then
                call writeatomicsymfunctions(symunit,i1,&
                  maxnum_funcvalues_short_atomic,&
                  num_funcvalues_short_atomic,symfunction_short_atomic_list)
              elseif(nn_type_short.eq.2)then
                call writepairsymfunctions(symunit,i1,pairs_charge_list,&
                  maxnum_funcvalues_short_pair,&
                  num_funcvalues_short_pair,symfunction_short_pair_list)
              endif
            endif ! lshort
!!
!!=============================================================================
!! write functione.data
!!=============================================================================
            if(lelec.and.(nn_type_elec.eq.1))then
              call writeatomicsymfunctions(symeunit,i1,&
                maxnum_funcvalues_elec,num_funcvalues_elec,symfunction_elec_list)
            endif ! lelec
!!
!!=============================================================================
!! write function_<>.data for nntb
!!=============================================================================
           if(lnntb) then
              if(nntb_flag(3)) then
                call writehextoffsymfunctions(symhextoffunit,i1,&
                maxnum_funcvalues_hextoff,num_funcvalues_hextoff,symfunction_hextoff_list)
              endif
            endif ! lnntb
!!=============================================================================
!! write trainstruct.data
!!=============================================================================
!! CHANGE ANDI: GFORTRAN: gfortran has different default width for logicals (1 instead of 2 for ifort),
!!                        this results in a missing space between i8 and l -> problem when reading file,
!!                        so I inserted an additional space. 
           !write(trainstructunit,'(i8,l)')numtrain+1,lperiodic_list(i1)
            write(trainstructunit,'(i8,tr1,l)')numtrain+1,lperiodic_list(i1)
!! END CHANGE
            if(lperiodic_list(i1))then
              do i2=1,3
                write(trainstructunit,'(3f20.14)')(lattice_list(i2,i3,i1),i3=1,3)
              enddo
            endif
            do i2=1,num_atoms_list(i1)
!! Here we write the total forces 
              write(trainstructunit,'(i3,8(x,f15.10))')zelem_list(i1,i2),&
                (xyzstruct_list(i3,i2,i1),i3=1,3),atomcharge_list(i1,i2),&
                atomenergy_list(i1,i2),(totalforce_list(i3,i2,i1),i3=1,3)
            enddo ! i2
!!
            if(luseforces)then
!!=============================================================================
!! write trainforces.data
!!=============================================================================
              if(lshort)then
                write(trainfunit,'(i8)')numtrain+1
                do i2=1,num_atoms_list(i1)
                  write(trainfunit,'(3f15.10)')(shortforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lshort
!!
!!=============================================================================
!! write trainforcese.data
!! FOR FITTING WE DO NOT NEED THIS FILE because electrostatic forces can be calculated exactly and are not fitted
!!=============================================================================
              if(lelec)then
                write(trainfeunit,'(i8)')numtrain+1
                do i2=1,num_atoms_list(i1)
                  write(trainfeunit,'(3f15.10)')(elecforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lelec
!!
            endif ! luseforces
!!
            numtrain=numtrain+1
            write(ounit,*)pointnumber,' Point is used for training ',numtrain
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else ! point is in test set
!!
!!=============================================================================
!! write testing.data
!!=============================================================================
            if(lshort)then
              if(nn_type_short.eq.1)then
                call writeatomicsymfunctions(tymunit,i1,&
                  maxnum_funcvalues_short_atomic,&
                  num_funcvalues_short_atomic,symfunction_short_atomic_list)
              elseif(nn_type_short.eq.2)then
                call writepairsymfunctions(tymunit,i1,pairs_charge_list,&
                  maxnum_funcvalues_short_pair,&
                  num_funcvalues_short_pair,symfunction_short_pair_list)
              endif 
            endif ! lshort
!!
!! write testinge.data
            if(lelec.and.(nn_type_elec.eq.1))then
              call writeatomicsymfunctions(tymeunit,i1,&
                maxnum_funcvalues_elec,num_funcvalues_elec,&
                symfunction_elec_list)
            endif ! lelec
!!=============================================================================
!! write function_<>.data for nntb
!!=============================================================================
            if(lnntb) then
              if(nntb_flag(3)) then
                call writehextoffsymfunctions(tymhextoffunit,i1,&
                maxnum_funcvalues_hextoff,num_funcvalues_hextoff,&
                symfunction_hextoff_list)
              endif
            endif ! lnntb
!!=============================================================================
!!
!! write teststruct.data
!! CHANGE ANDI: GFORTRAN: gfortran has different default width for logicals (1 instead of 2 for ifort),
!!                        this results in a missing space between i8 and l -> problem when reading file,
!!                        so I inserted an additional space. 
           !write(teststructunit,'(i8,l)')numtest+1,lperiodic_list(i1)
            write(teststructunit,'(i8,tr1,l)')numtest+1,lperiodic_list(i1)
!! END CHANGE

            if(lperiodic_list(i1))then
              do i2=1,3
                write(teststructunit,'(3f20.14)')(lattice_list(i2,i3,i1),i3=1,3)
              enddo
            endif ! lperiodic
            do i2=1,num_atoms_list(i1)
              write(teststructunit,'(i3,8(x,f15.10))')zelem_list(i1,i2),&
                (xyzstruct_list(i3,i2,i1),i3=1,3),atomcharge_list(i1,i2),&
                atomenergy_list(i1,i2),(totalforce_list(i3,i2,i1),i3=1,3)
            enddo ! i2
!!
!!
            if(luseforces)then
!!
!! write testforces.data
              if(lshort)then
                write(testfunit,'(i8)')numtest+1
                do i2=1,num_atoms_list(i1)
                  write(testfunit,'(3f15.10)')(shortforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lshort
!!
!! write testforcese.data
!! FOR FITTING WE DO NOT NEED THIS FILE
              if(lelec)then
                write(testfeunit,'(i8)')numtest+1
                do i2=1,num_atoms_list(i1)
                  write(testfeunit,'(3f15.10)')(elecforce_list(i3,i2,i1),i3=1,3)
                enddo ! i2
              endif ! lelec
            endif ! luseforces
!!
            numtest=numtest+1
            write(ounit,*)pointnumber,' Point is used for testing ',numtest
          endif ! random
!!
        else ! point is rejected because of lfitethres, lfitfthres or lrmin
!!
          numrej=numrej+1
          if(.not.lrmin(i1))then
            write(ounit,*)pointnumber,' Point is rejected (too short bond) ',numrej
          elseif(lfitethres.and.(totalenergy_list(i1).gt.(num_atoms_list(i1)*fitethres)))then
            write(ounit,*)pointnumber,' Point is rejected (high E) ',numrej
          elseif(lfitfthres.and.(.not.lforceok))then
            write(ounit,*)pointnumber,' Point is rejected (high F) ',&
              numrej,highfcomponent
          else
            write(ounit,*)'ERROR: point rejected for unknown reason ',numrej
            stop
          endif
!!
        endif ! lfitethres
!!
      enddo ! i1 loop over all structures
!!
      if(nn_type_short.eq.2)then
        deallocate(pairs_charge_list) 
        deallocate(pairs_charge)             
      endif
!!
      
      return
!!
      end
