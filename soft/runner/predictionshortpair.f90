!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - predict.f90
!!
      subroutine predictionshortpair(&
                 num_atoms,num_atoms_element,zelem,&
                 lattice,xyzstruct,&
                 minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
                 eshortmin,eshortmax,&
                 nntotalenergy,nnshortforce,&
                 nnatomenergy,nnpairenergy,nnshortenergy,&
                 nnstress_short,pairs_charge,&
                 atomenergysum,sens,lperiodic)
!!
      use mpi_mod
      use fileunits
      use predictionoptions 
      use nnflags
      use globaloptions
      use symfunctions
      use nnshort_pair 
      use timings
!!
      implicit none
!!
      integer zelem(max_num_atoms)                                     ! in
      integer num_atoms                                                ! in
      integer num_atoms_element(nelem)                                 ! in 
      integer num_pairs                                              
      integer npoints                                                  ! internal
      integer ncount                                                   ! internal
      integer ndone                                                    ! internal
      integer num_pairs_para                                           ! out 
      integer pairs_charge(2,max_num_pairs)                            ! in 
      integer i1,i2,i3,i4                                              ! internal
      integer, dimension(:), allocatable :: atomindex                  ! internal
      integer, dimension(:), allocatable :: pindex                     ! internal
      integer n_start                                                  ! internal
      integer n_end                                                    ! internal
      integer, allocatable :: lsta(:,:)                          ! numbers of neighbors
      integer, allocatable :: lstc(:)                            ! identification of atom
      integer, allocatable :: lste(:)                            ! nuclear charge of atom
      integer, allocatable :: num_neighbors(:)                          
      integer, allocatable :: neighboridx_short_pair(:,:)   
      integer, allocatable :: invneighboridx_short_pair(:,:)  
!!
      real*8 lattice(3,3)                                              ! in
      real*8 xyzstruct(3,max_num_atoms)                                ! in
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)                      ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)                      ! in
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)                       ! in
      real*8, dimension(:,:)  , allocatable :: symfunctionp
!! data predicted the NN
      real*8 nntotalenergy                                             ! out 
      real*8 nnshortforce(3,max_num_atoms)                               ! out 
      real*8 nnatomenergy(max_num_atoms)                               ! out 
      real*8 nnpairenergy(max_num_pairs)                               ! out 
!! symmetry function parameters
      real*8, dimension(:,:,:,:)  , allocatable :: dsfuncdxyz_pair 
      real*8, dimension(:,:,:,:)  , allocatable :: strsp 
      real*8 nnshortenergy                                              ! out 
      real*8 nnstress_short(3,3)                                        ! internal
      real*8, dimension(:,:)  , allocatable :: depairdsfunc
      real*8 sens(npairs,maxnum_funcvalues_short_pair)                            ! out 
      real*8 eshortmin                                                  ! in
      real*8 eshortmax                                                  ! in
      real*8 atomenergysum                                              ! out 
      real*8, allocatable :: lstb(:,:)                                  ! xyz and r_ij 
!!
      logical lperiodic                                                 ! in
      logical lrmin                                                     ! internal
      logical lextrapolation                                            ! internal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! General part: I/O, initializations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! initializations for NN data
      nntotalenergy    =0.0d0
      nnatomenergy(:)  =0.0d0
      nnpairenergy(:)  =0.0d0
      nnstress_short(:,:) = 0.0d0
      nnshortforce(:,:)  =0.0d0
!! further initializations
      nnshortenergy     =0.0d0
      atomenergysum     =0.0d0
!! initialization of sensitivities
      sens(:,:)         =0.0d0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! Short range energy part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! memory management: keep only nblock atoms in memory at once
!! => loop over blocks of atoms here
!! In parallel runs the block of atoms is further split among the processes.
!! In this first looping we get the short range energies and charges of all atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! initialize auxiliary counters for splitting of pairs 
      ncount  = max_num_pairs
      npoints = 0 ! number of pairs to be calculated in this loop step by all processes together
      ndone   = 0 ! number of pairs calculated in previous loops
!!
!! next block of points
 15   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! preparations for short-range parallel runs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! debug
      if((mpirank.eq.0).and.(.not.lmd))then
        write(ounit,*)'This cycle calculates total pairs ',npoints
      endif
!!
      call mpifitdistribution(npoints,num_pairs_para,n_start,n_end)
!! adjust position by atoms already done:
      n_start=n_start+ndone
      n_end  =n_end  +ndone
!!
      if(.not.lmd)then
        write(ounit,*)'process ',mpirank,' calculates pairs ',n_start,' to ',n_end
      endif
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
      if((mpirank.eq.0).and.(.not.lmd))then
        write(ounit,*)'-------------------------------------------------------------'
      endif !' mpirank.eq.0
!!
!! determine the pindex array
      num_pairs_para=max_num_pairs       !! FIXME
      allocate(pindex(num_pairs_para))
!! FIXME:
!!      do i2=1,num_pairs_para
!!        pindex(i2)=n_start+i2-1
!!      enddo
      do i2=1,num_pairs_para
        pindex(i2)=i2
      enddo

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end preparations for parallel runs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! Short range energy part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if(lshort)then
!!
        if(mpisize.gt.1)then
          write(ounit,*)'ERROR: predictionpair is not yet parallel'
          stop
        endif
        allocate(symfunctionp(maxnum_funcvalues_short_pair,max_num_pairs)) 
        symfunctionp(:,:)  =0.0d0
        allocate(dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3)) 
        dsfuncdxyz_pair(:,:,:,:)=0.0d0
        allocate(strsp(3,3,maxnum_funcvalues_short_pair,max_num_pairs))  
        strsp(:,:,:,:)=0.0d0
        allocate(depairdsfunc(max_num_pairs,maxnum_funcvalues_short_pair))  
        depairdsfunc(:,:) =0.0d0
!!
!! calculate the short range symmetry functions for pairs n_start to n_end
        call calconefunction_pair(n_start,n_end,num_pairs_para,&
          1,num_atoms,zelem,&
          num_pairs,pairs_charge,&
          lattice,xyzstruct,symfunctionp,&
          dsfuncdxyz_pair,strsp,&
          lperiodic,ldoforces,lrmin)
!!
        if(.not.lrmin)then
          write(ounit,*)'Error in predictionpair: lrmin=.false. (atoms too close in calconefunctionpair_para)'
          stop
        endif
!!
!! print unscaled symmetry functions
        if(ldebug.and.(mpisize.eq.1).and.(.not.lmd))then
          do i1=1,num_pairs_para 
            write(ounit,*)'                          struct  atom  function     symfunction'
            do i2=1,num_funcvalues_short_pair(pairindex(pairs_charge(1,i1),pairs_charge(2,i1)))  
              write(ounit,'(a27,3i6,f20.10)')'SYMFUNCTION SHORT UNSCALED ',1,i1,i2,symfunctionp(i2,i1)
            enddo
          enddo
          write(ounit,*)'-------------------------------------------------------------'
        endif !'
!!
!! check for extrapolation for atoms n_start to n_end
!! CHECK: is it right to do this here before scaling? Probably it is correct.
        call checkextrapolationpair(num_pairs,pindex,npairs,pairindex,&
          max_num_pairs,maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
          pairs_charge,symfunctionp,minvalue_short_pair,maxvalue_short_pair,'short',lextrapolation)

!!
!! scale symmetry functions for the short-range interaction
!! caution: internally nblock and npoints are set to 1 to avoid _list in zelem, symfunction and num_atoms
!! FIXME: pairs_charge is zelemp_list inside scalesympair_para???
      call scalesympair_para(num_pairs,&
         pindex,pairs_charge,symfunctionp,&
         minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)

!! we also need to scale the derivative terms dsfuncdxyz and strs 
        if(lscalesym)then
!! FIXME: pairs_charge is zelemp_list inside caledsfuncpair??
        call scaledsfuncpair(num_pairs_para,&
          minvalue_short_pair,maxvalue_short_pair,&
          pairs_charge,dsfuncdxyz_pair,strsp)
        endif
!!
!! print scaled symmetry functions
        if(ldebug.and.(mpisize.eq.1).and.(.not.lmd))then
          do i1=1,num_pairs_para 
            write(ounit,*)'                        struct  atom  function     symfunction'
            do i2=1,num_funcvalues_short_pair(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1)))) 
              write(ounit,'(a25,3i6,f20.10)')'SYMFUNCTION PAIR SCALED ',1,i1,i2,symfunctionp(i2,i1)
           enddo !'
          enddo
        endif
!!
        if((mpirank.eq.0).and.(.not.lmd))then
          write(ounit,*)'-------------------------------------------------------------'
        endif
!!'
!! predict the short-range atomic energies for atoms n_start to n_end
!!        call calconeshort_parapair(1,num_pairs_para,pindex,&
        call calconeshort_parapair(1,num_pairs,pindex,&
          pairs_charge,symfunctionp,nnpairenergy,&
          nntotalenergy)
!!
        if(ldoforces.or.lsens)then
!! calculation of short range forces 'nnshortforce'
        call getshortforces_parapair(pindex,&
          num_atoms,num_pairs,pairs_charge,&
!!          num_atoms,num_pairs_para,pairs_charge,&
          symfunctionp,&
          dsfuncdxyz_pair,depairdsfunc,nnshortforce)
!!
!! calculation of the sensitivity
          if(lsens)then
            do i1=1,num_pairs_para  
              do i2=1,num_funcvalues_short_pair(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))))
                sens(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))),i2)&
                =sens(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))),i2) + depairdsfunc(i1,i2)
              enddo
            enddo
          endif ! lsens
!!
        endif ! ldoforces
!!
!! calculation of short range stress 
        if(ldostress.and.lperiodic)then
          do i1=1,3
            do i2=1,3
!!              do i4=1,num_pairs_para
              do i4=1,num_pairs
!!                do i3=1,num_funcvaluesp(pairindex(pairs_charge(1,pindex(i4)),pairs_charge(2,pindex(i4))))
                do i3=1,num_funcvalues_short_pair(pairindex(pairs_charge(1,i4),pairs_charge(2,i4)))
                  nnstress_short(i1,i2)=nnstress_short(i1,i2)-strsp(i1,i2,i3,i4)*depairdsfunc(i4,i3)
                enddo ! i3
              enddo ! i4
            enddo ! i2
          enddo ! i1
        endif ! ldostress
!! 
        deallocate(symfunctionp) 
        deallocate(dsfuncdxyz_pair) 
        deallocate(strsp)  
        deallocate(depairdsfunc)  
      endif ! lshort

      ndone=ndone+npoints
      deallocate(pindex)
!!
!! if there are pairs left to be done go to back and do next block of pairs
      if(ncount.gt.0) goto 15
!!
!! combine the nnpairenergy arrays and the total short range energy of all processes
      if(lshort)then
        call mpi_allreduce(mpi_in_place,nnpairenergy,max_num_pairs,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,nntotalenergy,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(ldoforces)then
          call mpi_allreduce(mpi_in_place,nnshortforce,max_num_atoms*3,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif
        if(lsens)then
          call mpi_allreduce(mpi_in_place,sens,npairs*maxnum_funcvalues_short_pair,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif
        if(ldostress.and.lperiodic)then
          call mpi_allreduce(mpi_in_place,nnstress_short,9,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif 
      endif ! lshort
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! General output part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
      nnshortenergy=nntotalenergy
!! check for short range energy extrapolation
      if((mpirank.eq.0).and.lshort.and.(.not.lmd))then
        if((nntotalenergy/dble(max_num_atoms)).gt.eshortmax)then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a)')' WARNING: eshort is .gt. eshortmax '
          write(ounit,'(a,2f20.10)')' average eshort per atom, eshortmax ',&
            nntotalenergy/dble(max_num_atoms),eshortmax
          write(ounit,*)'-------------------------------------------------------------'
        endif
        if((nntotalenergy/dble(max_num_atoms)).lt.eshortmin)then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a)')' WARNING: eshort is .lt. eshortmin '
          write(ounit,'(a,2f20.10)')' average eshort per atom, eshortmin ',&
            nntotalenergy/dble(max_num_atoms),eshortmin
          write(ounit,*)'-------------------------------------------------------------'
        endif
      endif ! mpirank
!!
      return
      end
