!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - prediction.f90
!!
      subroutine predictelec(&
        num_atoms,zelem,lattice,xyzstruct,&
        minvalue_elec,maxvalue_elec,avvalue_elec,chargemin,chargemax,&
        nnelecenergy,nnatomcharge,nntotalcharge,nnelecforce,&
        nnstress_elec,sense,lperiodic)
!!
      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions
      use symfunctions
      use nnewald
      use timings
      use predictionoptions
!!
      implicit none
!!
      integer i1                                                   ! internal
      integer npoints                                              ! internal
      integer ncount                                               ! internal
      integer ndone                                                ! internal
      integer zelem(max_num_atoms)                                 ! in
      integer num_atoms                                            ! in
      integer n_start,n_end                                        ! internal
      integer natoms                                               ! internal
      integer, dimension(:), allocatable :: atomindex              ! internal
      integer, dimension(:,:), allocatable :: neighboridx_elec     ! internal
      integer, dimension(:,:), allocatable :: invneighboridx_elec  ! internal
      integer, allocatable :: lsta(:,:)                            ! numbers of neighbors
      integer, allocatable :: lstc(:)                              ! identification of atom
      integer, allocatable :: lste(:)                              ! nuclear charge of atom
      integer, allocatable :: num_neighbors_elec(:)                ! internal
      integer max_num_neighbors_elec                               ! internal
!!
      real*8 lattice(3,3)                                          ! in
      real*8 xyzstruct(3,max_num_atoms)                            ! in
      real*8, dimension(:,:)  , allocatable :: symfunctione        ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: dsfuncdxyze     ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: strse           ! internal
      real*8, dimension(:,:)  , allocatable :: dchargedsfunc       ! internal
      real*8, allocatable :: lstb(:,:)                             ! xyz and r_ij
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)           ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)           ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)            ! in
      real*8 nntotalcharge                                         ! out
      real*8 nnelecforce(3,max_num_atoms)                          ! out
      real*8 nnelecforcepart(3,max_num_atoms)                      ! internal 
      real*8 nnatomcharge(max_num_atoms)                           ! out
      real*8 nnstress_elec(3,3)                                    ! out 
      real*8 nnelecenergy ! total electrostatic energy             ! out
      real*8 erecip                                                ! internal 
      real*8 nnelecenergypart                                      ! internal
      real*8 sense(nelem,maxnum_funcvalues_elec)                   ! out
      real*8 chargemin(nelem)                                      ! in
      real*8 chargemax(nelem)                                      ! in
!!
      logical lperiodic                                            ! in
      logical lrmin                                                ! internal
      logical lextrapolation                                       ! internal
!!
!!======================================================================
!! memory management strategy: keep only 'nblock' atoms in memory at once
!! => loop step by step over these blocks of atoms 
!! In parallel runs the block of atoms is further split among the processes and each process does natoms atoms.
!!======================================================================
!! initialize auxiliary counters for splitting of atoms
        ncount  = num_atoms ! total number of atoms in structure
        npoints = 0 ! number of atoms to be calculated in this loop step by all processes together
        ndone   = 0 ! number of atoms calculated in previous loops
!!
!!======================================================================
!! Start atomic charge calculation for electrostatics for this block of atoms
!!======================================================================
!!
!! next block of atoms 
 12     continue
        if(ncount.gt.nblock)then
          npoints=nblock
          ncount=ncount-nblock
        else
          npoints=ncount
          ncount=ncount-npoints
        endif
!!
!!======================================================================
!! preparations for parallel runs 
!!======================================================================
        if((mpirank.eq.0).and.(.not.lmd))then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a,i8,a)')'This cycle calculates total ',npoints,' atoms'
        endif !'
!! each process now picks 'natoms' atoms from the total number of 'npoints' atoms of this block
        call mpifitdistribution(npoints,natoms,n_start,n_end)
!! adjust position by atoms already done:
        n_start=n_start+ndone
        n_end  =n_end  +ndone
        if(.not.lmd)then
          write(ounit,'(a,i6,a,i8,a,i8,a,i8)')&
            ' process ',mpirank,' calculates atoms ',n_start,' to ',n_end,' natoms: ',natoms
        endif
        call mpi_barrier(mpi_comm_world,mpierror)
!! determine the atomindex array '
        allocate(atomindex(natoms))
        do i1=1,natoms
          atomindex(i1)=n_start+i1-1
        enddo
!!======================================================================
!! end preparations for parallel runs 
!!======================================================================
!!
!!======================================================================
!! calculate NN charges only if they are not fixed
!!======================================================================
        if(nn_type_elec.eq.1)then
          allocate(lsta(2,max_num_atoms))
          allocate(lstc(listdim))
          allocate(lste(listdim))
          allocate(lstb(listdim,4))
          allocate(num_neighbors_elec(num_atoms))
!!======================================================================
!! get num_neighbors_elec and max_num_neighbors_elec for electrostatic part
!!======================================================================
          call getneighborsatomic(&
            num_atoms,num_neighbors_elec,zelem,&
            max_num_neighbors_elec,&
            lsta,lstc,lste,&
            maxcutoff_elec,lattice,xyzstruct,lstb,lperiodic)
!!======================================================================
!! get neighboridx_elec and invneighboridx_elec for electrostatic part 
!!======================================================================
          allocate(neighboridx_elec(natoms,0:max_num_neighbors_elec))
          allocate(invneighboridx_elec(natoms,max_num_atoms))
          call getneighboridxatomic_para(n_start,natoms,listdim,&
            max_num_atoms,max_num_neighbors_elec,&
            lsta,lstc,neighboridx_elec,invneighboridx_elec)
!!
!!======================================================================
!! calculation of short range atomic energies (nnatomenergy), 
!! forces (nnshortforce), and stress (nnstress_short)
!!======================================================================
          call getchargesatomic(n_start,natoms,atomindex,max_num_neighbors_elec,&
            invneighboridx_elec,zelem,&
            lsta,lstc,lste,lstb,xyzstruct,&
            sense,minvalue_elec,maxvalue_elec,avvalue_elec,&
            scmin_elec,scmax_elec,nnatomcharge,lextrapolation)
!!
!!======================================================================
!! for fixed charges set charges here
!!======================================================================
        elseif(nn_type_elec.eq.3)then  
          do i1=n_start,n_end
            nnatomcharge(i1)=fixedcharge(elementindex(zelem(atomindex(i1))))
          enddo ! i1    
        elseif(nn_type_elec.eq.4)then  
          call readcharges(max_num_atoms,num_atoms,nnatomcharge)
        endif 
!!
!!======================================================================
!! deallocate electrostatic arrays depending on natoms of this block of atoms 
!!======================================================================
        if(nn_type_elec.eq.1)then
          deallocate(lsta)
          deallocate(lstc)
          deallocate(lste)
          deallocate(lstb)
          deallocate(neighboridx_elec)
          deallocate(invneighboridx_elec)
          deallocate(num_neighbors_elec)
        endif
        deallocate(atomindex)
!!
!! update the number of finished atoms:
!!======================================================================
!! End atomic charge calculation for electrostatics for this block of atoms 
!! Now we have only the charges, but no electrostatic forces and stress
!!======================================================================
        ndone=ndone+npoints
!! if there are atoms left to be done go to back and do next block of atoms
        if(ncount.gt.0) goto 12
!!
!!======================================================================
!!======================================================================
!! Now block-wise looping is completed
!!======================================================================
!!======================================================================
!!
!!======================================================================
!! finalize electrostatic charge part  
!! combine the atomic charges of all processes
!!======================================================================
        if(lfinetime)then
          daycomm1=0
          call abstime(timecomm1start,daycomm1)
        endif ! lfinetime
!!
!! if we do not use fixed charges distribute the atomic charges
        if(nn_type_elec.eq.1)then
          call mpi_allreduce(mpi_in_place,nnatomcharge,max_num_atoms,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
          if(lsens)then
            call mpi_allreduce(mpi_in_place,sense,nelem*maxnum_funcvalues_elec,&
              mpi_real8,mpi_sum,mpi_comm_world,mpierror)
          endif ! lsens
        endif !  
!!
        if(lfinetime)then
          call abstime(timecomm1end,daycomm1)
          timecomm1=timecomm1+timecomm1end-timecomm1start
        endif ! lfinetime
!!======================================================================
!! calculate the NN total charge and check for extrapolation
!! TODO: in principle this subroutine should not write anything to runner.out
!!======================================================================
        nntotalcharge=0.0d0
        do i1=1,num_atoms
          nntotalcharge=nntotalcharge+nnatomcharge(i1)
          if((mpirank.eq.0).and.(.not.lmd).and.&
            (nn_type_elec.ne.3).and.(nn_type_elec.ne.4))then
            if(nnatomcharge(i1).gt.chargemax(elementindex(zelem(i1))))then
              write(ounit,*)'-------------------------------------------------------------'
              write(ounit,'(a,i8,a,x,a2)')' WARNING: charge .gt. chargemax for atom ',&
                i1,' element ',element(elementindex(zelem(i1)))
              write(ounit,'(a,2f20.10)')' chargemax, nnatomcharge ',&
                chargemax(elementindex(zelem(i1))),nnatomcharge(i1)
              write(ounit,*)'-------------------------------------------------------------'
            endif
            if(nnatomcharge(i1).lt.chargemin(elementindex(zelem(i1))))then
              write(ounit,*)'-------------------------------------------------------------'
              write(ounit,'(a,i8,a,x,a2)')' WARNING: charge .lt. chargemin for atom ',&
                i1,' element ',element(elementindex(zelem(i1)))
              write(ounit,'(a,2f20.10)')' chargemin, nnatomcharge ',&
                chargemin(elementindex(zelem(i1))),nnatomcharge(i1)
              write(ounit,*)'-------------------------------------------------------------'
            endif
          endif ! mpirank
        enddo
!!======================================================================
!! enforce total charge if requested 
!!======================================================================
        if(enforcetotcharge.eq.0)then
!! do nothing
        elseif(enforcetotcharge.eq.1)then
          call enforcecharge(nntotalcharge,nnatomcharge)
        else
          write(ounit,*)'Error: unknown enforce_totcharge flag ',enforcetotcharge
          stop
        endif ! enforcetotalcharge
!!
!!======================================================================
!!======================================================================
!! Now we are done with the first loop over all atoms.
!! We have now
!! - short range energy, forces and stress for all atoms
!! - charges on all atoms
!!======================================================================
!!======================================================================
!!
!!======================================================================
!!======================================================================
!! Now complete the electrostatic calculations: get electrostatic energy, forces and stress
!! This can only be done here (not in the loop above), because first 
!! each process needs to know all atomic charges before the derivatives
!! of the charges with respect to the atomic positions can be calculated
!!======================================================================
!!======================================================================
!!
!!======================================================================
!! Now we have to calculate the electrostatic energies and forces from the charges determined above
!! In the old version we needed to calculate dchargedxyz for all atoms (large array).
!! Now dchargedxyz is not calculated and stored explicitly anymore, instead its 
!! components dchargedsfunc and dsfuncdxyze are used directly to avoid the array dchargedxyz
!!======================================================================
!!
!!======================================================================
!! Now again loop block-wise over atoms and parallelize each block if requested
!!======================================================================
        ncount  = num_atoms ! total number of atoms in structure
        npoints = 0 ! number of atoms to be calculated in this loop step by all processes together
        ndone   = 0 ! number of atoms calculated in previous loops

!! start next block of atoms
 13     continue
        if(ncount.gt.nblock)then
          npoints=nblock
          ncount=ncount-nblock
        else
          npoints=ncount
          ncount=ncount-npoints
        endif
!!
!!======================================================================
!! preparations for parallel runs 
!!======================================================================
        call mpifitdistribution(npoints,natoms,n_start,n_end)
!! adjust position by atoms already done:
        n_start=n_start+ndone
        n_end  =n_end  +ndone
!! determine the atomindex array
        allocate(atomindex(natoms))
        do i1=1,natoms
          atomindex(i1)=n_start+i1-1
        enddo
        if(.not.lmd)then
          write(ounit,'(a,i6,a,i8,a,i8,a,i8)')&
            ' process ',mpirank,' calculates electrostatic E and F for atoms '&
            ,n_start,' to ',n_end,' natoms: ',natoms !'
        endif
        call mpi_barrier(mpi_comm_world,mpierror)
!!======================================================================
!! end preparations for parallel runs 
!!======================================================================
!!
        allocate(lsta(2,max_num_atoms))
        allocate(lstc(listdim))
        allocate(lste(listdim))
        allocate(lstb(listdim,4))
        allocate(num_neighbors_elec(num_atoms))
!!======================================================================
!! get num_neighbors_elec, max_num_neighbors_elec and neighborlists
!!======================================================================
        call getneighborsatomic(&
          num_atoms,num_neighbors_elec,zelem,&
          max_num_neighbors_elec,&
          lsta,lstc,lste,&
          maxcutoff_elec,lattice,xyzstruct,lstb,lperiodic)
!!======================================================================
!! get neighboridx_elec and invneighboridx_elec 
!!======================================================================
        allocate(neighboridx_elec(natoms,0:max_num_neighbors_elec))
        allocate(invneighboridx_elec(natoms,max_num_atoms))
        call getneighboridxatomic_para(n_start,natoms,listdim,&
          max_num_atoms,max_num_neighbors_elec,&
          lsta,lstc,neighboridx_elec,invneighboridx_elec)
!!
!!======================================================================
!! allocate arrays for this block of natoms atoms 
!!======================================================================
        allocate(symfunctione(maxnum_funcvalues_elec,natoms))
        allocate(dsfuncdxyze(maxnum_funcvalues_elec,natoms,0:max_num_neighbors_elec,3)) 
        allocate(strse(3,3,maxnum_funcvalues_elec,natoms))
        allocate(dchargedsfunc(natoms,maxnum_funcvalues_elec))
!!
!!======================================================================
!! calculate symfunctione, dsfuncdxyze and strse for natom atoms
!!======================================================================
        symfunctione(:,:)   =0.0d0
        dsfuncdxyze(:,:,:,:)=0.0d0
        strse(:,:,:,:)      =0.0d0
!!======================================================================
!! calculate symfunctione, dsfuncdxyze and strse for one atom, we need only dsfuncdxyze 
!!======================================================================
!! TODO: do not calculate dsfuncdxyze and dchargedsfunc if we don't want forces
        if(nn_type_elec.eq.1)then ! if we do not use fixed charges
          if(lfinetime)then
            daysymelec2=0
            call abstime(timesymelec2start,daysymelec2)
          endif ! lfinetime
          lrmin=.true.
          call calconefunction_atomic(cutoff_type,max_num_neighbors_elec,&
            max_num_atoms,n_start,natoms,atomindex,natoms,elementindex,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            nelem,zelem,listdim,&
            lsta,lstc,lste,invneighboridx_elec,&
            function_type_elec,symelement_elec,&
            xyzstruct,symfunctione,rmin,&
            funccutoff_elec,eta_elec,rshift_elec,lambda_elec,&
            zeta_elec,dsfuncdxyze,strse,lstb,&
            lperiodic,ldoforces,ldostress,lrmin)
          if(.not.lrmin)then
            write(ounit,*)'Error in predictelec: lrmin=.false. (atoms too close),rmin= ',rmin
            stop !'
          endif
!!======================================================================
!! scale the symmetry functions for the charge prediction
!! caution: internally nblock and npoints are set to 1 to avoid _list in zelem, symfunctione and num_atoms
!!======================================================================
          call scalesym_para(natoms,atomindex,&
            nelem,1,1,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            zelem,symfunctione,&
            minvalue_elec,maxvalue_elec,avvalue_elec,&
            scmin_elec,scmax_elec)
!!======================================================================
!! we also need to scale the derivative terms dsfuncdxyze and strse 
!!======================================================================
          if(lscalesym)then
            call scaledsfunc_para(natoms,atomindex,max_num_neighbors_elec,&
              maxnum_funcvalues_elec,num_funcvalues_elec,&
              nelem,minvalue_elec,maxvalue_elec,&
              scmin_elec,scmax_elec,&
              zelem,dsfuncdxyze,strse)
          endif
        endif ! 
        if(lfinetime)then
          call abstime(timesymelec2end,daysymelec2)
          timesymelec2=timesymelec2+timesymelec2end-timesymelec2start
        endif ! lfinetime
!!
!!======================================================================
!! Calculate dchargedsfunc for natoms atoms
!!======================================================================
        if(nn_type_elec.eq.1)then ! if we do not use fixed charges
          call getdchargedsfunc_para(natoms,atomindex,&
            zelem,symfunctione,dchargedsfunc)
        endif ! nn_type_elec.eq.1 
!!
!!======================================================================
!! Now calculate the electrostatic forces from charges, dsfuncdxyze and dchargedsfunc
!!======================================================================
        nnelecenergypart=0.0d0
        nnelecforcepart(:,:)=0.0d0
        if(nn_type_elec.eq.3)then ! if we use fixed charges
          dchargedsfunc(:,:)=0.0d0
        elseif(nn_type_elec.eq.4)then ! if we read fixed charges from file
          dchargedsfunc(:,:)=0.0d0
!!          call readcharges(max_num_atoms,num_atoms,nnatomcharge) !! commented because charges are read above already
        endif
        if(lperiodic)then
          if((nn_type_elec.eq.1).or.(nn_type_elec.eq.3).or.(nn_type_elec.eq.4))then
            call getewald(n_start,n_end,natoms,atomindex,&
              num_funcvalues_elec,zelem,num_atoms,&
              max_num_neighbors_elec,&
              invneighboridx_elec,&
              nnatomcharge,lattice,xyzstruct,dchargedsfunc,dsfuncdxyze,&
              maxcutoff_elec,nnelecenergypart,erecip,nnelecforcepart,ldoforces,lperiodic)
          else
            write(ounit,*)'FIX predictelec: nn_type_elec not implemented'
            stop !'
          endif
        else ! not periodic
          if((nn_type_elec.eq.1).or.(nn_type_elec.eq.3).or.(nn_type_elec.eq.4))then
            call coulomb_para(natoms,num_atoms,atomindex,zelem,&
              num_funcvalues_elec,max_num_neighbors_elec,invneighboridx_elec,&
              nnatomcharge,xyzstruct,dsfuncdxyze,dchargedsfunc,&
              maxcutoff_elec,nnelecenergypart,nnelecforcepart,ldoforces)
          else
            write(ounit,*)'FIX predictelec: nn_type_elec not implemented'
            stop
          endif
        endif ! lperiodic
!!
        ndone=ndone+npoints
!!
!!======================================================================
!! Now combine electrostatic energy and force contributions from all processes
!!======================================================================
        call mpi_allreduce(mpi_in_place,nnelecenergypart,&
          1,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        nnelecenergy=nnelecenergy+nnelecenergypart
        call mpi_allreduce(mpi_in_place,nnelecforcepart,&
          max_num_atoms*3,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        nnelecforce(:,:)=nnelecforce(:,:)+nnelecforcepart(:,:)
!!
!!======================================================================
!! deallocate arrays for this block of atoms
!!======================================================================
        deallocate(atomindex)
        deallocate(symfunctione)
        deallocate(dsfuncdxyze)
        deallocate(strse)
        deallocate(dchargedsfunc)
        deallocate(lsta)
        deallocate(lstc)
        deallocate(lste)
        deallocate(lstb)
        deallocate(neighboridx_elec)
        deallocate(invneighboridx_elec)
        deallocate(num_neighbors_elec)
!!
!!======================================================================
!! if there are atoms left go to next block of atoms
!!======================================================================
        if(ncount.gt.0) goto 13
!!
!! Now add reciprocal space energy contribution
        if(lperiodic)then
          nnelecenergy=nnelecenergy+erecip
        endif
!!TODO: check:
!!        nnewald=nnelecenergy/dble(num_atoms)
!!
!! calculation of electrostatic stress 
        if(ldostress.and.lperiodic)then
          write(ounit,*)'### WARNING ### electrostatic stress is not implemented'
          stop !'
        endif ! ldostress
!!
      return
      end
