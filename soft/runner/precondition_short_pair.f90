!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine precondition_short_pair(&
         ntrain,trainelempair,&
         minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
         eshortmin,eshortmax,eshortav,eshortstddev)
!!
      use mpi_mod
      use fileunits
      use fittingoptions 
      use nnflags
      use globaloptions 
      use structures
      use symfunctions
      use nnshort_pair
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5                              ! internal
      integer ntrain                                      ! in
      integer ndummy                                      ! internal
      integer trainelempair(npairs)                       ! in
      integer imaxerror_dummy                             ! internal
      integer nenergies                                   ! internal
      integer nforces                                     ! internal
      integer nforcest                                    ! internal
      integer ncount                                      ! internal
      integer netot                                       ! internal
      integer nstruct                                     ! internal
      integer n_start                                     ! internal
      integer n_end                                       ! internal
      integer icount                                      ! internal
      integer npoints                                     ! internal
      integer, dimension(:), allocatable :: num_pairs_mpi ! internal
      integer, dimension(:), allocatable :: num_atoms_mpi ! internal
      integer, dimension(:,:), allocatable :: zelem_mpi   ! internal
      integer, dimension(:,:,:), allocatable :: zelemp_mpi ! internal
!!
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)         ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)         ! in
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)          ! in
      real*8 eshortav                                     ! in
      real*8 eshortmin                                    ! in
      real*8 eshortmax                                    ! in
      real*8 eshortstddev                                 ! in
      real*8 minnode(maxnodes_short_pair,maxnum_layers_short_pair,npairs)            ! internal
      real*8 maxnode(maxnodes_short_pair,maxnum_layers_short_pair,npairs)            ! internal
      real*8 avnode(maxnodes_short_pair,maxnum_layers_short_pair,npairs)            ! internal
      real*8 symfunction_pair(maxnum_funcvalues_short_pair)                       ! internal
      real*8 weightsp(maxnum_weights_short_pair)                               ! internal
      real*8 nnoutput
      real*8 nodes_valuesp(maxnum_layers_short_pair,maxnodes_short_pair) ! 
      real*8 nodes_sump(maxnum_layers_short_pair,maxnodes_short_pair)    !
      real*8 nneshortmax                              ! internal
      real*8 nneshortmin                              ! internal
      real*8 nneshort_list(nblock)                    ! internal
      real*8 nnetot_list(nblock)                      ! internal
      real*8 edummy                                   ! internal
      real*8 nneshortsum                              ! internal
      real*8 nnstddev                                 ! internal
      real*8 nnstddevq(nelem)                         ! internal
      real*8, dimension(:), allocatable :: shortenergy_mpi       ! internal
      real*8, dimension(:,:,:), allocatable :: xyzstruct_mpi     ! internal
      real*8, dimension(:,:,:), allocatable :: symfunctione_mpi  ! internal
      real*8, dimension(:), allocatable ::  nneshort_mpi         ! internal
      real*8, dimension(:,:,:), allocatable :: nnshortforce_mpi  ! dummy here
      real*8 prefactor                                          ! internal
      real*8 bias                                               ! internal
      real*8 rmse_dummy                                         ! dummy
      real*8 mad_dummy                                          ! dummy
      real*8 maxerror_dummy                                     ! dummy
      real*8 avnneshort                                         ! out
!!
      logical, dimension(:), allocatable :: lperiodic_mpi       ! internal
!!
      ndummy = 0
!!
      if(mpisize.gt.1)then
        write(*,*)'ERROR: preconditioning does not yet work in parallel runs'
        stop !'
      endif
!!
!! write general header
      if(mpirank.eq.0)then
        write(ounit,*)'============================================================='
        write(ounit,*)'Weight Preconditioner:'     !'
        write(ounit,*)'----------------------'
      endif ! mpirank
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Preconditioner for the nodes in the hidden layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! short range part
      if(lshort)then
        write(ounit,*)'Short range part:'
        write(ounit,*)'-----------------'
!!!!!!!!!!!!!!!!
!! outermost loop: over all elemental pairs, because they can have different numbers of hidden layers
        avnode(:,:,:) =0.0d0
        minnode(:,:,:)=10000000.0d0
        maxnode(:,:,:)=-10000000.0d0
        do i3=1,npairs
          do i4=1,num_layers_short_pair(i3)-1       
            write(ounit,'(a,i4,a4,a2,x,a2)')' Preconditioning hidden layer ',&
              i4,' of ',element(elementindex(elempair(i3,1))),element(elementindex(elempair(i3,2)))

!! open files
            if(mpirank.eq.0)then
              open(symunit,file='function.data',form='formatted',status='old')
              rewind(symunit) !'
              open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
              rewind(trainstructunit) !'
            endif ! mpirank.eq.0
!!
            ncount = ntrain
 30         continue
            if(ncount.gt.nblock)then
              npoints=nblock
              ncount=ncount-nblock
            else
              npoints=ncount
              ncount=ncount-npoints
            endif
!!
!! read symmetry functions
            if(mpirank.eq.0)then
              call readfunctions(2,symunit,npoints,npairs,&
                max_num_pairs,maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
                symfunction_short_pair_list)
            endif ! mpirank.eq.0
!!
!! loop over all npoints structures and all atoms in these structures
            do i1=1,npoints
              do i2=1,num_pairs_list(i1)
!!
!! loop over all hidden layers, note: different elements can have different numbers of hidden layers!
!! identify the elemental pair
                if(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)).eq.i3)then
!!
                  symfunction_pair(:)=symfunction_short_pair_list(:,i2,i1)
                  weightsp(:)= weights_short_pair(:,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)))  
!!
!! calculate nodes_valuesp and nnoutput for pair i1
                  call calconenn(1,maxnum_funcvalues_short_pair,maxnodes_short_pair,maxnum_layers_short_pair,&
                    num_layers_short_pair(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))),&
                    maxnum_weights_short_pair,nodes_short_pair(0,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))),&
                    symfunction_pair,weights_short_pair,nodes_valuesp,nodes_sump,nnoutput,&
                    actfunc_short_pair(1,1,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))))
!!
!! loop over all nodes of this hidden layer
                  do i5=1,nodes_short_pair(i4,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)))
                    bias=weights_short_pair(windex_short_pair(i4,i3)+i5-1,i3)
                    avnode(i5,i4,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)))&
                         =avnode(i5,i4,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))) +nodes_sump(i4,i5)-bias
                    minnode(i5,i4,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)))&
                         =min(minnode(i5,i4,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))),nodes_sump(i4,i5)-bias)
                    maxnode(i5,i4,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)))&
                         =max(maxnode(i5,i4,pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))),nodes_sump(i4,i5)-bias)
                  enddo
                endif ! pair separation
              enddo ! i2 all independent pairs
            enddo ! i1 structure
            if(ncount.gt.0) goto 30
!!
!! analyze results for this hidden layer here
            do i5=1,nodes_short_pair(i4,i3)
              if(trainelempair(i3).gt.0)then !! just if pairs of this type are present!
                avnode(i5,i4,i3)=avnode(i5,i4,i3)/dble(trainelempair(i3))
              endif
            enddo ! i5
!! modify the connecting weights ending at layer i4 here:

!! TODO

!! shift the bias weights to the negative center of mass of the linear combination of input nodes
            icount=windex_short_pair(2*i4,i3)
            do i5=1,nodes_short_pair(i4,i3)
              weights_short_pair(icount,i3)=-avnode(i5,i4,i3) ! pre3
              icount=icount+1
            enddo ! i5
          enddo ! i4 hidden layer
        enddo ! i3 npairs
!!
      endif ! lshort
!!
!! end block of training points
!!
!! close files
      if(mpirank.eq.0)then
        close(symunit)
        close(trainstructunit)
      endif ! mpirank.eq.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Preconditioner for the connecting weights to the output node and 
!! the bias of the output node 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the initial training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! initializations
      nneshortmin              = 100000000.d0
      nneshortmax              = -100000000.d0
      nenergies                = 0
      nforces                  = 0
      nforcest                 = 0
      ncount                   = ntrain
      nneshort_list(:)         = 0.0d0
      nnetot_list(:)           = 0.0d0
      edummy                   = 1.d12
      netot                    = 0
      nneshortsum              = 0.0d0
      nnstddev                 = 0.0d0
      nnstddevq(:)             = 0.0d0
!!
!! open files
      if(mpirank.eq.0)then
        open(symunit,file='function.data',form='formatted',status='old')
        rewind(symunit)
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
      endif ! mpirank.eq.0
!!
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! determine which structures of this block should be calculated by this process
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! do all file reading for training error here at one place to allow for parallelization
!!--------------------------------------------------------------------------------------
!!
      if(mpirank.eq.0)then
        if(lshort)then
!! read npoint short range data sets     '
          call readfunctions(2,symunit,npoints,npairs,&
            max_num_pairs,maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
            symfunction_short_pair_list)
        endif
!!
!! read the structures needed for the calculation of the forces 
!! and get reference forces from DFT
!! must be called after readfunctions because it needs num_atoms_list
!!
        call getstructures(trainstructunit,npoints)
!!
      endif ! mpirank.eq.0
!!
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelemp_list,nblock*max_num_pairs,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_short_pair_list,nblock*max_num_pairs*maxnum_funcvalues_short_pair,&       !CHK
           mpi_real8,0,mpi_comm_world,mpierror)                                      ! CHK
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
!!
!! end of file reading for training error
!!--------------------------------------------------------------------
!!
!! prepare local input arrays for this process
      allocate(num_pairs_mpi(nstruct))
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
      allocate(zelemp_mpi(2,nstruct,max_num_pairs))
      allocate(shortenergy_mpi(nstruct))
      allocate(xyzstruct_mpi(3,max_num_atoms,nstruct))
      allocate(lperiodic_mpi(nstruct))

!! prepare local output arrays for this process
      allocate(nneshort_mpi(nstruct))
      nneshort_mpi(:)=0.0d0
      allocate(nnshortforce_mpi(3,max_num_atoms,nstruct))
      nnshortforce_mpi(:,:,:)=0.0d0
!!
!!    copy local arrays
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        num_pairs_mpi(icount)       = num_pairs_list(i1)
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        zelemp_mpi(:,icount,:)      = zelemp_list(:,i1,:)
        shortenergy_mpi(icount)     = shortenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
      enddo ! i1
!!
!! get the short range energies for this block of points 
      if(lshort)then
        call getshortenergies_parapair(nstruct,ndummy,&
          nenergies,imaxerror_dummy,&
          num_atoms_mpi,num_pairs_mpi,&
          zelem_mpi,zelemp_mpi,&
          minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
          symfunction_short_pair_list(1,1,n_start),lattice_list(1,1,n_start),&
          nnshortforce_mpi,xyzstruct_mpi,nneshort_mpi,&
          rmse_dummy,mad_dummy,maxerror_dummy,shortenergy_mpi,&
          lperiodic_mpi)
!!
!! copy results back on full array
        icount=0
        nneshort_list(:)       =0.0d0
        do i1=n_start,n_end
          icount=icount+1
          nneshort_list(i1)       = nneshort_mpi(icount)
        enddo
!!
!! distribute results to all processes
        call mpi_allreduce(mpi_in_place,nneshort_list,nblock,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
      endif ! lshort
!!
!!    sum up short range energies
      do i1=1,npoints
        if(shortenergy_list(i1).le.maxenergy)then
          nenergies  =nenergies+1
          nneshortsum=nneshortsum+nneshort_list(i1)
          nneshortmin=min(nneshortmin,nneshort_list(i1))
          nneshortmax=max(nneshortmax,nneshort_list(i1))
        endif
      enddo
!!
      deallocate(nnshortforce_mpi)
      deallocate(nneshort_mpi)
      deallocate(symfunctione_mpi)

      deallocate(lperiodic_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(shortenergy_mpi)
      deallocate(zelemp_mpi)
      deallocate(zelem_mpi)
      deallocate(num_atoms_mpi)
      deallocate(num_pairs_mpi)
!!
      if(ncount.gt.0) goto 10
!! end block of training points
!!
!! close files
      if(mpirank.eq.0)then
        close(symunit)
        close(trainstructunit)
      endif ! mpirank.eq.0
!!
      call mpi_allreduce(mpi_in_place,nenergies,1,&
           mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
!! now get averages
      if(lshort)then 
        avnneshort=nneshortsum/dble(nenergies)
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! Now get the standard deviations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!
      ncount                   = ntrain
!!
!! open files and write headers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mpirank.eq.0)then
        open(symunit,file='function.data',form='formatted',status='old')
        rewind(symunit)
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
      endif ! mpirank.eq.0
!!
 20   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! determine which structures of this block should be calculated by this process
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! do all file reading for training error here at one place to allow for parallelization
!!--------------------------------------------------------------------------------------
!!
      if(mpirank.eq.0)then
        if(lshort)then
!! read npoint short range data sets     '
          call readfunctions(2,symunit,npoints,npairs,&
            max_num_pairs,maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
            symfunction_short_pair_list)
        endif
!!
!! read the structures needed for the calculation of the forced 
!! and get reference forces from DFT
!! must be called after readfunctions because it needs num_atoms_list
!!
        call getstructures(trainstructunit,npoints)
!!
      endif ! mpirank.eq.0
!!
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelemp_list,nblock*max_num_pairs,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_short_pair_list,nblock*max_num_pairs*maxnum_funcvalues_short_pair,&       !CHK
           mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
!!
!! end of file reading for training error
!!--------------------------------------------------------------------
!!
!! prepare local input arrays for this process
      allocate(num_pairs_mpi(nstruct))
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
      allocate(zelemp_mpi(2,nstruct,max_num_pairs))
      allocate(shortenergy_mpi(nstruct))
      allocate(xyzstruct_mpi(3,max_num_atoms,nstruct))
      allocate(lperiodic_mpi(nstruct))

!! prepare local output arrays for this process
      allocate(nneshort_mpi(nstruct))
      nneshort_mpi(:)=0.0d0
      allocate(nnshortforce_mpi(3,max_num_atoms,nstruct))
      nnshortforce_mpi(:,:,:)=0.0d0
!!
!!    copy local arrays
      icount=0
      do i1=n_start,n_end
        icount=icount+1
        num_pairs_mpi(icount)       = num_pairs_list(i1)
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        zelemp_mpi(:,icount,:)      = zelemp_list(:,i1,:)
        shortenergy_mpi(icount)     = shortenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
      enddo ! i1
!!
!! get the short range energies for this block of points
      if(lshort)then
        call getshortenergies_parapair(nstruct,ndummy,&
          nenergies,imaxerror_dummy,&
          num_atoms_mpi,num_pairs_mpi,&
          zelem_mpi,zelemp_mpi,&
          minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
          symfunction_short_pair_list(1,1,n_start),lattice_list(1,1,n_start),&
          nnshortforce_mpi,xyzstruct_mpi,nneshort_mpi,&
          rmse_dummy,mad_dummy,maxerror_dummy,shortenergy_mpi,&
          lperiodic_mpi)
!!
!! copy results back on full array
        icount=0
        nneshort_list(:)    = 0.0d0
        do i1=n_start,n_end
          icount=icount+1
          nneshort_list(i1) = nneshort_mpi(icount)
        enddo
!!
!! distribute results to all processes
        call mpi_allreduce(mpi_in_place,nneshort_list,nblock,&
             mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! calculate part of the standard deviation   
        if(mpirank.eq.0)then
          do i1=1,npoints
            if(shortenergy_list(i1).le.maxenergy)then
              nnstddev=nnstddev+(nneshort_list(i1)-avnneshort)**2.d0
            endif
          enddo
        endif ! mpirank
      endif ! lshort
!!
      deallocate(nnshortforce_mpi)
      deallocate(nneshort_mpi)
      deallocate(symfunctione_mpi)
!!
      deallocate(lperiodic_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(shortenergy_mpi)
      deallocate(zelemp_mpi)
      deallocate(zelem_mpi)
      deallocate(num_atoms_mpi)
      deallocate(num_pairs_mpi)
!!
      if(ncount.gt.0) goto 20
!! end block of training points
!!
!! close files
      if(mpirank.eq.0)then
        close(symunit)
        close(trainstructunit)
      endif ! mpirank.eq.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! calculate final standard deviation of short range energies
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a)')' Final preconditioning of the output values:' 
        write(ounit,*)    '--------------------------------------------'
        if(lshort)then
          nnstddev=nnstddev/dble(nenergies)
          nnstddev=dsqrt(nnstddev)
          write(ounit,'(a,f14.6,a)')' Minimum NN Eshort    ',nneshortmin,' Ha/atom'
          write(ounit,'(a,f14.6,a)')' Minimum Ref Eshort   ',eshortmin,' Ha/atom'
          write(ounit,'(a,f14.6,a)')' Maximum NN Eshort    ',nneshortmax,' Ha/atom'
          write(ounit,'(a,f14.6,a)')' Maximum Ref Eshort   ',eshortmax,' Ha/atom'
          write(ounit,'(a,f14.6,a)')' Average NN Eshort    ',avnneshort,' Ha/atom'
          write(ounit,'(a,f14.6,a)')' Average Ref Eshort   ',eshortav,' Ha/atom'
          write(ounit,'(a,f14.6,a)')' Stddev NN Eshort     ',nnstddev,' Ha/atom'
          write(ounit,'(a,f14.6,a)')' Stddev Ref Eshort    ',eshortstddev,' Ha/atom'
!!'
!! scale the short range connecting weights between final hidden layer and output layer
          prefactor=eshortstddev/nnstddev
          write(ounit,'(a,f14.6)')' Factor for connecting short range weights: ',prefactor
          do i1=1,npairs !'
            do i2=windex_short_pair(2*num_layers_short_pair(i1)-1,i1),windex_short_pair(2*num_layers_short_pair(i1),i1)
              weights_short_pair(i2,i1)=weights_short_pair(i2,i1)*prefactor
            enddo
          enddo

!! shift the short range bias weights
          do i1=1,npairs
            weights_short_pair(num_weights_short_pair(i1),i1)&
              =weights_short_pair(num_weights_short_pair(i1),i1)-(prefactor*avnneshort-eshortav)   
          enddo
        endif ! lshort
!!
      endif ! mpirank
!!
      if(mpirank.eq.0)then
        write(ounit,*)'============================================================='
      endif ! mpirank
!!
      return
      end
