!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine precondition_short_atomic(ntrain,trainelem,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        eshortmin,eshortmax,eshortav,eshortstddev)
!!
      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions 
      use fittingoptions
      use symfunctions
      use structures
      use nnshort_atomic
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5  ! internal
      integer imaxerror_dummy ! internal
      integer nforces       ! internal
      integer nforcese      ! internal
      integer nforcest      ! internal
      integer ntrain        ! in
      integer ncount        ! internal
      integer npoints       ! internal
      integer nstruct       ! internal
      integer n_start       ! internal
      integer n_end         ! internal
      integer, dimension(:), allocatable :: num_atoms_mpi  ! internal
      integer, dimension(:,:), allocatable :: zelem_mpi    ! internal
      integer icount          ! internal
      integer nenergies                      ! internal
      integer ndummy                         ! internal
      integer netot                          ! internal
      integer ncharge(nelem)                 ! internal
      integer trainelem(nelem)               ! internal
!!
      real*8 prefactor                       ! internal
      real*8 nneshortmax                     ! internal
      real*8 nneshortmin                     ! internal
      real*8 eshortstddev                    ! in
      real*8 eshortmin                       ! in
      real*8 eshortmax                       ! in
      real*8 eshortav                        ! in
      real*8 nnstddev                        ! internal
      real*8 nnstddevq(nelem)                ! internal
      real*8 rmse_dummy                      ! dummy
      real*8 mad_dummy                       ! dummy
      real*8 maxerror_dummy                  ! dummy
      real*8 symfunction_atom(maxnum_funcvalues_short_atomic)             ! internal
      real*8, dimension(:), allocatable :: shortenergy_mpi       ! internal
      real*8, dimension(:,:,:), allocatable :: xyzstruct_mpi     ! internal
      real*8, dimension(:), allocatable ::  nneshort_mpi         ! internal
      real*8, dimension(:,:,:), allocatable :: nnshortforce_mpi  ! dummy here 
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)           ! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)           ! in
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)            ! in
      real*8 weights(maxnum_weights_short_atomic)        ! internal
      real*8 nneshort_list(nblock)                    ! internal
      real*8 nnetot_list(nblock)                      ! internal
      real*8 edummy                                   ! internal
      real*8 nneshortsum                              ! internal
      real*8 avnneshort                               ! out
      real*8 avnncharge(nelem)                        ! internal
      real*8 nodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic) ! 
      real*8 nodes_sum(maxnum_layers_short_atomic,maxnodes_short_atomic)    !
      real*8 nnoutput
      real*8 minnode(maxnodes_short_atomic,maxnum_layers_short_atomic,nelem)                  ! internal
      real*8 maxnode(maxnodes_short_atomic,maxnum_layers_short_atomic,nelem)                  ! internal
      real*8 avnode(maxnodes_short_atomic,maxnum_layers_short_atomic,nelem)                  ! internal
      real*8 bias                                         ! internal
!!
      logical, dimension(:), allocatable :: lperiodic_mpi    ! internal
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
        write(ounit,*)'Warning: Forces are not used for preconditioning'
        write(ounit,*)'----------------------'
      endif ! mpirank
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Preconditioner for the nodes in the hidden layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if(.true.)then ! temporary switch for development
      if(.false.)then ! temporary switch for development
!!!!!!!!!!!!!!!!
!! short range part
      if(lshort)then
        write(ounit,*)'Short range part:'
        write(ounit,*)'-----------------'
!!!!!!!!!!!!!!!!
!! outermost loop: over all elements, because they can have different numbers of hidden layers
        avnode(:,:,:)=0.0d0
        minnode(:,:,:)=10000000.0d0
        maxnode(:,:,:)=-10000000.0d0
        do i3=1,nelem
          do i4=1,num_layers_short_atomic(i3)-1         ! in
            write(ounit,'(a,i4,a4,a2)')' Preconditioning hidden layer ',i4,' of ',element(i3)

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
              call readfunctions(1,symunit,npoints,nelem,&
                max_num_atoms,maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
                symfunction_short_atomic_list)
            endif ! mpirank.eq.0
!!
!! loop over all npoints structures and all atoms in these structures
            do i1=1,npoints
              do i2=1,num_atoms_list(i1)
!!
!! loop over all hidden layers, note: different elements can have different numbers of hidden layers!
!! identify element
                if(elementindex(zelem_list(i1,i2)).eq.i3)then
                  symfunction_atom(:)=symfunction_short_atomic_list(:,i2,i1)
                  weights(:)= weights_short_atomic(:,elementindex(zelem_list(i1,i2)))  
!!
!! calculate nodes_values and nnoutput for atom i1
                  call calconenn(1,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
                    maxnum_layers_short_atomic,num_layers_short_atomic(elementindex(zelem_list(i1,i2))),&
                    maxnum_weights_short_atomic,nodes_short_atomic(0,elementindex(zelem_list(i1,i2))),&
                    symfunction_atom,weights,nodes_values,nodes_sum,&
                    nnoutput,actfunc_short_atomic(1,1,elementindex(zelem_list(i1,i2))))
!!
!! loop over all nodes of this hidden layer
                  do i5=1,nodes_short_atomic(i4,elementindex(zelem_list(i1,i2)))
                    bias=weights_short_atomic(windex_short_atomic(i4,i3)+i5-1,i3)
                    avnode(i5,i4,elementindex(zelem_list(i1,i2)))&
                      =avnode(i5,i4,elementindex(zelem_list(i1,i2)))&
                      +nodes_sum(i4,i5)-bias
                    minnode(i5,i4,elementindex(zelem_list(i1,i2)))&
                      =min(minnode(i5,i4,elementindex(zelem_list(i1,i2))),nodes_sum(i4,i5)-bias)
                    maxnode(i5,i4,elementindex(zelem_list(i1,i2)))&
                      =max(maxnode(i5,i4,elementindex(zelem_list(i1,i2))),nodes_sum(i4,i5)-bias)
!                    avnode(i5,i4,elementindex(zelem_list(i1,i2)))&
!                      =avnode(i5,i4,elementindex(zelem_list(i1,i2)))&
!                      +nodes_values(i4,i5)
!                    minnode(i5,i4,elementindex(zelem_list(i1,i2)))&
!                      =min(minnode(i5,i4,elementindex(zelem_list(i1,i2))),nodes_values(i4,i5))
!                    maxnode(i5,i4,elementindex(zelem_list(i1,i2)))&
!                      =max(maxnode(i5,i4,elementindex(zelem_list(i1,i2))),nodes_values(i4,i5))
                  enddo
                endif ! element separation
              enddo ! i2 atom
            enddo ! i1 structure
            if(ncount.gt.0) goto 30
!!
!! analyze results for this hidden layer here
            do i5=1,nodes_short_atomic(i4,i3)
              if(trainelem(i3).gt.0)then ! just if atoms of this type are really present
                avnode(i5,i4,i3)=avnode(i5,i4,i3)/dble(trainelem(i3))
                write(ounit,'(a2,x,2i5,3f14.6)')element(i3),i4,i5,&
                  minnode(i5,i4,i3),maxnode(i5,i4,i3),avnode(i5,i4,i3)
              endif
            enddo ! i5
!! modify the connecting weights ending at layer i4 here:

!! TODO

!! shift the bias weights to the negative center of mass of the linear combination of input nodes
            icount=windex_short_atomic(2*i4,i3)
            do i5=1,nodes_short_atomic(i4,i3)
              weights_short_atomic(icount,i3)=-avnode(i5,i4,i3) ! pre3
              icount=icount+1
            enddo ! i5
          enddo ! i4 hidden layer
        enddo ! i3 element
!! write results for this element here

      endif ! lshort
!!
!!
!! end block of training points
!!
!! close files
      if(mpirank.eq.0)then
        close(symunit)
        close(trainstructunit)
      endif ! mpirank.eq.0




      endif ! temporary switch for development

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
      nforcese                 = 0
      nforcest                 = 0
      ncount                   = ntrain
      nneshort_list(:)         = 0.0d0
      nnetot_list(:)           = 0.0d0
      edummy                   = 1.d12
      netot                    = 0
      nneshortsum              = 0.0d0
      nnstddev                 = 0.0d0
      nnstddevq(:)             = 0.0d0
      ncharge(:)               = 0
      avnncharge(:)            = 0.0d0
      ndummy                   = 0
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
          call readfunctions(1,symunit,npoints,nelem,&
            max_num_atoms,maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            symfunction_short_atomic_list)
        endif
!!
!! read the structures needed for the calculation of the electrostatic energy
!! and get reference forces from DFT
!! is needed for lshort (force fitting) and lelec (structure for electrostatics)
!! must be called after readfunctions because it needs num_atoms_list
        call getstructures(trainstructunit,npoints)
!!
      endif ! mpirank.eq.0
!!
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(symfunction_short_atomic_list,nblock*max_num_atoms*maxnum_funcvalues_short_atomic,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
!!
!! end of file reading for training error
!!--------------------------------------------------------------------
!!
!! prepare local input arrays for this process
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
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
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        shortenergy_mpi(icount)     = shortenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
      enddo ! i1
!!
!! get the short range energies for this block of points 
      if(lshort)then
        call getshortenergies_para(nstruct,ndummy,&
          ndummy,imaxerror_dummy,&
          num_atoms_mpi,zelem_mpi,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
          symfunction_short_atomic_list(1,1,n_start),lattice_list(1,1,n_start),&
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
!          write(ounit,*)' nneshort ',i1,nneshort_list(i1)
          nneshortsum=nneshortsum+nneshort_list(i1)
          nneshortmin=min(nneshortmin,nneshort_list(i1))
          nneshortmax=max(nneshortmax,nneshort_list(i1))
        endif
      enddo
!!
      deallocate(num_atoms_mpi)
      deallocate(zelem_mpi)
      deallocate(shortenergy_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(lperiodic_mpi)
      deallocate(nneshort_mpi)
      deallocate(nnshortforce_mpi)
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
      call mpi_allreduce(mpi_in_place,avnncharge,nelem,&
           mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,nenergies,1,&
           mpi_integer,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(mpi_in_place,ncharge,nelem,&
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
          call readfunctions(1,symunit,npoints,nelem,&
            max_num_atoms,maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            symfunction_short_atomic_list)
        endif
!!
!! read the structures needed for the calculation of the electrostatic energy
!! and get reference forces from DFT
!! is needed for lshort (force fitting) and lelec (structure for electrostatics)
!! must be called after readfunctions because it needs num_atoms_list
        call getstructures(trainstructunit,npoints)
!!
      endif ! mpirank.eq.0
!!
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(symfunction_short_atomic_list,nblock*max_num_atoms*maxnum_funcvalues_short_atomic,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
!!
!! end of file reading for training error
!!--------------------------------------------------------------------
!!
!! prepare local input arrays for this process
      allocate(num_atoms_mpi(nstruct))
      allocate(zelem_mpi(nstruct,max_num_atoms))
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
        num_atoms_mpi(icount)       = num_atoms_list(i1)
        zelem_mpi(icount,:)         = zelem_list(i1,:)
        shortenergy_mpi(icount)     = shortenergy_list(i1)
        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
        lperiodic_mpi(icount)       = lperiodic_list(i1)
      enddo ! i1
!!
!! get the short range energies for this block of points
      if(lshort)then
        call getshortenergies_para(nstruct,ndummy,&
          ndummy,imaxerror_dummy,&
          num_atoms_mpi,zelem_mpi,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
          symfunction_short_atomic_list(1,1,n_start),lattice_list(1,1,n_start),&
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
!! calculate part of the standard deviation   
        if(mpirank.eq.0)then
          do i1=1,npoints
            if(shortenergy_list(i1).le.maxenergy)then
!              write(ounit,*)i1,nneshort_list(i1),avnneshort
!              write(ounit,*)i1,(nneshort_list(i1)-avnneshort)**2.d0
              nnstddev=nnstddev+(nneshort_list(i1)-avnneshort)**2.d0
            endif
          enddo
        endif ! mpirank
      endif ! lshort
!!
      deallocate(num_atoms_mpi)
      deallocate(zelem_mpi)
      deallocate(shortenergy_mpi)
      deallocate(xyzstruct_mpi)
      deallocate(lperiodic_mpi)
      deallocate(nneshort_mpi)
      deallocate(nnshortforce_mpi)
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

!! scale the short range connecting weights between final hidden layer and output layer
          prefactor=eshortstddev/nnstddev
          write(ounit,'(a,f14.6)')' Factor for connecting short range weights: ',prefactor
          do i1=1,nelem !'
            do i2=windex_short_atomic(2*num_layers_short_atomic(i1)-1,i1),windex_short_atomic(2*num_layers_short_atomic(i1),i1)
              weights_short_atomic(i2,i1)=weights_short_atomic(i2,i1)*prefactor
            enddo
          enddo

!! shift the short range bias weights
          do i1=1,nelem
            weights_short_atomic(num_weights_short_atomic(i1),i1)&
              =weights_short_atomic(num_weights_short_atomic(i1),i1)-(prefactor*avnneshort-eshortav)   
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
