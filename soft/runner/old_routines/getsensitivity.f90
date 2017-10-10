!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! this is a multipurpose subroutine

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine getsensitivity(ndim1,ndim2,ntrain,train_local,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        minvalue_local,maxvalue_local,avvalue_local,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        scmin_local,scmax_local,sens,sense)
!!
      use mpi_mod
      use globaloptions
      use fileunits
      use symfunctions
      use nnshort_atomic
      use nnewald
      use structures
!!
      implicit none
!!
      integer ndim1                                                      ! nelem or npairs
      integer ndim2                                                      ! max_num_atoms or max_num_pairs
      integer maxnum_funcvalues_local                                    ! in
      integer num_funcvalues_local(ndim1)                                ! in
      integer i1,i2,i3                                                   ! internal
      integer ntrain                                                     ! in
      integer ndone                                                      ! internal
      integer ncount                                                     ! internal
      integer npoints                                                    ! internal
      integer num_atoms                                                  ! internal
      integer num_pairs                                                  ! internal
      integer zelem(max_num_atoms)                                       ! internal
      integer zelemp(2,max_num_pairs)                                    ! internal
      integer train_local(nelem)                                           ! in
!!
      real*8 sens(ndim1,maxnum_funcvalues_local)                         ! out
      real*8 sense(nelem,maxnum_funcvalues_elec)                             ! out
      real*8 symfunction(maxnum_funcvalues_local,ndim2)                  ! internal
      real*8 symfunction_local(maxnum_funcvalues_local,ndim2,nblock)     ! internal
      real*8 symfunctione(maxnum_funcvalues_elec,max_num_atoms)              ! internal
      real*8 deshortdsfunc(max_num_atoms,maxnum_funcvalues_local)        ! internal
      real*8 depairdsfunc(max_num_pairs,maxnum_funcvalues_local)         ! internal
      real*8 dchargedsfunc(max_num_atoms,maxnum_funcvalues_short_atomic)              ! internal
      real*8 minvalue_local(ndim1,maxnum_funcvalues_local)               ! in
      real*8 maxvalue_local(ndim1,maxnum_funcvalues_local)               ! in
      real*8 avvalue_local(ndim1,maxnum_funcvalues_local)                ! in
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)               ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)               ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)                ! in
      real*8 scmin_local                                                 ! in
      real*8 scmax_local                                                 ! in 
!!
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
      sens(:,:) = 0.0d0
      sense(:,:)= 0.0d0
      ncount    = ntrain
      ndone     = 0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! open files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mpirank.eq.0)then
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
        if(lshort)then
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit)
        endif ! lshort
        if(lelec.and.(nn_type_elec.eq.1))then
          open(symeunit,file='functione.data',form='formatted',status='old')
          rewind(symeunit) !'
        endif ! lelec
      endif ! mpirank.eq.0
!!
!! loop block-wise over all points in training set
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!!--------------------------------------------------------------------------------------
!! do all file reading for this block of points 
!!--------------------------------------------------------------------------------------
      if(mpirank.eq.0)then

!! FIXME ANDI: In case of scale_symmetry_functions the array symfunction_list must contain 
        if(lshort)then
!! read npoint short range data sets (train or test)
          if(nn_type_short.eq.1)then
            call readfunctions(1,symunit,npoints,nelem,&
              ndim2,maxnum_funcvalues_local,num_funcvalues_local,&
              symfunction_local)
!! scale symmetry functions for the short-range interaction
            call scalesymfit_para(1,nelem,npoints,1,npoints,&
              maxnum_funcvalues_local,num_funcvalues_local,&
              minvalue_local,maxvalue_local,avvalue_local,symfunction_local,&
              scmin_local,scmax_local)
          elseif(nn_type_short.eq.2)then
            call readfunctions(2,symunit,npoints,npairs,&
              ndim2,maxnum_funcvalues_local,num_funcvalues_local,&
              symfunction_local)
!! scale symmetry functions for the short-range interaction
            call scalesymfit_para(2,npairs,npoints,1,npoints,&
              maxnum_funcvalues_local,num_funcvalues_local,&
              minvalue_local,maxvalue_local,avvalue_local,symfunction_local,&
              scmin_local,scmax_local)
          else
            write(ounit,*)'ERROR in getsensitivity, unknown nn_type_short ',nn_type_short
            stop
          endif
        endif
        if(lelec.and.(nn_type_elec.eq.1))then
!! read the symmetry functions for the charge prediction (train or test)
          call readfunctions(1,symeunit,npoints,nelem,&
            max_num_atoms,maxnum_funcvalues_elec,num_funcvalues_elec,&
            symfunction_elec_list)
!! scale symmetry functions for the short-range interaction
          call scalesymfit_para(3,nelem,npoints,1,npoints,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            minvalue_elec,maxvalue_elec,avvalue_elec,symfunction_elec_list,&
            scmin_elec,scmax_elec)
        endif
!! END FIXME

!! read the structures
        call getstructures(trainstructunit,npoints)
      endif ! mpirank.eq.0
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_local,nblock*ndim2*maxnum_funcvalues_local,&
           mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_elec_list,nblock*max_num_atoms*maxnum_funcvalues_elec,&
           mpi_real8,0,mpi_comm_world,mpierror)
!!--------------------------------------------------------------------
!! end of file reading 
!!--------------------------------------------------------------------
!!
!!--------------------------------------------------------------------
!! calculate deshortdsfunc for all structures in block 
!!--------------------------------------------------------------------
      do i1=1,npoints
        num_atoms        =num_atoms_list(i1)
        num_pairs        =num_pairs_list(i1)
        zelem(:)         =zelem_list(i1,:)
        zelemp(:,:)      =zelemp_list(:,i1,:)
        symfunction(:,:) =symfunction_local(:,:,i1)
        symfunctione(:,:)=symfunction_elec_list(:,:,i1)
!!
        if(lshort)then
!! get deshortdsfunc
          if(nn_type_short.eq.1)then
            call getdeshortdsfunc(num_atoms,zelem,symfunction,&
              deshortdsfunc)
            do i2=1,num_atoms
              do i3=1,num_funcvalues_local(elementindex(zelem(i2)))

!!
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!
                !sens(elementindex(zelem(i2)),i3)&
                !=sens(elementindex(zelem(i2)),i3)+deshortdsfunc(i2,i3)
                sens(elementindex(zelem(i2)),i3)&
                =sens(elementindex(zelem(i2)),i3)+deshortdsfunc(i2,i3)**2.0d0
!! END CHANGE

              enddo ! i3
            enddo ! i2
          elseif(nn_type_short.eq.2)then
            call getdepairdsfunc(num_pairs,zelemp,symfunction,depairdsfunc)
            do i2=1,num_pairs
              do i3=1,num_funcvalues_local(pairindex(zelemp(1,i2),zelemp(2,i2)))

!!
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!
                !sens(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)&
                !=sens(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)+depairdsfunc(i2,i3)
                sens(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)&
                =sens(pairindex(zelemp(1,i2),zelemp(2,i2)),i3)+depairdsfunc(i2,i3)**2.0d0
!! END CHANGE
              enddo ! i3
            enddo ! i2
          else
            write(ounit,*)'ERROR in getsensitivity, unknown nn_type_short ',nn_type_short
            stop
          endif
        endif ! lshort
!!
        if(lelec.and.(nn_type_elec.eq.1))then 
!! get dchargedsfunc
          call getdchargedsfunc(num_atoms,zelem,symfunctione,&
            dchargedsfunc)
          do i2=1,num_atoms
            do i3=1,num_funcvalues_elec(elementindex(zelem(i2)))

!!
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!
             !sense(elementindex(zelem(i2)),i3)&
             !=sense(elementindex(zelem(i2)),i3)+dchargedsfunc(i2,i3)
              sense(elementindex(zelem(i2)),i3)&
              =sense(elementindex(zelem(i2)),i3)+dchargedsfunc(i2,i3)**2.0d0
!! END CHANGE

            enddo ! i3
          enddo ! i2
        endif ! lelec
      enddo ! i1
!!
      ndone=ndone+npoints
      if(ncount.gt.0) goto 10
!! end block of training points
!!
!! normalize
      if(lshort)then
        do i1=1,nelem

!!
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!
         !sens(i1,:)=sens(i1,:)/dble(train_local(i1))
          sens(i1,:)=dsqrt(sens(i1,:)/dble(train_local(i1)))
!! END CHANGE

        enddo
      endif
!!
!! FIXME: train_local is trainpair in case of nn_type_short=2, this does not make sense for electrostatics
      if(lelec.and.(nn_type_elec.eq.1))then 
        do i1=1,nelem

!!
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!
         !sense(i1,:)=sense(i1,:)/dble(train_local(i1))
          sense(i1,:)=dsqrt(sense(i1,:)/dble(train_local(i1)))
!! END CHANGE

        enddo
      endif
!!
      return
      end
