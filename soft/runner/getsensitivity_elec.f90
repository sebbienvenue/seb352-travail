!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################


!! called by:
!!
      subroutine getsensitivity_elec(ntrain,trainelem,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        sense)
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
      integer i1,i2,i3                                                   ! internal
      integer ntrain                                                     ! in
      integer ndone                                                      ! internal
      integer ncount                                                     ! internal
      integer npoints                                                    ! internal
      integer num_atoms                                                  ! internal
      integer num_pairs                                                  ! internal
      integer zelem(max_num_atoms)                                       ! internal
      integer trainelem(nelem)                                           ! in
!!
      real*8 sense(nelem,maxnum_funcvalues_elec)                             ! out
      real*8 symfunctione(maxnum_funcvalues_elec,max_num_atoms)              ! internal
      real*8 dchargedsfunc(max_num_atoms,maxnum_funcvalues_short_atomic)              ! internal
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)               ! in
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)               ! in
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)                ! in
!!
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
        open(symeunit,file='functione.data',form='formatted',status='old')
        rewind(symeunit) !'
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
!! read the symmetry functions for the charge prediction (train or test)
        call readfunctions(1,symeunit,npoints,nelem,&
          max_num_atoms,maxnum_funcvalues_elec,num_funcvalues_elec,&
          symfunction_elec_list)
!! scale symmetry functions for the short-range interaction
        call scalesymfit_para(3,nelem,npoints,1,npoints,&
          maxnum_funcvalues_elec,num_funcvalues_elec,&
          minvalue_elec,maxvalue_elec,avvalue_elec,symfunction_elec_list,&
          scmin_elec,scmax_elec)
!! read the structures
        call getstructures(trainstructunit,npoints)
      endif ! mpirank.eq.0
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
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
        symfunctione(:,:)=symfunction_elec_list(:,:,i1)
!!
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
      enddo ! i1
!!
      ndone=ndone+npoints
      if(ncount.gt.0) goto 10
!! end block of training points
!!
      do i1=1,nelem

!!
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!
       !sense(i1,:)=sense(i1,:)/dble(trainelem(i1))
        sense(i1,:)=dsqrt(sense(i1,:)/dble(trainelem(i1)))
!! END CHANGE
!!
      enddo
!!
      return
      end
