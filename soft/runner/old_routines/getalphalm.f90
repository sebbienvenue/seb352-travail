!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
!! This routine is called once for each block of points
!!
      subroutine getalphalm(nblock,npoints,dimlm,&
         nelem,maxnum_weightsewald,num_weightsewald,max_num_atoms,&
         maxnum_funcvaluese,num_funcvaluese,&
         maxnum_layersewald,num_layersewald,windexe,zelem_list,maxnodes_ewald,&
         num_atoms_list,nodes_ewald,elementindex,&
         symfunctione_list,weights_ewald,stddevcharge,&
         atomcharge_list,alphalm,betalm,chisq,&
         actfunc_ewald,&
         lnormnodes,ldebug)
!!
      use mpi_mod
      use fileunits
!!
      implicit none
!!
      integer elementindex(102)                     ! in
      integer i0,i1,i2,i3,i4                        ! internal
      integer max_num_atoms                         ! in
      integer maxnodes_ewald                        ! in
      integer nblock                                ! in
      integer nelem                                 ! in
      integer maxnum_layersewald                    ! in
      integer num_layersewald(nelem)                ! in
      integer nodes_ewald(0:maxnum_layersewald,nelem)     ! in
      integer npoints                               ! in
      integer num_atoms                             ! internal
      integer num_atoms_list(nblock)                ! in
      integer maxnum_funcvaluese                       ! in
      integer num_funcvaluese(nelem)                       ! in
      integer maxnum_weightsewald                      ! in
      integer num_weightsewald(nelem)                      ! in
      integer windexe(2*maxnum_layersewald,nelem)            ! in
      integer zelem_list(nblock,max_num_atoms)      ! in 
      integer dimlm                                 ! in
!!
      real*8 error                                                          ! internal
      real*8 weights_ewald(maxnum_weightsewald,nelem)                          ! in/out
      real*8 symfunctione_list(maxnum_funcvaluese,max_num_atoms,nblock)        ! in
!! CAUTION: just one output node is assumed here
      real*8 dedw(maxnum_weightsewald,1)            ! internal
      real*8 nnatomcharge                                                   ! internal
      real*8 atomcharge_list(nblock,max_num_atoms)                          ! in 
      real*8 nodes_values_dummy(maxnum_layersewald,maxnodes_ewald)          ! just dummy here 
      real*8 nodes_sum_dummy(maxnum_layersewald,maxnodes_ewald)             ! just dummy here 
      real*8 nnoutput                                                       ! internal
      real*8 alphalm(dimlm,dimlm,nelem)                                     ! in/out
      real*8 betalm(dimlm,nelem)                                            ! in/out
      real*8 temp                                                           ! internal
      real*8 stddevinv                                                      ! in
      real*8 chisq(nelem)                                                   ! in/out
      real*8 sigma
      real*8 stddevcharge(nelem)
!!
      character*1 actfunc_ewald(maxnodes_ewald,maxnum_layersewald,nelem)                            ! in
!!
      logical ldebug                                                        ! in
      logical lnormnodes                                                    ! in
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop over all points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i1=1,npoints
!!
!! initializations for this training point
        num_atoms           =num_atoms_list(i1)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop over all atoms i2 of point i1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i2=1,num_atoms
!!
!! predict the NN charge of that atom
          call calconeatom(maxnum_funcvaluese,maxnodes_ewald,&
            maxnum_layersewald,num_layersewald(elementindex(zelem_list(i1,i2))),&
            maxnum_weightsewald,nodes_ewald(0,elementindex(zelem_list(i1,i2))),&
            symfunctione_list(1,i2,i1),&
            weights_ewald(1,elementindex(zelem_list(i1,i2))),&
            nodes_values_dummy,nodes_sum_dummy,&
            nnoutput,actfunc_ewald(1,1,elementindex(zelem_list(i1,i2))),lnormnodes,ldebug)
!!         
          nnatomcharge=nnoutput
!!
          error = atomcharge_list(i1,i2)-nnatomcharge
          sigma = atomcharge_list(i1,i2)*stddevcharge(elementindex(zelem_list(i1,i2)))
          stddevinv = 1.d0/(sigma*sigma)
!!
!! calculate the derivate with respect to the weights
          dedw(:,:)=0.0d0 ! initialization for this atom
          call getonedeshortdw(&
            maxnum_funcvaluese,maxnum_weightsewald,&
            maxnodes_ewald,maxnum_layersewald,num_layersewald(elementindex(zelem_list(i1,i2))),&
            windexe(1,elementindex(zelem_list(i1,i2))),nodes_ewald(0,elementindex(zelem_list(i1,i2))),&
            symfunctione_list(1,i2,i1),&
            weights_ewald(1,elementindex(zelem_list(i1,i2))),&
            dedw,&
            actfunc_ewald(1,1,elementindex(zelem_list(i1,i2))),&
            ldebug)
!!
          do i3=1,num_weightsewald(elementindex(zelem_list(i1,i2)))
!!
            temp=dedw(i3,1)*stddevinv
!!
            do i4=1,i3
!!
              alphalm(i3,i4,elementindex(zelem_list(i1,i2)))&
                =alphalm(i3,i4,elementindex(zelem_list(i1,i2)))+temp*dedw(i4,1)
!!
            enddo ! i4
!!
            betalm(i3,elementindex(zelem_list(i1,i2)))&
              =betalm(i3,elementindex(zelem_list(i1,i2)))+error*temp
!!
          enddo ! i3
!!
          chisq(elementindex(zelem_list(i1,i2)))&
            =chisq(elementindex(zelem_list(i1,i2)))+error*error*stddevinv
!!
        enddo ! i2, loop over all atoms
!!
      enddo ! i1, loop over all points
!!
      return
      end
