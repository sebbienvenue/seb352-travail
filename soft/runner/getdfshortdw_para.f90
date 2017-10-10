!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called once for each individual force component of each training point.
!! calculates the derivative of the force with respect to the weights.
!!
!! called by:
!! - optimize_short_combined.f90
!!
      subroutine getdfshortdw_para(iatom,ixyz,max_num_neighbors_atomic,&
           num_weights_short_atomic_free,natoms,atomindex,&
           zelem,num_atoms,wconstraintidx,&
           symfunction,dsfuncdxyz,dfshortdw)
!!
      use mpi_mod
      use fileunits
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5                                                ! internal
      integer max_num_neighbors_atomic                                      ! in
      integer num_weights_short_atomic_free(nelem)                          ! in
      integer zelem(max_num_atoms)                                          ! in 
      integer num_atoms                                                     ! in
      integer iatom                                                         ! in
      integer ixyz                                                          ! in
      integer icount,jcount,kcount                                          ! internal
      integer itemp                                                         ! internal
      integer natoms                                                        ! in
      integer atomindex(natoms)                                             ! in
      integer day
      integer wconstraintidx(maxnum_weights_short_atomic,nelem)                    ! in
!!
      real*8 symfunction(maxnum_funcvalues_short_atomic,natoms)                      ! in
      real*8 dfshortdw(maxnum_weights_short_atomic,nelem)                              ! out 
      real*8 dfshortdw_temp(maxnum_weights_short_atomic,num_atoms)                     ! internal 
      real*8 dnodes_values_dg(maxnum_layers_short_atomic,maxnodes_short_atomic,maxnum_funcvalues_short_atomic)! internal
!! CAUTION: just one output node is assumed
      real*8 nnoutput                         ! internal 
      real*8 nodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic)                   ! internal 
      real*8 nodes_sum(maxnum_layers_short_atomic,maxnodes_short_atomic)                    ! internal
!! CAUTION: just one output node is assumed
      real*8 dedw(maxnum_weights_short_atomic,1)            ! internal
      real*8 dnodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic)                  ! internal
      real*8 ddnodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic)                 ! internal
!! CAUTION: just one output node is assumed
      real*8 dedgdw(maxnum_weights_short_atomic,maxnum_funcvalues_short_atomic,nelem,1) ! internal
      real*8 dsfuncdxyz(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_atomic,3)       ! in
      real*8 alphagaussian
!!
!!
!! initializations
      day=0
      dfshortdw(:,:)=0.0d0
      dfshortdw_temp(:,:)=0.0d0
      alphagaussian = 0.5d0
!!
!!      write(ounit,*)mpirank,' starts ',natoms
!!
!! the loop over all atoms MUST be here!!! We have to sum for each force component again over all atoms
!! because the dedgdw are different for each atom
      do i1=1,natoms
!!
!!
!! calculate nodes_values = y_i^j for atom i1
        nodes_values(:,:)=0.0d0
        nodes_sum(:,:)   =0.0d0
        itemp=elementindex(zelem(atomindex(i1)))
!!
!! calculate nodes_sum and nodes_values for atom i1 
        call calconenn(1,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic(itemp),&
          maxnum_weights_short_atomic,nodes_short_atomic(0,itemp),&
          symfunction(1,i1),weights_short_atomic(1,itemp),nodes_values,nodes_sum,&
          nnoutput,actfunc_short_atomic(1,1,itemp))
!!
!! calculate the first derivative at each node (\frac{\partial f_a(w)}{\partial G_i}) 
!! dnodes_values
        dnodes_values(:,:)=0.0d0
        call getdnodes_values(maxnum_layers_short_atomic,num_layers_short_atomic,maxnodes_short_atomic,&
          nodes_short_atomic,itemp,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_short_atomic)
!!
!! calculate the second derivative at each node (\frac{\partial f_a(w)}{\partial G_i}) 
!! ddnodes_values
        ddnodes_values(:,:)=0.0d0
        call getddnodes_values(maxnum_layers_short_atomic,num_layers_short_atomic,maxnodes_short_atomic,&
          nodes_short_atomic,itemp,nelem,nodes_sum,nodes_values,dnodes_values,&
          ddnodes_values,actfunc_short_atomic)
!!
!! calculate dedw for \frac{\partial E_m}{\partial b}
!! this subroutine is calculating only!!! the derivatives for the bias weights
!! calculation of dedw is correct
        dedw(:,:)          =0.0d0
        call getonededb(&
          maxnum_weights_short_atomic,&
          maxnodes_short_atomic,maxnum_layers_short_atomic,num_layers_short_atomic(itemp),&
          windex_short_atomic(1,itemp),nodes_short_atomic(0,itemp),&
          weights_short_atomic(1,itemp),dnodes_values,dedw)
!!
!! calculate dnodes_values_dg = \frac{\partial x_k^j}{\partial G_i}
!! calculation of dnodes_values_dg is correct
        dnodes_values_dg(:,:,:)=0.0d0
        call get_dnodes_values_dg(&
          maxnum_layers_short_atomic,num_layers_short_atomic(itemp),maxnodes_short_atomic,&
          maxnum_funcvalues_short_atomic,&
          maxnum_weights_short_atomic,windex_short_atomic(1,itemp),nodes_short_atomic(0,itemp),&
          weights_short_atomic(1,itemp),dnodes_values,dnodes_values_dg)
!!
!! calculate dedgdw
        dedgdw(:,:,:,:)=0.0d0
!!
!! calculate bias weight derivatives \frac{\partial}{\partial G_i}\frac{\partial E_m}{\partial b_{\mu}^{\nu}}
!! dedgdw for bias weights
!!        write(ounit,*)'calculating dedgdw for bias'
        do i2=num_layers_short_atomic(itemp),1,-1
          icount=windex_short_atomic(2*i2,itemp)-1  ! pointer to bias weight for output node
          if(i2.eq.num_layers_short_atomic(itemp))then ! output layer, we assume just one output node here
            do i3=1,nodes_short_atomic(num_layers_short_atomic(itemp),itemp) ! this is just one output node
              icount=icount+1  ! pointer to bias weight for output node
              dedgdw(icount,:,itemp,1)=ddnodes_values(i2,i3)*dnodes_values_dg(i2,i3,:)
            enddo
          else ! i2.lt.num_layers_short_atomic => hidden layer
            do i3=1,nodes_short_atomic(i2,itemp) ! one bias weight derivative for each node in previous hidden layer 
              kcount=windex_short_atomic(2*(i2+1),itemp)-1                        ! pointer to bias weight in target layer
              icount=icount+1                            ! pointer to bias weight in previous layer 
              jcount=windex_short_atomic(2*i2+1,itemp)+nodes_short_atomic(i2+1,itemp)*(i3-1)-1 ! pointer to connecting weight
              do i4=1,nodes_short_atomic(i2+1,itemp) ! nodes in target layer
                jcount=jcount+1
                kcount=kcount+1
                dedgdw(icount,:,itemp,1)=dedgdw(icount,:,itemp,1) &
                + weights_short_atomic(jcount,itemp)*dnodes_values(i2,i3)*dedgdw(kcount,:,itemp,1) &
                + dedw(kcount,1)*weights_short_atomic(jcount,itemp)*ddnodes_values(i2,i3)*dnodes_values_dg(i2,i3,:)
              enddo ! i4
            enddo ! i3
          endif
        enddo ! i2
!!
!! calculate connection weight derivatives \frac{\partial}{\partial G_i}\frac{\partial E_m}{\partial a_{ij}^{kl}}
!! dedgdw for connection weights
        do i2=1,num_layers_short_atomic(itemp) ! loop over target layers
          if(i2.eq.1)then ! 1st hidden layer
            do i4=1,num_funcvalues_short_atomic(itemp) ! over all input nodes
              kcount=windex_short_atomic(2*i2,itemp)-1                            ! pointer to bias weight in target layer
              icount=windex_short_atomic(2*i2-1,itemp)+nodes_short_atomic(i2,itemp)*(i4-1)-1   ! pointer to connecting weight
              do i3=1,nodes_short_atomic(i2,itemp) ! over all nodes in 1st target hidden layer
                kcount=kcount+1
                icount=icount+1
!! bug detected by Andi
!! old wrong line:
!!                dedgdw(icount,i4,itemp,1)=symfunction(i4,i1)*dedgdw(kcount,i4,itemp,1)+dedw(kcount,1) 
!! correct:
                do i5=1,num_funcvalues_short_atomic(itemp) ! over all input nodes
                  dedgdw(icount,i5,itemp,1)=symfunction(i4,i1)*dedgdw(kcount,i5,itemp,1) 
                  if(i4.eq.i5)then
                    dedgdw(icount,i5,itemp,1)=dedgdw(icount,i5,itemp,1)+dedw(kcount,1) 
                  endif
                enddo ! i5
              enddo ! i4
            enddo ! i3
          else ! i2.gt.1 => not first hidden layer
            do i4=1,nodes_short_atomic(i2-1,itemp) ! over all nodes in previous layer
              kcount=windex_short_atomic(2*i2,itemp)-1                                ! pointer to bias weight in target layer
              icount=windex_short_atomic(2*i2-1,itemp)+nodes_short_atomic(i2,itemp)*(i4-1)-1       ! pointer to connection weight
              do i3=1,nodes_short_atomic(i2,itemp)     ! over all nodes in target layer
                icount=icount+1
                kcount=kcount+1                                ! pointer to bias weight in target layer
                do i5=1,num_funcvalues_short_atomic(itemp)  ! over all input nodes
                  dedgdw(icount,i5,itemp,1)=nodes_values(i2-1,i4)*dedgdw(kcount,i5,itemp,1)&
                  +dedw(kcount,1)*dnodes_values(i2-1,i4)*dnodes_values_dg(i2-1,i4,i5)
                enddo ! i5
              enddo ! i3
            enddo ! i4
          endif
        enddo ! i2
!!
!! calculate derivative of force with respect to weights \frac{\partial F_{\alpha}}{\partial w}
!! dfshortdw
!! calculate dfshortdw for atom iatom and component ixyz
!! CAUTION: This is critical because in parallel jobs the order of summation matters!!!
        do i2=1,num_funcvalues_short_atomic(itemp)
          do i3=1,num_weights_short_atomic_free(itemp)
            dfshortdw(i3,itemp)=dfshortdw(i3,itemp) &
            + dsfuncdxyz(i2,atomindex(i1),iatom,ixyz)*dedgdw(wconstraintidx(i3,itemp),i2,itemp,1)
          enddo
        enddo ! i2
!!
!! debug
!!
!! alternative (part1), which is hopefully independent of the number of processes
!! CAUTION: With the alternative the results are independent of the number of processors, but different, 
!! probably because sum summation of dedgdw is messed up
!!        do i2=1,num_funcvalues
!!          do i3=1,num_weightsshortfree
!!          dfshortdw_temp(i3,atomindex(i1))=dsfuncdxyz(i2,atomindex(i1),iatom,ixyz)*dedgdw(wconstraintidx(i3),i2,1)
!!          enddo ! i3
!!        enddo ! i2
!! end alternative (part1)
!!
      enddo ! i1 loop over atoms
!!
!! alternative (part2), which is hopefully independent of the number of processes
!! now we do the summation for ALL atoms here in a unique way
!!      call mpi_reduce(dfshortdw_temp,dfshortdw_temp,maxnum_weightsshort*num_atoms,&
!!       mpi_real8,mpi_sum,0,mpi_comm_world,mpierror)
!!      call mpi_allreduce(mpi_in_place,dfshortdw_temp,maxnum_weightsshort*num_atoms,&
!!       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!      if(mpirank.eq.0)then
!!        do i1=1,num_atoms
!!          dfshortdw(:)=dfshortdw(:)+dfshortdw_temp(:,i1)
!!        enddo ! i1
!!      endif
!!      call mpi_bcast(dfshortdw,maxnum_weightsshort,mpi_real8,0,mpi_comm_world,mpierror)
!! end alternative (part2)

!!      call mpi_barrier(mpi_comm_world,mpierror)
!!      stop
!! combine all processes
!!
      call mpi_allreduce(mpi_in_place,dfshortdw,maxnum_weights_short_atomic*nelem,&
       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! CAUTION: The alternative part has a serious bug causing explosions of the fit
!!
      return
      end
