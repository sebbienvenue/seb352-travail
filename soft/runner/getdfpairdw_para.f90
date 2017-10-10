!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called once for each force component of each training point
!!
!! called by:
!! - optimize_short_combinedpair.f90 
!!

      subroutine getdfpairdw_para(iatom,ixyz,&
           num_weightspairfree,num_pairsp,pindex,&
           zelemp,num_pairs,wconstraintpidx,&
           symfunctionp,dsfuncdxyz_pair,dfpairdw)
!!
      use mpi_mod
      use fileunits
      use globaloptions
      use nnshort_pair
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5                                                ! internal
      integer num_weightspairfree(npairs)                                   ! in
      integer zelemp(2,max_num_pairs)                                       ! in
      integer num_pairs                                                     ! in
      integer iatom                                                         ! in
      integer ixyz                                                          ! in
      integer icount,jcount,kcount                                          ! internal
      integer itemp                                                         ! internal
      integer num_pairsp                                                    ! in
      integer pindex(num_pairsp)                                            ! in
      integer day                                                           ! internal
      integer wconstraintpidx(maxnum_weights_short_pair,npairs)                     ! in
!!
      real*8 symfunctionp(maxnum_funcvalues_short_pair,num_pairsp)                      ! in
      real*8 dfpairdw(maxnum_weights_short_pair,npairs)                            ! out 
      real*8 dfpairdw_temp(maxnum_weights_short_pair,num_pairs)                    ! internal 
      real*8 dnodes_values_dg(maxnum_layers_short_pair,maxnodes_short_pair,maxnum_funcvalues_short_pair)! internal
!! CAUTION: just one output node is assumed
      real*8 nnoutput                         ! internal 
      real*8 nodes_values(maxnum_layers_short_pair,maxnodes_short_pair)                   ! internal 
      real*8 nodes_sum(maxnum_layers_short_pair,maxnodes_short_pair)                    ! internal
!! CAUTION: just one output node is assumed
      real*8 dedw(maxnum_weights_short_pair,1)            ! internal
      real*8 dnodes_values(maxnum_layers_short_pair,maxnodes_short_pair)                  ! internal
      real*8 ddnodes_values(maxnum_layers_short_pair,maxnodes_short_pair)                 ! internal
!! CAUTION: just one output node is assumed
      real*8 dedgdw(maxnum_weights_short_pair,maxnum_funcvalues_short_pair,npairs,1) ! internal
      real*8 dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3)       ! in
      real*8 alphagaussian
!!
!!
!! initializations
      day=0
      dfpairdw(:,:)=0.0d0
      dfpairdw_temp(:,:)=0.0d0
      alphagaussian = 0.5d0
!!
!!      write(ounit,*)mpirank,' starts ',num_pairsp
!!
!! the loop over all atoms MUST be here!!! We have to sum for each force component again over all atoms
!! because the dedgdw are different for each atom
      do i1=1,num_pairsp
!!
!!
!! calculate nodes_values = y_i^j for atom i1
        nodes_values(:,:)=0.0d0
        nodes_sum(:,:)   =0.0d0
        itemp=pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1)))  
!!
        call calconenn(1,maxnum_funcvalues_short_pair,maxnodes_short_pair,&
          maxnum_layers_short_pair,num_layers_short_pair(itemp),&
          maxnum_weights_short_pair,nodes_short_pair(0,itemp),&
          symfunctionp(1,i1),weights_short_pair(1,itemp),nodes_values,nodes_sum,&
          nnoutput,actfunc_short_pair(1,1,itemp))
!!
!! calculate the first derivative at each node (\frac{\partial f_a(w)}{\partial G_i}) 
!! dnodes_values
        dnodes_values(:,:)=0.0d0
!!
        call getdnodes_values(maxnum_layers_short_pair,num_layers_short_pair,maxnodes_short_pair,&
          nodes_short_pair,itemp,npairs,nodes_sum,nodes_values,dnodes_values,actfunc_short_pair)
!!
!! calculate the second derivative at each node (\frac{\partial f_a(w)}{\partial G_i}) 
!! ddnodes_values
        ddnodes_values(:,:)=0.0d0
!!
        call getddnodes_values(maxnum_layers_short_pair,num_layers_short_pair,maxnodes_short_pair,&
          nodes_short_pair,itemp,npairs,nodes_sum,nodes_values,dnodes_values,&
          ddnodes_values,actfunc_short_pair)
!!
!! calculate dedw for \frac{\partial E_m}{\partial b}
!! this subroutine is calculating only!!! the derivatives for the bias weights
!! calculation of dedw is correct
        dedw(:,:)          =0.0d0
        call getonededb(&
          maxnum_weights_short_pair,&
          maxnodes_short_pair,maxnum_layers_short_pair,num_layers_short_pair(itemp),&
          windex_short_pair(1,itemp),nodes_short_pair(0,itemp),&
          weights_short_pair(1,itemp),dnodes_values,dedw)
!!
!! calculate dnodes_values_dg = \frac{\partial x_k^j}{\partial G_i}
!! calculation of dnodes_values_dg is correct
        dnodes_values_dg(:,:,:)=0.0d0
        call get_dnodes_values_dg(&
          maxnum_layers_short_pair,num_layers_short_pair(itemp),maxnodes_short_pair,maxnum_funcvalues_short_pair,&
          maxnum_weights_short_pair,windex_short_pair(1,itemp),nodes_short_pair(0,itemp),&
          weights_short_pair(1,itemp),dnodes_values,dnodes_values_dg)
!!
!! calculate dedgdw
        dedgdw(:,:,:,:)=0.0d0

!! calculate bias weight derivatives \frac{\partial}{\partial G_i}\frac{\partial E_m}{\partial b_{\mu}^{\nu}}
!! dedgdw for bias weights
!!        write(ounit,*)'calculating dedgdw for bias'
        do i2=num_layers_short_pair(itemp),1,-1
          if(i2.eq.num_layers_short_pair(itemp))then ! output layer, we assume just one output node here
            do i3=1,nodes_short_pair(num_layers_short_pair(itemp),itemp) ! this is just one output node
              icount=windex_short_pair(2*i2,itemp)+i3-1  ! pointer to bias weight for output node
              dedgdw(icount,:,itemp,1)=ddnodes_values(i2,i3)*dnodes_values_dg(i2,i3,:)
            enddo
          else ! i2.lt.num_layerspair => hidden layer
            do i3=1,nodes_short_pair(i2,itemp) ! one bias weight derivative for each node in previous hidden layer 
              do i4=1,nodes_short_pair(i2+1,itemp) ! nodes in target layer
                icount=windex_short_pair(2*i2,itemp)+i3-1                            ! pointer to bias weight in previous layer 
                jcount=windex_short_pair(2*i2+1,itemp)+nodes_short_pair(i2+1,itemp)*(i3-1)+i4-1 ! pointer to connecting weight
                kcount=windex_short_pair(2*(i2+1),itemp)+i4-1                        ! pointer to bias weight in target layer
                dedgdw(icount,:,itemp,1)=dedgdw(icount,:,itemp,1) &
                + weights_short_pair(jcount,itemp)*dnodes_values(i2,i3)*dedgdw(kcount,:,itemp,1) &
                + dedw(kcount,1)*weights_short_pair(jcount,itemp)*ddnodes_values(i2,i3)*dnodes_values_dg(i2,i3,:)
              enddo ! i4
            enddo ! i3
          endif
!!
        enddo ! i2
!!
!! calculate connection weight derivatives \frac{\partial}{\partial G_i}\frac{\partial E_m}{\partial a_{ij}^{kl}}
!! dedgdw for connection weights
        do i2=1,num_layers_short_pair(itemp) ! loop over target layers
          if(i2.eq.1)then ! 1st hidden layer
            do i3=1,nodes_short_pair(i2,itemp) ! over all nodes in 1st target hidden layer
              do i4=1,num_funcvalues_short_pair(itemp) ! over all input nodes
                icount=windex_short_pair(2*i2-1,itemp)+nodes_short_pair(i2,itemp)*(i4-1)+i3-1   ! pointer to connecting weight
                kcount=windex_short_pair(2*i2,itemp)+i3-1                            ! pointer to bias weight in target layer
                dedgdw(icount,i4,itemp,1)=symfunctionp(i4,i1)*dedgdw(kcount,i4,itemp,1)*dedw(kcount,1)
              enddo ! i4
            enddo ! i3
          else ! i2.gt.1 => not first hidden layer
            do i3=1,nodes_short_pair(i2,itemp)     ! over all nodes in target layer
              do i4=1,nodes_short_pair(i2-1,itemp) ! over all nodes in previous layer
                do i5=1,num_funcvalues_short_pair(itemp)  ! over all input nodes
            icount=windex_short_pair(2*i2-1,itemp)+nodes_short_pair(i2,itemp)*(i4-1)+i3-1       ! pointer to connection weight
            kcount=windex_short_pair(2*i2,itemp)+i3-1                                ! pointer to bias weight in target layer
            dedgdw(icount,i5,itemp,1)=nodes_values(i2-1,i4)*dedgdw(kcount,i5,itemp,1)&
              +dedw(kcount,1)*dnodes_values(i2-1,i4)*dnodes_values_dg(i2-1,i4,i5)
                enddo ! i5
              enddo ! i4
            enddo ! i3
          endif
        enddo ! i2
!!
!! calculate derivative of force with respect to weights \frac{\partial F_{\alpha}}{\partial w}
!! dfpairdw
!! calculate dfpairdw for atom iatom and component ixyz
!! CAUTION: This is critical because in parallel jobs the order of summation matters!!!
        do i2=1,num_funcvalues_short_pair(itemp)
          do i3=1,num_weightspairfree(pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1)))) 
            dfpairdw(i3,itemp)=dfpairdw(i3,itemp) &
            + dsfuncdxyz_pair(i2,pindex(i1),iatom,ixyz)&
            *dedgdw(wconstraintpidx(i3,pairindex(zelemp(1,pindex(i1)),zelemp(2,pindex(i1)))),i2,itemp,1)
          enddo
        enddo ! i2
!!
      enddo ! i1 loop over pairs 
!!
      call mpi_allreduce(mpi_in_place,dfpairdw,maxnum_weights_short_pair*npairs,&
       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
      return
      end

