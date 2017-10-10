!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate nnatomcharge() and sense()

!! called by: 
!!
      subroutine getchargesatomic(n_start,natoms,atomindex,max_num_neighbors_elec,&
        invneighboridx_elec,zelem,&
        lsta,lstc,lste,lstb,xyzstruct,&
        sense,minvalue_local,maxvalue_local,&
        avvalue_local,scmin_local,scmax_local,nnatomcharge,&
        lextrapolation)
!!
      use fileunits
      use globaloptions
      use nnewald
      use symfunctions
      use timings
      use predictionoptions
!!
      implicit none
!!
      integer n_start                                                   ! in
      integer max_num_neighbors_elec                                    ! in
      integer ielem                                                     ! internal
      integer iindex                                                    ! internal
      integer zelem(max_num_atoms)                                      ! in
      integer i1,i2,i3,i4,i5                                            ! internal
      integer icount                                                    ! internal
      integer jcount                                                    ! internal
      integer natoms                                                    ! in
      integer atomindex(natoms)                                         ! in
      integer invneighboridx_elec(natoms,max_num_atoms)                 ! in
      integer lsta(2,max_num_atoms)                                     ! in
      integer lstc(listdim)                                             ! in, identification of atom
      integer lste(listdim)                                             ! in, nuclear charge of atom
!!
      real*8 xyzstruct(3,max_num_atoms)                                 ! in
      real*8 lstb(listdim,4)                                            ! in, xyz and r_ij
      real*8 minvalue_local(nelem,maxnum_funcvalues_elec)                   ! in 
      real*8 maxvalue_local(nelem,maxnum_funcvalues_elec)                   ! in 
      real*8 avvalue_local(nelem,maxnum_funcvalues_elec)                    ! in 
      real*8 scmin_local                                                ! in
      real*8 scmax_local                                                ! in
      real*8 strs_dummy(3,3,maxnum_funcvalues_elec) 
      real*8 dsfuncdxyz_dummy(0:max_num_neighbors_elec,3)               ! internal
      real*8 dchargedsfunc(maxnum_funcvalues_elec)                          ! internal 
      real*8 symfunctione(maxnum_funcvalues_elec)                           ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                                   ! internal 
      real*8 nodes_values(maxnum_layers_elec,maxnodes_elec)            ! internal
      real*8 nodes_sum(maxnum_layers_elec,maxnodes_elec)               ! internal
      real*8 dnodes_values(maxnum_layers_elec,maxnodes_elec)           ! internal
      real*8 tempderivative(maxnum_layers_elec,maxnodes_elec,maxnum_funcvalues_elec) ! internal
      real*8 nnatomcharge(max_num_atoms)                                ! out 
      real*8 threshold                                                  ! internal
      real*8 sense(nelem,maxnum_funcvalues_elec)                            ! out 
!!  
      logical lrmin                                                     ! internal 
      logical lextrapolation                                            ! out
!!
!! initialization
      lrmin           =.true.
      threshold       = 0.0001d0
!!
!!====================================================================================
!! loop over all atoms 
!!====================================================================================
      jcount=n_start
      do i1=1,natoms
!!====================================================================================
!! initializations for this atom 
!!====================================================================================
        dchargedsfunc(:)= 0.0d0
        symfunctione(:)  = 0.0d0
        ielem=elementindex(zelem(atomindex(i1)))
        iindex=elementindex(zelem(jcount))
!!====================================================================================
!! get the symmetry functions for atom jcount 
!!====================================================================================
        do i2=1,num_funcvalues_elec(ielem) ! over all symmetry functions
          call getatomsymfunctions(i1,i2,iindex,natoms,atomindex,natoms,&
            max_num_atoms,max_num_neighbors_elec,&
            invneighboridx_elec,jcount,listdim,lsta,lstc,lste,symelement_elec,maxnum_funcvalues_elec,&
            cutoff_type,nelem,function_type_elec,&
            lstb,funccutoff_elec,xyzstruct,symfunctione,dsfuncdxyz_dummy,strs_dummy,&
            eta_elec,zeta_elec,lambda_elec,rshift_elec,rmin,&
            .false.,.false.)
          if(.not.lrmin)then
            write(ounit,*)'Error in prediction: lrmin=.false. (atoms too close)'
            stop !'
          endif
!!====================================================================================
!! check for extrapolation atoms n_start to n_end
!! This needs to be done before scaling 
!!====================================================================================
          if(lfinetime)then
            dayextrapolationewald=0
            call abstime(timeextrapolationewaldstart,dayextrapolationewald)
          endif ! lfinetime
          lextrapolation=.false.
          if((symfunctione(i2)-maxvalue_local(ielem,i2)).gt.threshold)then
            write(ounit,'(a,a5,x,2i5,x,a,2f18.8)')&
            '### EXTRAPOLATION WARNING ### ','Ewald',&
            atomindex(i1),i2,'too large ',symfunctione(i2),maxvalue_local(ielem,i2)
            lextrapolation=.true.
          elseif((-symfunctione(i2)+minvalue_local(ielem,i2)).gt.threshold)then
            write(ounit,'(a,a5,x,2i5,x,a,2f18.8)')&
            '### EXTRAPOLATION WARNING ### ','Ewald',&
            atomindex(i1),i2,'too small ',symfunctione(i2),minvalue_local(ielem,i2)
            lextrapolation=.true.
          endif
          if(lfinetime)then
            call abstime(timeextrapolationewaldend,dayextrapolationewald)
            timeextrapolationewald=timeextrapolationewald+timeextrapolationewaldend-timeextrapolationewaldstart
          endif ! lfinetime
!!====================================================================================
!! scale symmetry functions for the electrostatic interaction
!! caution: internally nblock and npoints are set to 1 to avoid _list in zelem, symfunctione and num_atoms
!!====================================================================================
          if(lfinetime)then
            dayscalesymewald=0
            call abstime(timescalesymewaldstart,dayscalesymewald)
          endif ! lfinetime
          if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
            symfunctione(i2)=symfunctione(i2)-avvalue_local(ielem,i2)
          elseif(lscalesym.and..not.lcentersym)then
!! Scale each symmetry function value for the respective element
            symfunctione(i2)=(symfunctione(i2)&
           -minvalue_local(ielem,i2))/ &
           (maxvalue_local(ielem,i2)-minvalue_local(ielem,i2))&
            *(scmax_local-scmin_local) + scmin_local
          elseif(lscalesym.and.lcentersym)then
            symfunctione(i2)=(symfunctione(i2)&
            -avvalue_local(ielem,i2))/ &
            (maxvalue_local(ielem,i2)-minvalue_local(ielem,i2))
          else
          endif
          if(lfinetime)then
            call abstime(timescalesymewaldend,dayscalesymewald)
            timescalesymewald=timescalesymewald+timescalesymewaldend-timescalesymewaldstart
          endif ! lfinetime
        enddo ! i2 loop over all symmetry functions
!!
!!====================================================================================
!! now we have all symmetry functions of atom i1/jcount, now calculate the atom energy 
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!!====================================================================================
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
        call calconenn(1,maxnum_funcvalues_elec,maxnodes_elec,&
          maxnum_layers_elec,num_layers_elec(ielem),&
          maxnum_weights_elec,nodes_elec(0,ielem),&
          symfunctione,weights_elec(1,ielem),nodes_values,nodes_sum,&
          nnoutput,actfunc_elec(1,1,ielem))
        nnatomcharge(jcount)=nnoutput
!!
        if(lsens)then
!!====================================================================================
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
!!====================================================================================
          call getdnodes_values(maxnum_layers_elec,num_layers_elec,maxnodes_elec,&
            nodes_elec,ielem,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_elec)
!!
!!====================================================================================
!! calculate dchargedsfunc for this atom i1  
!! calculate the full derivative of E_i with respect to G_j^i dchargedsfunc
!!====================================================================================
          tempderivative(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_ewald(1) values for each input node
          do i2=1,num_funcvalues_elec(ielem)
            do i4=1,nodes_elec(1,ielem) ! over all nodes in layer 1 ("target layer") 
              icount=(i2-1)*nodes_elec(1,ielem)+i4 ! set pointer in weights array, don't need windexe for first weights
              tempderivative(1,i4,i2)=dnodes_values(1,i4)*weights_elec(icount,ielem)
            enddo ! i4
          enddo ! i2
!! for layers 2 and beyond (if present) 
          if(num_layers_elec(ielem).gt.1)then
            do i2=1,num_funcvalues_elec(ielem)
              do i5=2,num_layers_elec(ielem) ! over all hidden and output layers
                do i3=1,nodes_elec(i5,ielem) ! over all nodes in the target layer
!! we have to sum over the nodes in the previous layer (i4)
                  do i4=1,nodes_elec(i5-1,ielem) ! sum over all nodes in previous layer
                    icount=windex_elec(2*i5-1,ielem)+(i4-1)*nodes_elec(i5,ielem)+i3-1 
                    tempderivative(i5,i3,i2)=tempderivative(i5,i3,i2) + &
                    dnodes_values(i5,i3)*weights_elec(icount,ielem)*tempderivative(i5-1,i4,i2)
                  enddo ! i4
                enddo ! i3
              enddo ! i5
            enddo ! i2
          endif
          do i2=1,num_funcvalues_elec(ielem)
            dchargedsfunc(i2)=tempderivative(num_layers_elec(ielem),1,i2)
          enddo ! i2
!!
!!====================================================================================
!! calculation of the sensitivity
!!====================================================================================
          do i2=1,num_funcvalues_elec(ielem)
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!          sense(ielem,i2)=sense(ielem,i2)+dchargedsfunc(i2)
            sense(ielem,i2)=sense(ielem,i2)+dchargedsfunc(i2)**2.0
!! END CHANGE ANDI
          enddo
        endif ! lsens
!!
        jcount=jcount+1
      enddo ! i1 ! natoms
!!
      return
      end
