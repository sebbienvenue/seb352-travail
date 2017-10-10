!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - predictionshortatomic.f90
!!
      subroutine getshortatomic(n_start,natoms,atomindex,&
        max_num_neighbors_short_atomic,num_atoms,&
        invneighboridx_short_atomic,num_neighbors_short_atomic,&
        neighboridx_short_atomic,zelem,&
        lsta,lstc,lste,lstb,xyzstruct,&
        sens,nnshortforce,nnstress_short,minvalue_local,maxvalue_local,&
        avvalue_local,scmin_local,scmax_local,nnatomenergy,&
        lextrapolation,lperiodic)
!!
      use fileunits
      use globaloptions
      use nnshort_atomic
      use symfunctions
      use timings
      use predictionoptions
      use wrapper, only : extrapolated_atoms
!!
      implicit none
!!
      integer n_start                                                   ! in
      integer num_atoms                                                 ! in
      integer max_num_neighbors_short_atomic                            ! in
      integer ielem                                                     ! internal
      integer iindex                                                    ! internal
      integer zelem(max_num_atoms)                                      ! in
      integer i1,i2,i3,i4,i5                                            ! internal
      integer icount                                                    ! internal
      integer jcount                                                    ! internal
      integer natoms                                                    ! in
      integer atomindex(natoms)                                         ! in
      integer neighboridx_short_atomic(natoms,0:max_num_neighbors_short_atomic) ! in
      integer num_neighbors_short_atomic(num_atoms)                     ! in
      integer invneighboridx_short_atomic(natoms,max_num_atoms)         ! in
      integer lsta(2,max_num_atoms)                                     ! in
      integer lstc(listdim)                                             ! in, identification of atom
      integer lste(listdim)                                             ! in, nuclear charge of atom
!!
      real*8 xyzstruct(3,max_num_atoms)                                 ! in
      real*8 lstb(listdim,4)                                            ! in, xyz and r_ij
      real*8 minvalue_local(nelem,maxnum_funcvalues_short_atomic)       ! in 
      real*8 maxvalue_local(nelem,maxnum_funcvalues_short_atomic)       ! in 
      real*8 avvalue_local(nelem,maxnum_funcvalues_short_atomic)        ! in 
      real*8 scmin_local                                                ! in
      real*8 scmax_local                                                ! in
      real*8 strs(3,3,maxnum_funcvalues_short_atomic)                   ! internal???                   
      real*8 dsfuncdxyz_temp(maxnum_funcvalues_short_atomic,0:max_num_neighbors_short_atomic,3)   ! internal
      real*8 dsfuncdxyz_local(0:max_num_neighbors_short_atomic,3)       ! internal
      real*8 nnshortforce(3,max_num_atoms)                              ! out 
      real*8 deshortdsfunc(maxnum_funcvalues_short_atomic)              ! internal 
      real*8 symfunction(maxnum_funcvalues_short_atomic)                ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                                   ! internal 
      real*8 nodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic)            ! internal
      real*8 nodes_sum(maxnum_layers_short_atomic,maxnodes_short_atomic)               ! internal
      real*8 dnodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic)           ! internal
      real*8 tempderivative(maxnum_layers_short_atomic,maxnodes_short_atomic,maxnum_funcvalues_short_atomic) ! internal
      real*8 nnatomenergy(max_num_atoms)                                ! out 
      real*8 nnstress_short(3,3)                                        ! out 
      real*8 threshold                                                  ! internal
      real*8 sens(nelem,maxnum_funcvalues_short_atomic)                 ! out 
!!  
      logical lrmin 
      logical lextrapolation
      logical lperiodic
!!
!! initialization
      lrmin           =.true.
      threshold       = 0.0001d0
!!
!! TODO
!! rm iindex


!!====================================================================================
!! loop over all atoms 
!!====================================================================================
      jcount=n_start
      do i1=1,natoms
!!====================================================================================
!! initializations for this atom 
!!====================================================================================
        deshortdsfunc(:)= 0.0d0
        symfunction(:)  = 0.0d0
        ielem=elementindex(zelem(atomindex(i1)))
        iindex=elementindex(zelem(jcount))
!!====================================================================================
!! get the symmetry functions for atom jcount 
!!====================================================================================
        do i2=1,num_funcvalues_short_atomic(ielem) ! over all symmetry functions
          call getatomsymfunctions(i1,i2,iindex,natoms,atomindex,natoms,&
            max_num_atoms,max_num_neighbors_short_atomic,&
            invneighboridx_short_atomic,jcount,listdim,lsta,lstc,lste,&
            symelement_short_atomic,maxnum_funcvalues_short_atomic,&
            cutoff_type,nelem,function_type_short_atomic,&
            lstb,funccutoff_short_atomic,xyzstruct,symfunction,dsfuncdxyz_local,strs,&
            eta_short_atomic,zeta_short_atomic,lambda_short_atomic,rshift_short_atomic,rmin,&
            ldoforces,ldostress)
          dsfuncdxyz_temp(i2,:,:)=dsfuncdxyz_local(:,:)
          if(.not.lrmin)then
            write(ounit,*)'Error in prediction: lrmin=.false. (atoms too close)'
            stop !'
          endif
!!====================================================================================
!! check for extrapolation atoms n_start to n_end
!! This needs to be done before scaling 
!!====================================================================================
          if(lfinetime)then
            dayextrapolationshort=0
            call abstime(timeextrapolationshortstart,dayextrapolationshort)
          endif ! lfinetime
          lextrapolation=.false.
          if((symfunction(i2)-maxvalue_local(ielem,i2)).gt.threshold)then
            extrapolated_atoms(atomindex(i1)) = .True.
            lextrapolation=.true.
          elseif((-symfunction(i2)+minvalue_local(ielem,i2)).gt.threshold)then
            extrapolated_atoms(atomindex(i1)) = .True.
            lextrapolation=.true.
          endif
          if(lfinetime)then
            call abstime(timeextrapolationshortend,dayextrapolationshort)
            timeextrapolationshort=timeextrapolationshort+timeextrapolationshortend-timeextrapolationshortstart
          endif ! lfinetime
!!====================================================================================
!! scale symmetry functions for the short-range interaction
!! caution: internally nblock and npoints are set to 1 to avoid _list in zelem, symfunction and num_atoms
!!====================================================================================
          if(lfinetime)then
            dayscalesymshort=0
            call abstime(timescalesymshortstart,dayscalesymshort)
          endif ! lfinetime
          if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
            symfunction(i2)=symfunction(i2)-avvalue_local(ielem,i2)
          elseif(lscalesym.and..not.lcentersym)then
!! Scale each symmetry function value for the respective element
            symfunction(i2)=(symfunction(i2)&
           -minvalue_local(ielem,i2))/ &
           (maxvalue_local(ielem,i2)-minvalue_local(ielem,i2))&
            *(scmax_local-scmin_local) + scmin_local
          elseif(lscalesym.and.lcentersym)then
            symfunction(i2)=(symfunction(i2)&
            -avvalue_local(ielem,i2))/ &
            (maxvalue_local(ielem,i2)-minvalue_local(ielem,i2))
          else
          endif
          if(lfinetime)then
            call abstime(timescalesymshortend,dayscalesymshort)
            timescalesymshort=timescalesymshort+timescalesymshortend-timescalesymshortstart
          endif ! lfinetime
!!====================================================================================
!! scale dsfuncdxyz for forces if requested
!!====================================================================================
          if(ldoforces.and.lscalesym)then
            if(lfinetime)then
              dayscaledsfuncshort=0
              call abstime(timescaledsfuncshortstart,dayscaledsfuncshort)
            endif ! lfinetime
            do i3=0,num_neighbors_short_atomic(jcount)  ! over all atoms in structure
              do i4=1,3 ! x,y,z
                dsfuncdxyz_temp(i2,i3,i4)=dsfuncdxyz_temp(i2,i3,i4)/&
                (maxvalue_local(ielem,i2)-minvalue_local(ielem,i2))&
                *(scmax_local-scmin_local)
              enddo ! i4
            enddo ! i3
            if(lfinetime)then
              call abstime(timescaledsfuncshortend,dayscaledsfuncshort)
              timescaledsfuncshort=timescaledsfuncshort+timescaledsfuncshortend-timescaledsfuncshortstart
            endif ! lfinetime
          endif ! ldoforces
!!====================================================================================
!! scale stress components, CHECK IF THIS IS RIGHT!!!
!!====================================================================================
          if(ldostress)then
            strs(:,:,i2)=strs(:,:,i2)/&
            (maxvalue_local(ielem,i2)-minvalue_local(ielem,i2))*(scmax_local-scmin_local)
          endif
        enddo ! i2 loop over all symmetry functions
!!
!!====================================================================================
!! now we have all symmetry functions of atom i1/jcount, now calculate the atom energy 
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!!====================================================================================
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
        call calconenn(1,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic(ielem),&
          maxnum_weights_short_atomic,nodes_short_atomic(0,ielem),&
          symfunction,weights_short_atomic(1,ielem),nodes_values,nodes_sum,&
          nnoutput,actfunc_short_atomic(1,1,ielem))
        nnatomenergy(jcount)=nnoutput
!!
!!====================================================================================
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
!!====================================================================================
        call getdnodes_values(maxnum_layers_short_atomic,num_layers_short_atomic,maxnodes_short_atomic,&
          nodes_short_atomic,ielem,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_short_atomic)
!!
!!====================================================================================
!! calculate deshortdsfunc for this atom i1  
!! calculate the full derivative of E_i with respect to G_j^i deshortdsfunc
!!====================================================================================
        tempderivative(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_short_atomic(1) values for each input node
        do i2=1,num_funcvalues_short_atomic(ielem)
          do i4=1,nodes_short_atomic(1,ielem) ! over all nodes in layer 1 ("target layer") 
            icount=(i2-1)*nodes_short_atomic(1,ielem)+i4 ! set pointer in weights array, don't need windex_short_atomic for first weights
            tempderivative(1,i4,i2)=dnodes_values(1,i4)*weights_short_atomic(icount,ielem)
          enddo ! i4
        enddo ! i2
!! for layers 2 and beyond (if present) 
        if(num_layers_short_atomic(ielem).gt.1)then
          do i2=1,num_funcvalues_short_atomic(ielem)
            do i5=2,num_layers_short_atomic(ielem) ! over all hidden and output layers
              do i3=1,nodes_short_atomic(i5,ielem) ! over all nodes in the target layer
!! we have to sum over the nodes in the previous layer (i4)
                do i4=1,nodes_short_atomic(i5-1,ielem) ! sum over all nodes in previous layer
                  icount=windex_short_atomic(2*i5-1,ielem)+(i4-1)*nodes_short_atomic(i5,ielem)+i3-1 
                  tempderivative(i5,i3,i2)=tempderivative(i5,i3,i2) + &
                  dnodes_values(i5,i3)*weights_short_atomic(icount,ielem)*tempderivative(i5-1,i4,i2)
                enddo ! i4
              enddo ! i3
            enddo ! i5
          enddo ! i2
        endif
        do i2=1,num_funcvalues_short_atomic(ielem)
          deshortdsfunc(i2)=tempderivative(num_layers_short_atomic(ielem),1,i2)
        enddo ! i2
!!
!!====================================================================================
!! calculate force
!!====================================================================================
        do i2=1,num_funcvalues_short_atomic(ielem)
          do i3=0,num_neighbors_short_atomic(jcount)  ! over all atoms in structure
            do i4=1,3 ! x,y,z
              nnshortforce(i4,neighboridx_short_atomic(i1,i3))=nnshortforce(i4,neighboridx_short_atomic(i1,i3)) & 
                -deshortdsfunc(i2)*dsfuncdxyz_temp(i2,i3,i4)
            enddo ! i4
          enddo ! i3
        enddo ! i2
!!
!!====================================================================================
!! calculation of the sensitivity
!!====================================================================================
        if(lsens)then
          do i2=1,num_funcvalues_short_atomic(ielem)
!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!          sens(ielem,i2)=sens(ielem,i2)+deshortdsfunc(i2)
            sens(ielem,i2)=sens(ielem,i2)+deshortdsfunc(i2)**2.0
!! END CHANGE ANDI
          enddo
        endif ! lsens
!!
!!====================================================================================
!! calculation of short range stress 
!!====================================================================================
        if(lfinetime)then
          daysshort=0
          call abstime(timesshortstart,daysshort)
        endif ! lfinetime
        if(ldostress.and.lperiodic)then
          do i2=1,num_funcvalues_short_atomic(ielem)
            do i4=1,3
              do i3=1,3
                nnstress_short(i3,i4)=nnstress_short(i3,i4)&
                -strs(i3,i4,i2)*deshortdsfunc(i2)
              enddo ! i3
            enddo ! i4
          enddo ! i2
        endif ! ldostress
        if(lfinetime)then
          call abstime(timesshortend,daysshort)
          timesshort=timesshort+timesshortend-timesshortstart
        endif ! lfinetime
!!
        jcount=jcount+1
      enddo ! i1 ! natoms
!!
      return
      end
