!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine
!! CAUTION: This routine is not yet prepared for pair NN, some nelems still have to be replaced (not all!)
!! get fixed weights for short range and electrostatic NN

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!! - fitting_batch.f90
!!
      subroutine getfixedweights(iswitch,ndim,&
          maxnum_layers_local,num_layers_local,&
          nodes_local,maxnum_weights_local,num_weights_local,windex_local,&
          num_weights_localfree,num_weights_localfixed,wconstraint_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2                     ! internal
      integer ndim                      ! in
      integer iswitch                   ! in
      integer maxnum_layers_local        ! in
      integer num_layers_local(ndim)    ! in
      integer maxnum_weights_local          ! in
      integer num_weights_local(ndim)          ! in
      integer num_weights_localfree(ndim)      ! out 
      integer num_weights_localfixed(ndim)     ! out 
      integer nodes_local(0:maxnum_layers_local,ndim) ! in
      integer layer1                    ! internal
      integer layer2                    ! internal
      integer node1                     ! internal
      integer node2                     ! internal
      integer icount                    ! internal
      integer wconstraint_local(maxnum_weights_local,ndim) !
      integer windex_local(2*maxnum_layers_local,ndim) ! in
      integer isum(ndim)               ! internal
      integer ztemp                     ! internal

      character*40 dummy                ! internal
      character*40 keyword              ! internal
      character*40 reference            ! internal
      character*10 cswitch              ! internal
      character*2 elementtemp           ! internal
      character*2 elementtemp1          ! internal
      character*2 elementtemp2          ! internal
      character*2 elementtemp3          ! internal

!! initializations
      wconstraint_local(:,:)=0
      if(iswitch.eq.0)then
        reference='weight_constraint'
      elseif(iswitch.eq.1)then
        reference='weighte_constraint'
      elseif(iswitch.eq.2)then
        if(nntb_flag(3))then
          reference='weighthextoff_constraint'
        endif
      endif
!!
!! read constraints from file
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 10   continue
      read(nnunit,*,END=20)keyword     
      layer1=0
      layer2=0
      node1=0
      node2=0 
      if(keyword.eq.reference)then
        backspace(nnunit)
        if(nntb_flag(3)) then
           write(*,*) 'read fixed hextoff'
           STOP !! CMH this all needs to be finished but for now assume we never use constraints and fixed weights.
           read(nnunit,*,ERR=99)keyword,elementtemp1,elementtemp2,elementtemp3,dummy
!! check elements match the hextoff_training_triplet
           if((elementtemp1.eq.element(elementindex(hextoff_training_triplet(1)))).and.&
              (elementtemp2.eq.element(elementindex(hextoff_training_triplet(2)))).and.&
              (elementtemp3.eq.element(elementindex(hextoff_training_triplet(3))))) then
             goto 40
           endif
           write(ounit,*) 'Error: unknown triplet in constraint', elementtemp1,elementtemp2,elementtemp3
        else
          read(nnunit,*,ERR=99)keyword,elementtemp,dummy
!! chekc if element is known
          do i1=1,nelem
            if(elementtemp.eq.element(i1))then
              goto 40
            endif
          enddo
          write(ounit,*)'Error: unknown element in constraint ',elementtemp
          stop
        endif
 40     continue
!!        write(ounit,*)'getweightconstraints: ',keyword,elementtemp,dummy
!!
!! set all weights adjacent to a certain node 
        STOP
        if(dummy.eq.'node')then
          backspace(nnunit)
          read(nnunit,*,END=30)keyword,elementtemp,dummy,layer1,node1,cswitch
          call nuccharge(elementtemp,ztemp)
          if(layer1.gt.num_layers_local(elementindex(ztemp)))then
            write(ounit,*)'ERROR in weight constraints (node): layer1 too large ',layer1,num_layers_local(elementindex(ztemp))
            stop
          endif
          if(node1.gt.nodes_local(layer1,elementindex(ztemp)))then
            write(ounit,*)'ERROR in weight constraints (node): node1 too large ',node1,nodes_local(layer1,elementindex(ztemp))
            stop
          endif
!! set connecting weights to node1
          if(layer1.gt.0)then
            do i1=1,nodes_local(layer1-1,elementindex(ztemp))
              icount=windex_local(2*layer1-1,elementindex(ztemp))
              icount=icount+node1-1
              icount=icount+(i1-1)*nodes_local(layer1,elementindex(ztemp))
              if(cswitch.eq.'fixed')then
                wconstraint_local(icount,elementindex(ztemp))=1
              elseif(cswitch.eq.'free')then
                wconstraint_local(icount,elementindex(ztemp))=0
              else
                write(ounit,*)'ERROR: weights can either be fixed or free'
                stop
              endif
            enddo ! i1'
          endif !layer1.gt.0
!! set connecting weights from node1
          if(layer1.lt.num_layers_local(elementindex(ztemp)))then
            do i1=1,nodes_local(layer1+1,elementindex(ztemp))
              icount=windex_local(2*(layer1+1)-1,elementindex(ztemp))
              icount=icount+(node1-1)*nodes_local(layer1+1,elementindex(ztemp))
              icount=icount+i1-1
              if(cswitch.eq.'fixed')then
                wconstraint_local(icount,elementindex(ztemp))=1
              elseif(cswitch.eq.'free')then
                wconstraint_local(icount,elementindex(ztemp))=0
              else
                write(ounit,*)'ERROR: weights can either be fixed or free'
                stop !'
              endif
            enddo ! i1
          endif ! layer.lt.num_layers_local
!! set also bias weight
          icount=windex_local(2*layer1,elementindex(ztemp))+node1-1
          if(cswitch.eq.'fixed')then
            wconstraint_local(icount,elementindex(ztemp))=1
          elseif(cswitch.eq.'free')then
            wconstraint_local(icount,elementindex(ztemp))=0
          else
            write(ounit,*)'ERROR: weights can either be fixed or free'
            stop
          endif
!!
!! set all connecting weights between 2 layers
        elseif(dummy.eq.'interlayer')then
          backspace(nnunit)
          read(nnunit,*,END=30)keyword,elementtemp,dummy,layer1,layer2,cswitch
          call nuccharge(elementtemp,ztemp)
!!          write(ounit,*)dummy,layer1,layer2,cswitch
          if((layer1+1).ne.layer2)then
            write(ounit,*)'ERROR in weight constraints (interlayer)'
            stop
          endif
          if(layer2.gt.num_layers_local(elementindex(ztemp)))then
            write(ounit,*)'ERROR in weight constraints (interlayer): layer2 too large'
            stop !'
          endif
          do i1=windex_local(2*layer2-1,elementindex(ztemp)),windex_local(2*layer2,elementindex(ztemp))-1
            icount=i1
            if(cswitch.eq.'fixed')then
              wconstraint_local(icount,elementindex(ztemp))=1
            elseif(cswitch.eq.'free')then
              wconstraint_local(icount,elementindex(ztemp))=0
            else
              write(ounit,*)'ERROR: weights can either be fixed or free'
              stop
            endif
          enddo ! i1
!!
!! set a specific connecting weight 
        elseif(dummy.eq.'weight')then
          backspace(nnunit)
          read(nnunit,*,END=30)keyword,elementtemp,dummy,layer1,node1,layer2,node2,cswitch
          call nuccharge(elementtemp,ztemp)
          if((layer1+1).ne.layer2)then
            write(ounit,*)'ERROR in weight constraints (weight)'
            stop
          endif
          if(node1.gt.nodes_local(layer1,elementindex(ztemp)))then
            write(ounit,*)'ERROR in weight constraints (weight): node1 too large'
            stop
          endif
          if(node2.gt.nodes_local(layer2,elementindex(ztemp)))then
            write(ounit,*)'ERROR in weight constraints (weight): node2 too large'
            stop
          endif
          if(layer2.gt.num_layers_local(elementindex(ztemp)))then
            write(ounit,*)'ERROR in weight constraints (weight): layer2 too large'
            stop !'
          endif
          icount=windex_local(2*(layer1+1)-1,elementindex(ztemp))
          icount=icount+(node1-1)*nodes_local(layer2,elementindex(ztemp))
          icount=icount+node2-1
          if(cswitch.eq.'fixed')then
            wconstraint_local(icount,elementindex(ztemp))=1
          elseif(cswitch.eq.'free')then
            wconstraint_local(icount,elementindex(ztemp))=0
          else
            write(ounit,*)'ERROR: weights can either be fixed or free'
            stop
          endif
!!
!! set all weights 
        elseif(dummy.eq.'all')then
          backspace(nnunit)
          read(nnunit,*,END=30)keyword,elementtemp,dummy,cswitch
          call nuccharge(elementtemp,ztemp)
          if(cswitch.eq.'fixed')then
            wconstraint_local(:,elementindex(ztemp))=1
          elseif(cswitch.eq.'free')then
            wconstraint_local(:,elementindex(ztemp))=0
          else
            write(ounit,*)'ERROR: weights can either be fixed or free'
            stop
          endif
!!
!! set a bias weight 
        elseif(dummy.eq.'bias')then
          backspace(nnunit)
          read(nnunit,*,END=30)keyword,elementtemp,dummy,layer1,node1,cswitch
          call nuccharge(elementtemp,ztemp)
          icount=windex_local(2*layer1,elementindex(ztemp))+node1-1
          if(cswitch.eq.'fixed')then
            wconstraint_local(icount,elementindex(ztemp))=1
          elseif(cswitch.eq.'free')then
            wconstraint_local(icount,elementindex(ztemp))=0
          else
            write(ounit,*)'ERROR: weights can either be fixed or free'
            stop
          endif
!!
        else
          write(ounit,*)'ERROR: unknown weight constraint ',elementtemp,dummy
          stop
        endif
      endif !(keyword.eq.reference)then
      goto 10
!!
 20   continue
      close(nnunit)
!!
!! check if there are still free weights for optimization
      isum(:)=0
      do i1=1,nelem
        do i2=1,num_weights_local(i1)
          isum(i1)=isum(i1)+wconstraint_local(i2,i1)
        enddo
        if(isum(i1).eq.num_weights_local(i1))then
          write(ounit,*)'ERROR: all weights are constrained ',iswitch
          write(ounit,*)'for element ',i1
          stop
        endif
        num_weights_localfree(i1)=num_weights_local(i1)-isum(i1)
        num_weights_localfixed(i1)=isum(i1)
      enddo ! i1
!!
!!
      return
!!
!! incomplete constraints file:
 30   continue
      write(ounit,*)'ERROR in weight constraints file'
      stop
!!
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop

      end
