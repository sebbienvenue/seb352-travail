!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine getwconstraintidx(iswitch,ndim,windex_local,&
        maxnum_layers_local,num_layers_local,maxnum_weights_local,&
        num_weights_local,num_weightsfree_local,num_weightsfixed_local,&
        nodes_local,wconstraintidx_local)

      use mpi_mod
      use fileunits
      use nnflags
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer i1,i2,i0,ldim,mdim
      integer icount
      integer ndim                                            ! in
      integer iswitch                                         ! in
      integer wconstraint(maxnum_weights_local,ndim)          ! internal 
      integer num_weights_local(ndim)                         ! in
      integer num_weightsfree_local(ndim)                     ! out
      integer num_weightsfixed_local(ndim)                    ! out
      integer maxnum_weights_local                            ! in
      integer maxnum_layers_local                             ! in
      integer num_layers_local(ndim)                          ! in
      integer windex_local(2*maxnum_layers_local,ndim)        ! in
      integer nodes_local(0:maxnum_layers_local,ndim)         ! in 
      integer wconstraintidx_local(maxnum_weights_local,ndim) ! out 
!!
      wconstraint(:,:)         =0
      num_weightsfree_local(:) =num_weights_local(:)
      num_weightsfixed_local(:)=0
!!
!!    CMH REMEMBER!!!!!!!!
!!    if we are training for hextoff we only train one triplet at a time
!!    the lable of the triplet i.e. tripletindex is the important thing
!!
      if((lnntb).and.(nntb_flag(3)))then
        i0 = tripletindex(hextoff_training_triplet(1),hextoff_training_triplet(2),hextoff_training_triplet(3))
      endif
      if(mpirank.eq.0)then
!!
        if(.not.lfixweights) goto 10
!!
!! FIXME: it is not clear if this works for nn_type_short 2:
        call getfixedweights(iswitch,ndim,&
          maxnum_layers_local,num_layers_local,&
          nodes_local,maxnum_weights_local,num_weights_local,windex_local,&
          num_weightsfree_local,num_weightsfixed_local,wconstraint)
!!
        
        if((iswitch.eq.0).and.(nn_type_short.eq.1))then
          do i1=1,ndim
            write(ounit,*)element(i1),'  Total number of short range weights       ',num_weights_local(i1)
            write(ounit,*)element(i1),'  Number of optimized short range weights   ',num_weightsfree_local(i1)
            write(ounit,*)element(i1),'  Number of constrained short range weights ',num_weightsfixed_local(i1)
          enddo
        elseif((iswitch.eq.0).and.(nn_type_short.eq.2))then
          do i1=1,ndim
            write(ounit,*)element(pairindex(1,i1)),element(pairindex(2,i1)),&
              '  Total number of short range weights       ',num_weights_local(i1)
            write(ounit,*)element(pairindex(1,i1)),element(pairindex(2,i1)),&
              '  Number of optimized short range weights   ',num_weightsfree_local(i1)
            write(ounit,*)element(pairindex(1,i1)),element(pairindex(2,i1)),&
              '  Number of constrained short range weights ',num_weightsfixed_local(i1)
          enddo
        elseif(iswitch.eq.1)then
          do i1=1,ndim
            write(ounit,*)element(i1),'  Total number of charge weights            ',num_weights_local(i1)
            write(ounit,*)element(i1),'  Number of optimized charge weights        ',num_weightsfree_local(i1)
            write(ounit,*)element(i1),'  Number of constrained charge weights      ',num_weightsfixed_local(i1)
          enddo
        elseif(iswitch.eq.2)then
          if(nntb_flag(3))then
            i1 = i0
              write(ounit,*) hextoff_training_triplet,'  Total number of hextoff weights       ',num_weights_local(i1)
              write(ounit,*) hextoff_training_triplet,'  Number of optimized hextoff weights   ',num_weightsfree_local(i1)
              write(ounit,*) hextoff_training_triplet,'  Number of constrained hextoff weights ',num_weightsfixed_local(i1)
            
          endif
        endif
        
 10     continue
      endif ! mpirank.eq.0
!!
      
      call mpi_bcast(num_weightsfree_local,ndim,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsfixed_local,ndim,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!
!! get wconstraintidx_local
      wconstraintidx_local = 0
      mdim=1
      ldim=ndim
      if((lnntb).and.(nntb_flag(3)))then
        mdim=1
        ldim=1
      endif
      do i1=mdim,ldim
        icount=0
        do i2=1,num_weights_local(i1)
          if(wconstraint(i2,i1).eq.0)then
            icount=icount+1
            wconstraintidx_local(icount,i1)=i2
          endif
        enddo
      enddo
      
!! 
      return
      end
