!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - geterror.f90
!! - precondition.f90
!! - fitting_batch.f90
!! - fittingpair.f90 
!! - geterrorpair.f90
!! - preconditionpair.f90
!!
      subroutine readfunctions_hextoff_mixed(iswitch,npoints,ndim,ntrain,&
         pointindex,block_counter,max_num,maxnum_funcvalues_local,&
         num_funcvalues_local,&
         symfunction_list_local)
!!               call readfunctions_hextoff_mixed(0,npoints,1,ntrain,&
!!            max_num_atoms,maxnum_funcvalues_hextoff,num_funcvalues_hextoff,&
!!            symfunction_hextoff_list)

!!
      use fileunits
      use globaloptions
      use structures
      use basismod
!!
      implicit none
!!     
      integer funit                                                 ! internal
      integer ndim                                                  ! in
      integer ntrain                                                ! in
      integer iswitch                                               ! in
      integer npoints                                               ! in
      integer max_num                                               ! in
      integer pointindex(ntrain)                                    ! in
      integer icount, jcount                                        ! internal      
      integer block_counter
      integer maxnum_funcvalues_local                               ! in
      integer num_funcvalues_local(ndim)                            ! in
      integer i1,i2,i3,i4                                           ! internal
      real*8  dummy,dummy1,dummy2,dummy3,dummy4                     ! internal
!!
      integer matrixsize
      real*8, dimension(:), allocatable :: hextoffin                ! internal
      real*8 symfunction_list_local(maxnum_funcvalues_local,nblock)     ! out
!!
      matrixsize  = num_basis(elementindex(hextoff_training_triplet(1)))*&
                         num_basis(elementindex(hextoff_training_triplet(2)))
      allocate(hextoffin(matrixsize))
      rewind(symhextoffunit)
      funit = symhextoffunit
      icount = block_counter
      jcount = 1
      do i1=1,ntrain
        if(i1.eq.pointindex(icount)) then
          num_atoms_list(i1) = 3
          read(funit,*)(symfunction_list_local(i3,i1),i3=1,num_funcvalues_local(1))
          read(funit,*)(hextoffin(i2),i2=1,matrixsize)
          i4 = 0
          do i2=1,num_basis(elementindex(hextoff_training_triplet(1)))
            do i3=1,num_basis(elementindex(hextoff_training_triplet(2)))
              i4 = i4 + 1
              hextoff_list(i1,i2,i3) = hextoffin(i4)
            enddo
          enddo
          if(icount.gt.ntrain) goto 10
          if(jcount.gt.npoints) goto 10
        else
          read(funit,*) dummy
          read(funit,*) dummy
        endif

      enddo ! i1
10    continue
      close(funit)
!!
      deallocate(hextoffin)
      return
      end
