!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purposes:
!! - read zelem_list for nntype 2 case since in this case it is not available from function.data
!! - read structures

!! called by: 
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine getstructures_mixed(npoints,&
           ntrain,block_counter,pointindex,num_atoms_all)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!
      integer num_atoms_all(ntrain)                 ! in
      integer npoints                               ! in
      integer zelem_local(max_num_atoms)            ! internal 
      integer i1,i2,i3                              ! internal
      integer ntrain                                ! in
      integer pointindex(ntrain)                    ! in
      integer block_counter                         ! in
      integer jcount                                ! internal
      integer idummy                                ! internal
!!
      real*8 force_local(3,max_num_atoms)           ! internal dummy, no output!
      real*8 lattice_local(3,3)                     ! internal
      real*8 atomenergy_local(max_num_atoms)        ! internal
      real*8 atomcharge_local(max_num_atoms)        ! internal
      real*8 xyzstruct_local(3,max_num_atoms)       ! internal
!!
      logical lperiodic_local                       ! internal
!!
!!      write(ounit,*)'getstructures_mixed npoints ',npoints,block_counter
!!
      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      rewind(trainstructunit)
!!
!! here we read only a subset of strutures
      do i1=1,ntrain
!! read structure data
        lperiodic_local=.false.
        read(trainstructunit,*)idummy,lperiodic_local
        if(lperiodic_local)then
          do i2=1,3
            read(trainstructunit,*)(lattice_local(i2,i3),i3=1,3)
          enddo
        endif
        do i2=1,num_atoms_all(i1)
          read(trainstructunit,*)zelem_local(i2),&
            (xyzstruct_local(i3,i2),i3=1,3),&
            atomcharge_local(i2),&
            atomenergy_local(i2),&
            (force_local(i3,i2),i3=1,3)
        enddo
!! check if we need this structure now
        jcount=block_counter
        do i2=1,npoints
          if(i1.eq.pointindex(jcount))then
!!            write(ounit,'(a,i5,a,i5,a,i5)')'JBgetstructures structure ',i1,' is number ',jcount,' and goes to field ',i2
!! copy the data of this structure to the list arrays
            do i3=1,num_atoms_all(i1)
              xyzstruct_list(:,i3,i2)=xyzstruct_local(:,i3)
              atomcharge_list(i2,i3) =atomcharge_local(i3)
              atomenergy_list(i2,i3) =atomenergy_local(i3)
              lattice_list(:,:,i2)   =lattice_local(:,:)
              lperiodic_list(i2)     =lperiodic_local
!! zelem_local is not output here, comes from readfunctions_mixed
!! forces_local are no output here, comes from readforces_mixed
            enddo ! i3
          endif
          jcount=jcount+1
        enddo ! i2
      enddo ! i1
 10   continue
!!
      close(trainstructunit)
!!
      return
      end
