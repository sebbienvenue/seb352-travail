!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: fit the atomic short range part and if requested atomic charges as second output nodes

!! called by:
!!
      subroutine readcharges(max_num_atoms,num_atoms,nnatomcharge)

      use fileunits
      use mpi_mod

      implicit none

      integer icount                                 ! internal
      integer max_num_atoms                          ! in
      integer num_atoms                              ! in
      integer i1                                     ! internal
      
      real*8 rdummy                                  ! internal
      real*8 nnatomcharge(max_num_atoms)             ! out
!!
      logical lexist                                 ! internal

!!======================================================================
!! file reading for parallel processes only for mpirank 0
!!======================================================================
      if(mpirank.eq.0)then
!!======================================================================
!! check if file charges.in exists
!!======================================================================
        inquire(file='charges.in',exist=lexist)
        if(.not.lexist)then
          write(ounit,*)'ERROR: file charges.in not found'
          stop
        endif

!!======================================================================
!! check if file charges.in contains the right number of atoms
!!======================================================================
        icount=0
        open(chargeunit,file='charges.in',status='old')
        rewind(chargeunit)
 10     continue
        read(chargeunit,*,END=11)rdummy
        icount=icount+1
        goto 10
 11     continue
        close(chargeunit)
        if(icount.ne.num_atoms)then
          write(ounit,*)'ERROR: charges.in does not contain the right number of atoms ',icount,num_atoms
          stop !'
        endif

!!======================================================================
!!  read charges
!!======================================================================
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Reading charges from file charges.in:'
        nnatomcharge(:)=0.0d0
        open(chargeunit,file='charges.in',status='old')
        rewind(chargeunit)
        do i1=1,num_atoms
          read(chargeunit,*)nnatomcharge(i1)
          write(ounit,'(i5,f14.8)')i1,nnatomcharge(i1)
        enddo
        close(chargeunit)
!!
      endif ! mpirank.eq.0 '
!!
!!======================================================================
!! distribute charges to all processes in case of parallel run 
!!======================================================================
      if(mpisize.gt.1)then
        call mpi_bcast(nnatomcharge,&
          max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      return
      end
