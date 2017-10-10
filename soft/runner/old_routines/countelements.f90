!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine countelements(max_num_atoms,zelem,nelem)
      implicit none
!!
      integer max_num_atoms
      integer nelem          ! number of elements
      integer nelemtemp 
      integer zelem(max_num_atoms)
      integer i,j
!!
!!
      nelemtemp = 1
!!
      do i=2,max_num_atoms
        do j=1,i-1
          if(zelem(i).eq.zelem(j)) goto 20 
        enddo ! j
        nelemtemp=nelemtemp+1
 20     continue
      enddo ! i
!!
      nelem=max(nelem,nelemtemp)
!!      write(*,*)'nelem = ',nelem
!!
      return
      end
