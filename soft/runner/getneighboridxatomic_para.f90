!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! Goal: determine max_num_neighbors as dimension for large array dsfuncdxyz
!! calculate lsta,lstb,lstc,lste

!! called by:
!!
      subroutine getneighboridxatomic_para(n_start,natoms,listdim,&
        max_num_atoms,max_num_neighbors,&
        lsta,lstc,neighboridx_local,invneighboridx_local)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2                                          ! internal
      integer n_start                                        ! in
      integer natoms                                         ! in
      integer iatom                                          ! internal
      integer listdim                                        ! in
      integer lsta(2,max_num_atoms)                          ! in, numbers of neighbors
      integer lstc(listdim)                                  ! in, identification of atom
      integer neighboridx_local(natoms,0:max_num_neighbors)  ! out, neighboridx(*,0) is number of central atom itself
      integer invneighboridx_local(natoms,max_num_atoms)     ! out, yields number of neighbor (invneighboridx(X)<=max_num_neighbors)
      integer icount                                         ! internal
      integer max_num_atoms                                  ! in 
      integer max_num_neighbors                              ! in
!!
!! initialize as -1 (not allowed value) to cause crash in case of bugs
      neighboridx_local(:,:)   =-1
      invneighboridx_local(:,:)=-1 
!!
!! There is no unique correspondence between neighbor number and atom number
!!
      do i1=1,natoms
        icount=0
!! central atom itself:        
        iatom=n_start+i1-1
        invneighboridx_local(i1,iatom)=0
        neighboridx_local(i1,icount)=iatom
        do i2=lsta(1,n_start+i1-1),lsta(2,n_start+i1-1) ! loop over neighbors
          icount=icount+1 ! counter for neighbors
!! get atom number for each neighbor icount of each atom i1
!! value range is between 0 and max_num_neighbors
          invneighboridx_local(i1,lstc(i2))=icount 
!! get neighbor index for each atom i1
!! value range is between 1 and max_num_atoms
          neighboridx_local(i1,icount)=lstc(i2)
!!
        enddo
!!
      enddo ! i1
!!
      return
      end
