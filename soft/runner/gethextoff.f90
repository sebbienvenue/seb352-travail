!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!!
      subroutine gethextoff(matrixsize,npoints,&
        symfunction_hextoff_local,nnhextoff_list)
!!
      use fileunits
      use globaloptions
      use basismod
      use nnham
!!
      implicit none
!!
      integer ndim                                                    ! in
      integer npoints                                                 ! in
      integer i1,i2,i3,i4
      integer matrixsize                                              ! in
!!
      real*8 symfunction_hextoff_local(maxnum_funcvalues_hextoff,nblock)  ! in
      real*8 nnhextoff_list(nblock,&
               num_basis(elementindex(hextoff_training_triplet(1))),&
               num_basis(elementindex(hextoff_training_triplet(2)))) ! out
      real*8 nodes_sum_dummy(maxnum_layers_hextoff,maxnodes_hextoff)      ! internal 
      real*8 nodes_values_dummy(maxnum_layers_hextoff,maxnodes_hextoff)   ! internal 
      real*8 nnoutput(matrixsize)
!!
!!
!!
!!
!!
      do i1=1,npoints
!!
!! calculate nodes_values and nnoutput for atom i1
        call calconenn(matrixsize,maxnum_funcvalues_hextoff,maxnodes_hextoff,&
          maxnum_layers_hextoff,num_layers_hextoff,&
          maxnum_weights_hextoff,nodes_hextoff,&
          symfunction_hextoff_local(:,i1),weights_hextoff(:,1),&
          nodes_values_dummy,nodes_sum_dummy,&
          nnoutput,actfunc_hextoff(:,:,1))
!!
        i4 = 0
        do i2=1,num_basis(elementindex(hextoff_training_triplet(1)))
          do i3=1,num_basis(elementindex(hextoff_training_triplet(2)))
            i4 = i4 + 1
            nnhextoff_list(i1,i2,i3)=nnoutput(i4)
          enddo
        enddo
!!
      enddo ! i1
!!
      return
      end
