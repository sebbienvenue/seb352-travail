!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

      module nnshort_atomic 
!!
      implicit none
!!
      integer, dimension(:)  , allocatable :: num_layers_short_atomic
      integer, dimension(:,:), allocatable :: nodes_short_atomic
      integer, dimension(:,:), allocatable :: windex_short_atomic
      integer, dimension(:)  , allocatable :: num_weights_short_atomic
      integer, dimension(:)  , allocatable :: num_funcvalues_short_atomic
      integer maxnodes_short_atomic

      real*8, dimension(:,:)   , allocatable :: weights_short_atomic
      real*8, dimension(:,:,:) , allocatable :: symfunction_short_atomic_list
      real*8 scmin_short_atomic
      real*8 scmax_short_atomic

      character*1, dimension(:,:,:), allocatable :: actfunc_short_atomic

      end module nnshort_atomic 

