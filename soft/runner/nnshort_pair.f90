!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

      module nnshort_pair 
!!
      implicit none
!!
      integer, dimension(:)  , allocatable :: num_layers_short_pair
      integer, dimension(:,:), allocatable :: nodes_short_pair
      integer, dimension(:,:), allocatable :: windex_short_pair
      integer, dimension(:)  , allocatable :: num_weights_short_pair
      integer, dimension(:)  , allocatable :: num_funcvalues_short_pair
      integer maxnodes_short_pair

      real*8, dimension(:,:)   , allocatable :: weights_short_pair
      real*8, dimension(:,:,:) , allocatable :: symfunction_short_pair_list
      real*8 scmin_short_pair
      real*8 scmax_short_pair

      character*1, dimension(:,:,:), allocatable :: actfunc_short_pair

      end module nnshort_pair 

