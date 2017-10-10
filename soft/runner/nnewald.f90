!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

      module nnewald 
!!
      implicit none
!!
      integer, dimension(:)  , allocatable :: num_layers_elec
      integer, dimension(:,:), allocatable :: nodes_elec
      integer, dimension(:,:), allocatable :: windex_elec
      integer, dimension(:)  , allocatable :: num_weights_elec
      integer, dimension(:)  , allocatable :: num_funcvalues_elec
      integer maxnodes_elec

      real*8, dimension(:,:)   , allocatable :: weights_elec
      real*8, dimension(:,:,:) , allocatable :: symfunction_elec_list
      real*8 scmin_elec
      real*8 scmax_elec

      character*1, dimension(:,:,:), allocatable :: actfunc_elec

      end module nnewald 

