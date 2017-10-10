!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

      module nnham 
!!
      implicit none
!!

      integer, dimension(:)  , allocatable :: num_layers_s
      integer, dimension(:,:), allocatable :: nodes_s
      integer, dimension(:,:), allocatable :: windex_s  
      integer, dimension(:)  , allocatable :: num_weights_s  
      integer, dimension(:)  , allocatable :: num_funcvalues_s  
      integer maxnodes_s  

      real*8, dimension(:,:)   , allocatable :: weights_s  
      real*8, dimension(:,:,:) , allocatable :: symfunction_s_list
      real*8 scmin_s  
      real*8 scmax_s  

      character*1, dimension(:,:,:), allocatable :: actfunc_s  

      integer, dimension(:)  , allocatable :: num_layers_hextoff
      integer, dimension(:,:), allocatable :: nodes_hextoff
      integer, dimension(:,:), allocatable :: windex_hextoff
      integer, dimension(:)  , allocatable :: num_weights_hextoff
      integer, dimension(:)  , allocatable :: num_funcvalues_hextoff
      integer maxnodes_hextoff

      real*8, dimension(:,:)   , allocatable :: weights_hextoff
      real*8, dimension(:,:) , allocatable :: symfunction_hextoff_list
      real*8 scmin_hextoff
      real*8 scmax_hextoff

      character*1, dimension(:,:,:), allocatable :: actfunc_hextoff

      integer, dimension(:)  , allocatable :: num_layers_hexton
      integer, dimension(:,:), allocatable :: nodes_hexton
      integer, dimension(:,:), allocatable :: windex_hexton
      integer, dimension(:)  , allocatable :: num_weights_hexton
      integer, dimension(:)  , allocatable :: num_funcvalues_hexton
      integer maxnodes_hexton

      real*8, dimension(:,:)   , allocatable :: weights_hexton
      real*8, dimension(:,:,:) , allocatable :: symfunction_hexton_list
      real*8 scmin_hexton
      real*8 scmax_hexton

      character*1, dimension(:,:,:), allocatable :: actfunc_hexton

      integer, dimension(:)  , allocatable :: num_layers_dens
      integer, dimension(:,:), allocatable :: nodes_dens
      integer, dimension(:,:), allocatable :: windex_dens
      integer, dimension(:)  , allocatable :: num_weights_dens
      integer, dimension(:)  , allocatable :: num_funcvalues_dens
      integer maxnodes_dens

      real*8, dimension(:,:)   , allocatable :: weights_dens
      real*8, dimension(:,:,:) , allocatable :: symfunction_dens_list
      real*8 scmin_dens
      real*8 scmax_dens

      character*1, dimension(:,:,:), allocatable :: actfunc_dens

      integer, dimension(:)  , allocatable :: num_layers_ham
      integer, dimension(:,:), allocatable :: nodes_ham
      integer, dimension(:,:), allocatable :: windex_ham
      integer, dimension(:)  , allocatable :: num_weights_ham
      integer, dimension(:)  , allocatable :: num_funcvalues_ham
      integer maxnodes_ham

      real*8, dimension(:,:)   , allocatable :: weights_ham
      real*8, dimension(:,:,:) , allocatable :: symfunction_ham_list
      real*8 scmin_ham
      real*8 scmax_ham

      character*1, dimension(:,:,:), allocatable :: actfunc_ham

!!      integer :: hextoff_training_triplet(3)

      end module nnham

