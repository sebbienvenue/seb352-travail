!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

      module structures 
!!
      use globaloptions
!!
      implicit none

      integer, dimension(:)      , allocatable :: num_atoms_list 
      integer, dimension(:)      , allocatable :: num_pairs_list
!!
      integer, dimension(:,:)    , allocatable :: zelem_list 
      integer, dimension(:,:,:)  , allocatable :: zelemp_list 
      integer, dimension(:,:,:)  , allocatable :: zelemtrip_list

      real*8, dimension(:,:,:)  , allocatable :: lattice_list 
      real*8, dimension(:,:,:)  , allocatable :: xyzstruct_list  
      real*8, dimension(:)      , allocatable :: totalcharge_list  
      real*8, dimension(:)      , allocatable :: totalenergy_list  
      real*8, dimension(:)      , allocatable :: shortenergy_list  
      real*8, dimension(:)      , allocatable :: elecenergy_list  
      real*8, dimension(:,:,:)  , allocatable :: totalforce_list  
      real*8, dimension(:,:,:)  , allocatable :: totforce_list  
      real*8, dimension(:,:,:)  , allocatable :: elecforce_list  
      real*8, dimension(:,:,:)  , allocatable :: nntbforce_list  
      real*8, dimension(:,:,:)  , allocatable :: shortforce_list  
      real*8, dimension(:,:)    , allocatable :: atomcharge_list  
      real*8, dimension(:,:)    , allocatable :: atomenergy_list  
      real*8, dimension(:,:,:)  , allocatable :: hextoff_list

      logical, dimension(:)     , allocatable :: lperiodic_list  

      character*2, dimension(:,:) , allocatable :: elementsymbol_list  

!!      integer num_atoms_list(nblock)
!!      integer num_pairs_list(nblock)
!!      integer zelem_list(nblock,max_num_atoms)
!!      integer zelemp_list(2,nblock,max_num_pairs)

!!      real*8 lattice_list(3,3,nblock)
!!      real*8 xyzstruct_list(3,max_num_atoms,nblock)
!!      real*8 totalcharge_list(nblock)
!!      real*8 totalenergy_list(nblock)
!!      real*8 shortenergy_list(nblock)
!!      real*8 elecenergy_list(nblock)
!!      real*8 totalforce_list(3,max_num_atoms,nblock)
!!      real*8 elecforce_list(3,max_num_atoms,nblock)
!!      real*8 totforce_list(3,max_num_atoms,nblock)
!!      real*8 atomcharge_list(nblock,max_num_atoms)
!!      real*8 atomenergy_list(nblock,max_num_atoms)

!!      logical lperiodic_list(nblock)

!!      character*2 elementsymbol_list(nblock,max_num_atoms)

      end module structures 

