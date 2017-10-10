      module globaloptions 

      implicit none

!! cutoff_type = functional form of cutoff function (cutoff_type=1: f_c= 0.5(cos(pi*R_ij/R_c)+1)  , cutoff_type=1: f_c=(tanh(1.d0-rij/funccutoff_local(i2,iindex)))**3 )
      integer cutoff_type

!! nran = integer switch for selecting the type of random number generator
      integer nran

!! nblock = number of structure in memory at the same time (modes 1 + 2) or number of atoms in memory at the same time (mode 3)
      integer nblock

!! max_num_atoms: number of atoms in the largest structure of all structures in input.data
      integer max_num_atoms 
      integer max_num_pairs
      integer max_num_triplets ! max number of unique triplets (is this even useful.... it's here consistency with max_num_pairs!!! CMH

!! elementindex: yields number of element in system when nuclear charge of atom is provided. Example: for pure ZnO we have elementindex(6)=1 and elementindex(30)=2, all others are 0. 
      integer elementindex(102)

!! pairindex
      integer pairindex(102,102)
!! tripletindex
      integer tripletindex(102,102,102)  !! FIXME> should be allocated
!!      integer, dimension(:,:,:), allocatable :: tripletindex
      integer, dimension(:,:), allocatable :: storetriplets
      integer tripletmode12(3) !! for modes 1 and 2 training of Hextoff we need to know which triplet we are using
      integer triplettag       !! of all the possible triplets this flag is the tag for it
      integer hextoff_training_triplet(3)
!! maxnum_weightsshort:
      integer maxnum_weights_short_atomic
      integer maxnum_weights_elec
      integer maxnum_weights_short_pair
      integer maxnum_weights_s  
      integer maxnum_weights_hexton
      integer maxnum_weights_hextoff
      integer maxnum_weights_dens
      integer maxnum_weights_ham



!! maxnum_layers_short_atomic = number of short range hidden layers + 1 (atomic NN)
!! determined in main/initialization/getdimensions
      integer maxnum_layers_short_atomic

!! maxnum_layersewald = number of electrostatic NN hidden layers + 1 (atomic and pair NN)
!! determined in main/initialization/getdimensions
      integer maxnum_layers_elec

!! maxnum_layerspair = number of short range hidden layers + 1 (pair NN)
!! determined in main/initialization/getdimensions
      integer maxnum_layers_short_pair

!! maxnum_layersham = number of hidden layers + 1 (ham NN)
!! determined in main/initialization/getdimensions
      integer maxnum_layers_s
      integer maxnum_layers_hexton
      integer maxnum_layers_hextoff
      integer maxnum_layers_dens
      integer maxnum_layers_ham

      integer maxnum_funcvalues_short_atomic
      integer maxnum_funcvalues_elec
      integer maxnum_funcvalues_short_pair
      integer maxnum_funcvalues_s
      integer maxnum_funcvalues_hexton
      integer maxnum_funcvalues_hextoff
      integer maxnum_funcvalues_dens
      integer maxnum_funcvalues_ham


      integer paramode
      integer listdim
      integer ewaldkmax

!! determined in getdimensions
      integer npairs
      integer ntriplets
!! enforced value for max_num_neighbors
      integer max_num_neighbors_atomic_input

!! determined in getdimensions
      integer nelem ! number of elements in input.nn

!! determined in main/initialization/structurecount
      integer totnum_structures ! total number of structures in input.data file

      integer, dimension(:)  , allocatable :: nucelem
      integer, dimension(:,:), allocatable :: elempair
      integer, dimension(:,:), allocatable :: elemtriplet

      real*8 cutoffvdw
      real*8 ewaldalpha
      real*8 ewaldcutoff
      real*8 rscreen_cut
      real*8 rscreen_onset
      real*8, dimension(:) , allocatable :: fixedcharge
      real*8 rmin
      real*8, dimension(:)     , allocatable :: atomrefenergies

      integer, dimension(:,:)  , allocatable :: vdw_param 

      logical lscalesym
      logical lcentersym
      logical lremoveatomenergies
      logical lnormnodes
      logical lscreen
      logical lreadunformatted
      logical lwriteunformatted
      logical lfinetime
      logical lfinetimeepoch
      logical lsilent
      logical ldebug
      logical lompmkl
      logical lcheckf
      logical luseatomcharges
      logical luseatomenergies
      logical luseforces
      logical lenforcemaxnumneighborsatomic
      logical lmd

!! calculate correlation of symmetry functions
      logical lpearson_correlation

!! calculate the sensitivity
      logical lsens

      logical ldostress

      logical ldetect_saturation
      real*8 saturation_threshold

!! give statistics of atomic environments
      logical lenvironmentanalysis

      character*2, dimension(:), allocatable :: element
      character*20 pstring

!! CMH Hamiltonian global options

      logical nntb_flag(0:4) ! this is a flag container for turning on and off the hamiltonian symmfunction determination for
                           ! overlap, hexton, hextoff, dens

      end module globaloptions 

