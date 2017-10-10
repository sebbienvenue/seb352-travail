!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine allocatestructures()
     
      use globaloptions
      use structures
      use nnflags
      use basismod

      implicit none

      allocate (num_atoms_list(nblock))
      allocate (num_pairs_list(nblock))
      allocate (zelem_list(nblock,max_num_atoms))
      allocate (zelemp_list(2,nblock,max_num_pairs))
      allocate (zelemtrip_list(3,nblock,max_num_triplets))
      allocate (lattice_list(3,3,nblock))
      allocate (xyzstruct_list(3,max_num_atoms,nblock))
      allocate (totalcharge_list(nblock))
      allocate (totalenergy_list(nblock))
      allocate (shortenergy_list(nblock))
      allocate (elecenergy_list(nblock))
      allocate (totalforce_list(3,max_num_atoms,nblock))
      allocate (totforce_list(3,max_num_atoms,nblock))
      allocate (elecforce_list(3,max_num_atoms,nblock))
      allocate (nntbforce_list(3,max_num_atoms,nblock))
      allocate (shortforce_list(3,max_num_atoms,nblock))
      allocate (atomcharge_list(nblock,max_num_atoms))
      allocate (atomenergy_list(nblock,max_num_atoms))
      allocate (lperiodic_list(nblock))
      allocate (elementsymbol_list(nblock,max_num_atoms))
      if(lnntb) allocate(hextoff_list(nblock,maxnum_basis,maxnum_basis))

      return
      end
