!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine deallocatestructures()
     
      use globaloptions
      use structures 
      use nnflags
      implicit none

      deallocate (num_atoms_list)
      deallocate (num_pairs_list)
      deallocate (zelem_list)
      deallocate (zelemp_list)
      deallocate (lattice_list)
      deallocate (xyzstruct_list)
      deallocate (totalcharge_list)
      deallocate (totalenergy_list)
      deallocate (shortenergy_list)
      deallocate (elecenergy_list)
      deallocate (totalforce_list)
      deallocate (elecforce_list)
      deallocate (nntbforce_list)
      deallocate (shortforce_list)
      deallocate (totforce_list)
      deallocate (atomcharge_list)
      deallocate (atomenergy_list)
      deallocate (lperiodic_list)
      deallocate (elementsymbol_list)
      if(lnntb) deallocate(hextoff_list)

      return
      end
