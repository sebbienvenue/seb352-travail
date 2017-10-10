!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine distribute_globaloptions()
!!
      use mpi_mod
      use nnflags
      use globaloptions
!!
      implicit none
!!
      call mpi_bcast(nran,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(paramode,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nelem,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(cutoff_type,1,mpi_integer,0,mpi_comm_world,mpierror)
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(maxnum_weights_short_atomic,1,mpi_integer,0,mpi_comm_world,mpierror)
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(maxnum_weights_short_pair,1,mpi_integer,0,mpi_comm_world,mpierror)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(maxnum_weights_elec,1,mpi_integer,0,mpi_comm_world,mpierror)
      endif
      if(lnntb)then
       if(maxnum_weights_ham.ge.1)then
         call mpi_bcast(maxnum_weights_s,1,mpi_integer,0,mpi_comm_world,mpierror)
       endif
       if(nntb_flag(1))then
         call mpi_bcast(maxnum_weights_s,1,mpi_integer,0,mpi_comm_world,mpierror)
       endif
       if(nntb_flag(2))then
         call mpi_bcast(maxnum_weights_hexton,1,mpi_integer,0,mpi_comm_world,mpierror)
       endif
       if(nntb_flag(3))then
         call mpi_bcast(maxnum_weights_hextoff,1,mpi_integer,0,mpi_comm_world,mpierror)
       endif
       if(nntb_flag(4))then
         call mpi_bcast(maxnum_weights_dens,1,mpi_integer,0,mpi_comm_world,mpierror)
       endif

      endif

      call mpi_bcast(elempair,2*npairs,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(elemtriplet,3*ntriplets,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nucelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nblock,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(listdim,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ewaldkmax,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(max_num_neighbors_atomic_input,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(elementindex,102,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(pairindex,102*102,mpi_integer,0,mpi_comm_world,mpierror)

      call mpi_bcast(rscreen_cut,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(rscreen_onset,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(rmin,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(ewaldalpha,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(ewaldcutoff,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomrefenergies,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(fixedcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(saturation_threshold,1,mpi_real8,0,mpi_comm_world,mpierror)

      call mpi_bcast(luseatomcharges,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseatomenergies,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lscalesym,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lcentersym,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lremoveatomenergies,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfinetime,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lreadunformatted,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lwriteunformatted,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lcheckf,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lompmkl,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseforces,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(ldebug,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(ldostress,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lscreen,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lsilent,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lnormnodes,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lsens,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lpearson_correlation,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lenvironmentanalysis,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lenforcemaxnumneighborsatomic,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lmd,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(ldetect_saturation,1,mpi_logical,0,mpi_comm_world,mpierror)

!! check if pstring is distributed correctly!
      call mpi_bcast(pstring,1,mpi_character,0,mpi_comm_world,mpierror)
!! check if element is distributed correctly!
      call mpi_bcast(element,1,mpi_character,0,mpi_comm_world,mpierror)

!! NNTB stuff
      call mpi_bcast(tripletindex,102*102*102,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(triplettag,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(tripletmode12,3,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(hextoff_training_triplet,3,mpi_integer,0,mpi_comm_world,mpierror)
!!
      return
      end
