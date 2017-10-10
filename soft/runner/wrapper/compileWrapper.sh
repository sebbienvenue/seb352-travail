#!/bin/bash
source settings.sh

# MPI has to be selected manually
mpi_mod=(../mpi_module.f90)
# mpi_mod=(../mpi_dummy.f90
#          ../mpi_dummy_routines.f90)

# all files which create a .mod
data_mods=(../basismod.f90
           ../fileunits.f90
           ../fittingoptions.f90
           ../globaloptions.f90
           ../inputnncounters.f90
           ../mode1options.f90
           ../nnewald.f90
           ../nnflags.f90
           ../nnham.f90
           ../nnshort_atomic.f90
           ../nnshort_pair.f90
           ../predictionoptions.f90
           ../symfunctions.f90
           ../timings.f90)

# subroutines directly taken from RuNNer (parent directory)
subroutines=(../abstime.f90
             ../addatoms.f90
             ../atomsymfunction1.f90
             ../atomsymfunction2.f90
             ../atomsymfunction3.f90
             ../atomsymfunction4.f90
             ../atomsymfunction5.f90
             ../atomsymfunction6.f90
             ../atomsymfunction8.f90
             ../atomsymfunction9.f90
             ../calconenn.f90
             ../checkelement.f90
             ../distribute_fittingoptions.f90
             ../distribute_globaloptions.f90
             ../distribute_predictionoptions.f90
             ../getatomsymfunctions.f90
             ../getdimensions.f90
             ../getdnodes_values.f90
             ../getlistdim.f90
             ../getneighboridxatomic_para.f90
             ../getneighborsatomic_para.f90
             ../getnumpairs.f90
             ../getvolume.f90
             ../initializecounters.f90
             ../initmode3.f90
             ../inputnndefaults.f90
             ../mpifitdistribution.f90
             ../neighbor_para_short.f90
             ../nuccharge.f90
             ../paircount.f90
             ../predictionshortatomic.f90
             ../preparemd.f90
             ../readbasis.f90
             ../readelementlayersatomic.f90
             ../readelementlayerspair.f90
             ../readsymfunctionatomic.f90
             ../readsymfunctionelementatomic.f90
             ../readsymfunctionelementpair.f90
             ../readsymfunctionglobalatomic.f90
             ../readsymfunctionglobalpair.f90
             ../readsymfunctionpair.f90
             ../readweights.f90
             ../setglobalactivation.f90
             ../sortelements.f90
             ../sortpairsymfunctions.f90
             ../sortsymfunctions.f90
             ../structurecount.f90
             ../translate.f90
             ../zerotime.f90)

# subroutines which differ from RuNNer (thus in this directory)
local_subroutines=(atomsymfunction3Andi.f90
                   checkfiles.f90
                   checkinputnn.f90
                   getshortatomic.f90
                   initialization.f90
                   initnn.f90
                   readatomenergies.f90
                   readkeywords.f90
                   readinput.f90
                   readscale.f90)

f2py3 --f90exec=${mpi_compiler} \
      --fcompiler=${compiler_type} \
      --f90flags="${fflags}" \
      --opt="${opt}" \
      -c wrapper.pyf \
      ${mpi_mod[@]} \
      ${data_mods[@]} \
      ${subroutines[@]} \
      ${local_subroutines[@]} \
      wrapper.f90
