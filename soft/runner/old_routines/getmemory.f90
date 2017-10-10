!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine getmemory(ounit,max_num_atoms,nelem,&
        num_layersshort,num_layersewald,&
        num_funcvalues,num_funcvaluese,nblock,&
        num_weightsshort,num_weightsewald,&
        maxnodes_short,maxnodes_ewald,ldebug)
!!
      implicit none
!!
      integer ounit
      integer max_num_atoms
      integer nelem
      integer num_layersshort
      integer num_layersewald
      integer num_funcvalues
      integer num_funcvaluese
      integer nblock
      integer num_weightsshort
      integer num_weightsewald
      integer maxnodes_short
      integer maxnodes_ewald
      integer memory      
!!
      logical ldebug
!!
      memory=0
!!
!! xyzstruct_list, xyzforce_list
!!      memory=memory+3*nblock*max_num_atoms*2
!!
!! weight arrays
!!
!!
!!
!!
!! final conversion to a useful unit
!!      memory=memory*8
!!      
      write(ounit,*)'estimated memory: ',memory,' bytes'
      write(ounit,*)'-------------------------------------------------------------'
!!
      return
!!
      end
