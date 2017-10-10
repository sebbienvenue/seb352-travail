!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine getenergies(ounit,nelem,nblock,npoints,&
        max_num_atoms,num_funcvalues,num_funcvaluese,elementindex,&
        nnindex,nnindexe,function_type,function_typee,&
        num_atoms_list,windex,windexe,maxnodes_short,maxnodes_ewald,&
        num_layersshort,num_layersewald,nodes_short,nodes_ewald,&
        num_weightsshort,num_weightsewald,ewaldkmax,zelem_list,&
        num_functions,num_functionse,listdim,nucelem,ncharges,&     
        funccutoff,funccutoffe,eta,etae,lambda,lambdae,&
        zeta,zetae,rshift,rshifte,maxcutoff,maxcutoffe,&
        minvalue,maxvalue,avvalue,minvaluee,maxvaluee,avvaluee,&
        symfunction_list,symfunctione_list,lattice_list,&
        nnxyzforce_list,nnewaldforce_list,xyzstruct_list,&
        weights_short,weights_ewald,nneshort_list,nnewald_list,&
        rmse_short,rmse_charge,rmse_totalcharge,rmse_ewald,&
        shortenergy_list,ewaldenergy_list,&
        atomcharge_list,nnatomcharge_list,&
        totalcharge_list,nnchargesum_list,&
        actfunc_short,actfunc_ewald,&
        ewaldalpha,ewaldcutoff,&
        lperiodic_list,lscalesym,lcentersym,luseforces,&
        lshort,lewald,ldebug)
!!
      implicit none
!!
      integer ounit                                      ! in
      integer nelem                                      ! in
      integer nblock                                     ! in
      integer npoints                                    ! in
      integer max_num_atoms                              ! in
      integer num_funcvalues                             ! in 
      integer num_funcvaluese                            ! in 
      integer num_atoms_list(nblock)                     ! in
      integer zelem_list(nblock,max_num_atoms)           ! in
      integer elementindex(102)                          ! in
      integer num_layersshort                            ! in
      integer num_layersewald                            ! in
      integer maxnodes_short                             ! in
      integer maxnodes_ewald                             ! in
      integer windex(2*num_layersshort)                  ! in
      integer windexe(2*num_layersewald)                 ! in
      integer num_weightsshort                           ! in
      integer num_weightsewald                           ! in
      integer nodes_short(0:num_layersshort)             ! in
      integer nodes_ewald(0:num_layersewald)             ! in
      integer num_functions                              ! in
      integer num_functionse                             ! in
      integer listdim                                    ! in
      integer nucelem(nelem)                             ! in                 
      integer ncharges                                   ! out
      integer ewaldkmax                                  ! in
      integer day                                        ! internal
      integer nnindex(num_functions,nelem,nelem)         ! in
      integer nnindexe(num_functionse,nelem,nelem)       ! in
      integer function_type(num_functions)               ! in
      integer function_typee(num_functionse)             ! in
!!
      real*8 minvalue(nelem,num_funcvalues)              ! in
      real*8 minvaluee(nelem,num_funcvaluese)            ! in
      real*8 maxvalue(nelem,num_funcvalues)              ! in
      real*8 maxvaluee(nelem,num_funcvaluese)            ! in
      real*8 avvalue(nelem,num_funcvalues)               ! in
      real*8 avvaluee(nelem,num_funcvaluese)             ! in
      real*8 symfunction_list(num_funcvalues,max_num_atoms,nblock)     ! in/out
      real*8 symfunctione_list(num_funcvaluese,max_num_atoms,nblock)   ! in/out
      real*8 weights_short(num_weightsshort,nelem)       ! in
      real*8 weights_ewald(num_weightsewald,nelem)       ! in
      real*8 nneshort_list(nblock)                       ! out
      real*8 nnewald_list(nblock)                        ! out
      real*8 nnxyzforce_list(nblock,max_num_atoms,3)     ! out
      real*8 nnewaldforce_list(nblock,max_num_atoms,3)   ! out
      real*8 lattice_list(3,3,nblock)                    ! in
      real*8 xyzstruct_list(3,max_num_atoms,nblock)      ! in
      real*8 rmse_short                                  ! out
      real*8 rmse_charge                                 ! out
      real*8 rmse_totalcharge                            ! out
      real*8 rmse_ewald                                  ! out
      real*8 shortenergy_list(nblock)                    ! in 
      real*8 ewaldenergy_list(nblock)                    ! in
      real*8 nnatomcharge_list(nblock,max_num_atoms)     ! out
      real*8 totalcharge_list(nblock)                    ! in
      real*8 atomcharge_list(nblock,max_num_atoms)       ! in
      real*8 nnchargesum_list(nblock)                    ! out
      real*8 ewaldalpha                                  ! in
      real*8 ewaldcutoff                                 ! in
      real*8 timestart                                   ! internal
      real*8 timeend                                     ! internal
!! symmetry function parameters
      real*8 funccutoff(num_functions)                   ! in
      real*8 funccutoffe(num_functionse)                 ! in
      real*8 eta(num_functions)                          ! in
      real*8 etae(num_functionse)                        ! in
      real*8 lambda(num_functions)                       ! in
      real*8 lambdae(num_functionse)                     ! in
      real*8 zeta(num_functions)                         ! in
      real*8 zetae(num_functionse)                       ! in
      real*8 rshift(num_functions)                       ! in
      real*8 rshifte(num_functionse)                     ! in
      real*8 maxcutoff                                   ! in
      real*8 maxcutoffe                                  ! in
!!
      character*1 actfunc_short(num_layersshort)         ! in
      character*1 actfunc_ewald(num_layersewald)         ! in
!!
      logical lshort                                     ! in
      logical lewald                                     ! in
      logical luseforces                                 ! in
      logical lperiodic_list(nblock)                     ! in
      logical lscalesym                                  ! in
      logical lcentersym                                 ! in
      logical ldebug                                     ! in
      logical lfinetime                                  ! internal only here for testing
!!
!! fine timing just for debugging here
      lfinetime=.false.
!!
      if(lshort)then
!! scale symmetry functions for the short-range interaction
        call scalesym(ounit,nelem,nblock,npoints,&
          max_num_atoms,num_funcvalues,num_atoms_list,&
          zelem_list,elementindex,symfunction_list,&
          minvalue,maxvalue,avvalue,&
          lscalesym,lcentersym,ldebug)
!!
!! predict the short range NN output for npoint data sets
        call geteshort(ounit,nelem,nblock,npoints,num_weightsshort,&
          maxnodes_short,&
          num_layersshort,windex,&
          num_funcvalues,max_num_atoms,zelem_list,&
          num_atoms_list,nodes_short,elementindex,&
          symfunction_list,weights_short,&
          nneshort_list,&
          actfunc_short,&
          ldebug)
!!
!! predict the short range NN forces for the test points here
        if(luseforces)then
          call getallshortforces(ounit,nelem,nblock,npoints,&
            num_funcvalues,max_num_atoms,num_weightsshort,&
            num_layersshort,maxnodes_short,num_functions,listdim,&
            num_atoms_list,zelem_list,nucelem,nnindex,function_type,&
            elementindex,windex,nodes_short,&
            symfunction_list,weights_short,nnxyzforce_list,&
            lattice_list,xyzstruct_list,minvalue,maxvalue,&
            funccutoff,eta,rshift,zeta,lambda,maxcutoff,&
            actfunc_short,lperiodic_list,lscalesym,ldebug)
        endif ! luseforces
!!
!! calculate rmse_short: in/out rmse_short
        call calcrmse_energy(nblock,npoints,rmse_short,&
          shortenergy_list,nneshort_list,ldebug)
!!
      endif ! lshort
!!
!! predict the charges and the electrostatic energy
!!-------------------------------------------------
!!
      if(lewald)then
!! scale the symmetry functions for the charge prediction
        call scalesym(ounit,nelem,nblock,npoints,&
          max_num_atoms,num_funcvaluese,num_atoms_list,&
          zelem_list,elementindex,symfunctione_list,&
          minvaluee,maxvaluee,avvaluee,&
          lscalesym,lcentersym,ldebug)
!!
!! calculate the charges on the atoms
        call getcharges(ounit,nelem,nblock,npoints,num_weightsewald,&
          maxnodes_ewald,&
          num_layersewald,windexe,&
          num_funcvaluese,max_num_atoms,zelem_list,&
          num_atoms_list,nodes_ewald,elementindex,&
          symfunctione_list,weights_ewald,&
          nnatomcharge_list,&
          actfunc_ewald,&
          ldebug)
!!
!! calculate the RMSE for the atomic charges and the total charge
        call calcrmse_charge(nblock,npoints,ncharges,&
          max_num_atoms,num_atoms_list,rmse_charge,&
          totalcharge_list,rmse_totalcharge,&
          atomcharge_list,nnatomcharge_list,nnchargesum_list)
!!
!! calculate the electrostatic energy for npoints points
        call getallelectrostatic(ounit,nblock,npoints,&
          nelem,num_funcvaluese,num_layersewald,num_functionse,&
          maxnodes_ewald,elementindex,windexe,nucelem,&
          nodes_ewald,num_weightsewald,nnindexe,function_typee,&
          max_num_atoms,num_atoms_list,ewaldkmax,listdim,&
          zelem_list,lattice_list,minvaluee,maxvaluee,&
          nnatomcharge_list,xyzstruct_list,nnewald_list,&
          weights_ewald,symfunctione_list,&
          ewaldalpha,ewaldcutoff,nnewaldforce_list,&
          funccutoffe,etae,lambdae,zetae,rshifte,maxcutoffe,&
          actfunc_ewald,lperiodic_list,&
          luseforces,lscalesym,ldebug)
!!
!! calculate the RMSE for the electrostatic energy: in/out rmse_ewald 
        call calcrmse_energy(nblock,npoints,rmse_ewald,&
          ewaldenergy_list,nnewald_list,ldebug)
!!
      endif  ! lewald
!!
      return
      end
