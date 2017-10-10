!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine readinput(ielem,iseed,lelement)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use mode1options
      use predictionoptions
      use symfunctions
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use nnham
      use inputnncounters
      use basismod
!!    
      implicit none
!!
      integer icount                                    ! internal
      integer jcount                                    ! internal
      integer nodes_short_atomic_temp(0:maxnum_layers_short_atomic) ! internal
      integer nodes_elec_temp(0:maxnum_layers_elec)                 ! internal
      integer nodes_short_pair_temp(0:maxnum_layers_short_pair)     ! internal
      integer nodes_ham_temp(0:maxnum_layers_ham)                   ! internal
      integer iseed                                     ! out, seed for weight initialization
      integer ielem                                     ! in (read also before from input.nn) 
      integer wcount                                    ! internal 
      integer i,j
      integer i1,i2,i3                                  ! internal
      integer ztemp                                     ! internal
      integer ztemp1                                    ! internal
      integer ztemp2                                    ! internal
      integer itemp                                     ! internal
      integer layer                                     ! internal
      integer node                                      ! internal
      integer sym_short_atomic_count(nelem)             ! internal
      integer sym_elec_count(nelem)                     ! internal
      integer sym_short_pair_count(npairs)              ! internal
      integer sym_ham_count(npairs)                     ! internal
!!
      real*8 kalmanlambda_local                         ! internal
      real*8 kalmanlambdae_local                        ! internal
      real*8 chargetemp                                 ! internal
!!
      character*40 dummy                                ! internal
      character*40 keyword                              ! internal
      character*1 actfunc_short_atomic_dummy(maxnum_layers_short_atomic) ! internal
      character*1 actfunc_elec_dummy(maxnum_layers_elec)                 ! internal
      character*1 actfunc_short_pair_dummy(maxnum_layers_short_pair)     ! internal
      character*1 actfunc_ham_dummy(maxnum_layers_ham)                    ! internal
      character*1 actfunc                               ! internal
      character*2 elementtemp                           ! internal
      character*2 elementtemp1                          ! internal
      character*2 elementtemp2                          ! internal
!!
      logical lelement(102)                             ! internal
!!
!! initializations
!! all inputs must be initialized here, except for the ones read before in getdimensions.f90
!! initialization of counters
      call initializecounters()

!! initializations of input quantities
      if(lshort.and.(nn_type_short.eq.1))then
        nodes_short_atomic_temp(:)   =0
        actfunc_short_atomic_dummy(:)=' '
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        nodes_short_pair_temp(:)     =0
        actfunc_short_pair_dummy(:)  =' '
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        nodes_elec_temp(:)           =0
        actfunc_elec_dummy(:)        =' '
      endif
      if(lnntb)then
        nodes_ham_temp(:)            =0
        actfunc_ham_dummy(:)         =' '
      endif
      kalmanlambda_local =0.98000d0
      kalmanlambdae_local=0.98000d0
      iseed=200
!! read all other defaults
      call inputnndefaults()
!!
!! initializations of derived quantities
      if(lshort.and.(nn_type_short.eq.1))then
        windex_short_atomic(:,:)    =0
        num_weights_short_atomic(:) =0 
        maxnum_weights_short_atomic =0
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        windex_short_pair(:,:)      =0
        num_weights_short_pair(:)   =0 
        maxnum_weights_short_pair   =0
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        windex_elec(:,:)            =0
        num_weights_elec(:)         =0 
        maxnum_weights_elec         =0
      endif
      if(lnntb)then
        windex_ham(:,:)             =0
        num_weights_ham(:)          =0 
        maxnum_weights_ham          =0
      endif
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! start regular reading of input.nn '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      call readkeywords(iseed,&
        nodes_short_atomic_temp,nodes_elec_temp,nodes_short_pair_temp,nodes_ham_temp,&
        kalmanlambda_local,kalmanlambdae_local)
!!
!! TODO: move this to a better place:
!! set output nodes
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          nodes_short_atomic(maxnum_layers_short_atomic,i1)=1 ! initialize as one output node
!! add second output node for atomic charge
          if(lelec.and.(nn_type_elec.eq.2))then
            nodes_short_atomic(maxnum_layers_short_atomic,i1)&
              =nodes_short_atomic(maxnum_layers_short_atomic,i1)+1
          endif
        enddo
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        do i1=1,npairs
          nodes_short_pair(maxnum_layers_short_pair,i1)=1
        enddo
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          nodes_elec(maxnum_layers_elec,i1)=1
        enddo
      endif
!! FIXME: we need to determine the proper number of output nodes here
      if(lnntb)then
        do i1=1,npairs
          nodes_ham(maxnum_layers_ham,i1)=1
        enddo
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read global activation functions 
!! this must be done after reading the global NN architecture
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 50   continue
      read(nnunit,*,END=51) keyword 
!!
      if((keyword.eq.'global_activation_short').and.lshort.and.(nn_type_short.eq.1))then
        if(count_global_activation_short_atomic.gt.0)then
          write(ounit,*)'Error: global_activation_short specified twice'
          stop
        endif
        call setglobalactivation(nelem,count_global_activation_short_atomic,&
          maxnum_layers_short_atomic,maxnodes_short_atomic,nodes_short_atomic,actfunc_short_atomic,&
          actfunc_short_atomic_dummy,keyword)
!!
      elseif((keyword.eq.'global_activation_electrostatic').and.lelec.and.(nn_type_elec.eq.1))then
        if(count_global_activation_elec.gt.0)then
          write(ounit,*)'Error: global_activation_electrostatic specified twice'
          stop !'
        endif
        call setglobalactivation(nelem,count_global_activation_elec,&
          maxnum_layers_elec,maxnodes_elec,nodes_elec,actfunc_elec,&
          actfunc_elec_dummy,keyword)
!!
      elseif((keyword.eq.'global_activation_pair').and.lshort.and.(nn_type_short.eq.2))then
        if(count_global_activation_short_pair.gt.0)then
          write(ounit,*)'Error: global_activation_pair specified twice'
          stop
        endif
        call setglobalactivation(npairs,count_global_activation_short_pair,&
          maxnum_layers_short_pair,maxnodes_short_pair,nodes_short_pair,actfunc_short_pair,&
          actfunc_short_pair_dummy,keyword)
!!
      elseif((keyword.eq.'global_activation_ham').and.lnntb)then
        if(count_global_activation_ham.gt.0)then
          write(ounit,*)'Error: global_activation_ham specified twice'
          stop
        endif
        call setglobalactivation(npairs,count_global_activation_ham,&
          maxnum_layers_ham,maxnodes_ham,nodes_ham,actfunc_ham,&
          actfunc_ham_dummy,keyword)
!!
      endif
      goto 50
 51   continue
      close(nnunit)
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get nuclear charges and sort arrays element and nucelem according to nuclear charge
!! this must be done before reading any element-specific input 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      do i=1,nelem
        call nuccharge(element(i),nucelem(i))
      enddo
      call sortelements()
!! we have to check here if the elements specified in input.nn are the same as the ones
!! occuring in the data set. Every element in input.data must be in input.nn
!! WARNING: The following check will not work for mode 2, because in mode 2 the array
!! lelement is not determined and always false
      do i=1,102
        if(lelement(i)) then ! element with nuclear charge i is present in input.data
          do j=1,nelem
            if(i.eq.nucelem(j)) then
              goto 11
            endif
          enddo
          write(ounit,*)'Error: element ',i,' in input.data is not in input.nn'
          stop
        endif !'
 11     continue
      enddo
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read fixed charges if requested 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(lelec.and.(nn_type_elec.eq.3))then
        open(nnunit,file='input.nn',form='formatted',status='old')
        rewind(nnunit)
 52     continue
        read(nnunit,*,END=53) keyword 
        if(keyword.eq.'fixed_charge')then
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,elementtemp,chargetemp
          call nuccharge(elementtemp,ztemp)
          fixedcharge(elementindex(ztemp))=chargetemp
        endif
        goto 52
 53     continue
        close(nnunit)
!! final check if all elements have been specified
        do i1=1,nelem
          if(fixedcharge(i1).gt.10.0d0)then
            write(ounit,*)'ERROR: No fixed charge specified for element ',element(i1)
            stop !'
          endif
        enddo
      endif 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get the unique element pairs for nn_type_short=2 and nn_type_nntb 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      icount=0
      do i=1,nelem
        do j=i,nelem
          icount=icount+1
          elempair(icount,1)=nucelem(i)
          elempair(icount,2)=nucelem(j)
        enddo
      enddo
!!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read element-specific numbers of layers 
!! this must be done after reading the global NN architecture and activation functions 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! read element-specific numbers of layers, at this point we already need to know all chemical species 
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 40   continue
      read(nnunit,*,END=42) keyword 
!!
      if(keyword.eq.'element_hidden_layers_short')then
        call readelementlayersatomic(maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic,nodes_short_atomic,actfunc_short_atomic)
!!
      elseif(keyword.eq.'element_hidden_layers_electrostatic')then
        call readelementlayersatomic(maxnodes_elec,&
          maxnum_layers_elec,num_layers_elec,nodes_elec,actfunc_elec)
!!
      elseif(keyword.eq.'element_hidden_layers_pair')then
        call readelementlayerspair(maxnodes_short_pair,&
          maxnum_layers_short_pair,num_layers_short_pair,nodes_short_pair,actfunc_short_pair)
!!
      elseif(keyword.eq.'element_hidden_layers_ham')then
        call readelementlayerspair(maxnodes_ham,&
          maxnum_layers_ham,num_layers_ham,nodes_ham,actfunc_ham)
!!
      endif
      goto 40
 42   continue
      close(nnunit)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read element-specific numbers of nodes 
!! this must be done after reading the global NN architecture and activation functions 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 43   continue
      read(nnunit,*,END=44) keyword 
!!
      if(keyword.eq.'element_nodes_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp
        call checkelement(elementtemp)
        call nuccharge(elementtemp,ztemp)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,layer,node
        if(layer.eq.num_layers_short_atomic(elementindex(ztemp)))then
          write(ounit,*)'Error: do not modifiy the number of output nodes'
          stop !'
        endif
        if(node.gt.maxnodes_short_atomic)then
          write(ounit,*)'Error: too many nodes requested for element ',element(elementindex(ztemp))
          stop 
        endif
        nodes_short_atomic(layer,elementindex(ztemp))=node
!! delete all other activation functions
        do i1=nodes_short_atomic(layer,elementindex(ztemp))+1,maxnodes_short_atomic
          actfunc_short_atomic(i1,layer,elementindex(ztemp))=' '
        enddo
!!
      elseif(keyword.eq.'element_nodes_electrostatic')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp
        call checkelement(elementtemp)
        call nuccharge(elementtemp,ztemp)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,layer,node
        if(layer.eq.num_layers_elec(elementindex(ztemp)))then
          write(ounit,*)'Error: do not modifiy the number of output nodes'
          stop !'
        endif
        if(node.gt.maxnodes_elec)then
          write(ounit,*)'Error: too many nodes requested for element ',element(elementindex(ztemp))
          stop
        endif
        nodes_elec(layer,elementindex(ztemp))=node
!! delete all other activation functions
        do i1=nodes_elec(layer,elementindex(ztemp))+1,maxnodes_elec
          actfunc_elec(i1,layer,elementindex(ztemp))=' '
        enddo
!!
      elseif(keyword.eq.'element_nodes_pair')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2
        call checkelement(elementtemp1)
        call checkelement(elementtemp2)
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,layer,node
        icount=0
        jcount=0
        do i1=1,nelem
          do i2=i1,nelem
            jcount=jcount+1
            if((ztemp1.eq.elempair(jcount,1)).and.(ztemp2.eq.elempair(jcount,2)))then
              icount=jcount
            elseif((ztemp2.eq.elempair(jcount,1)).and.(ztemp1.eq.elempair(jcount,2)))then
              icount=jcount
            endif
          enddo
        enddo
        if(layer.eq.num_layers_short_pair(icount))then
          write(ounit,*)'Error: do not modifiy the number of output nodes'
          stop !'
        endif
        if(node.gt.maxnodes_short_pair)then
          write(ounit,*)'Error: too many nodes requested for element pair ',&
            element(elementindex(elempair(icount,1))),element(elementindex(elempair(icount,2)))
          stop !'
        endif
        nodes_short_pair(layer,icount)=node
!! delete all other activation functions
        do i1=nodes_short_pair(layer,icount)+1,maxnodes_short_pair
          actfunc_short_pair(i1,layer,icount)=' '
        enddo
!!
      elseif(keyword.eq.'element_nodes_ham')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2
        call checkelement(elementtemp1)
        call checkelement(elementtemp2)
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,layer,node
        icount=0
        jcount=0
        do i1=1,nelem
          do i2=i1,nelem
            jcount=jcount+1
            if((ztemp1.eq.elempair(jcount,1)).and.(ztemp2.eq.elempair(jcount,2)))then
              icount=jcount
            elseif((ztemp2.eq.elempair(jcount,1)).and.(ztemp1.eq.elempair(jcount,2)))then
              icount=jcount
            endif
          enddo
        enddo
        if(layer.eq.num_layers_ham(icount))then
          write(ounit,*)'Error: do not modifiy the number of output nodes'
          stop !'
        endif
        if(node.gt.maxnodes_ham)then
          write(ounit,*)'Error: too many nodes requested for element pair ',&
            element(elementindex(elempair(icount,1))),element(elementindex(elempair(icount,2)))
          stop !'
        endif
        nodes_ham(layer,icount)=node
!! delete all other activation functions
        do i1=nodes_ham(layer,icount)+1,maxnodes_ham
          actfunc_ham(i1,layer,icount)=' '
        enddo
!!
      endif
      goto 43
 44   continue
      close(nnunit)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read element-specific activation functions 
!! this must be done after reading the NN architecture  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 45   continue
      read(nnunit,*,END=46) keyword 
!!
      if(keyword.eq.'element_activation_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp
        call checkelement(elementtemp)
        call nuccharge(elementtemp,ztemp)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,layer,node,actfunc
        if(layer.gt.num_layers_short_atomic(elementindex(ztemp)))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_short_atomic(layer,elementindex(ztemp)))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_short_atomic(node,layer,elementindex(ztemp))=actfunc 
!!
      elseif(keyword.eq.'element_activation_electrostatic')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp
        call checkelement(elementtemp)
        call nuccharge(elementtemp,ztemp)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,layer,node,actfunc
        if(layer.gt.num_layers_elec(elementindex(ztemp)))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_elec(layer,elementindex(ztemp)))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_elec(node,layer,elementindex(ztemp))=actfunc
!!
     elseif(keyword.eq.'element_activation_pair')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2
        call checkelement(elementtemp1)
        call checkelement(elementtemp2)
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,&
          layer,node,actfunc
        icount=0
        jcount=0
        do i1=1,nelem
          do i2=i1,nelem
            jcount=jcount+1
            if((ztemp1.eq.elempair(jcount,1)).and.(ztemp2.eq.elempair(jcount,2)))then
              icount=jcount
            elseif((ztemp2.eq.elempair(jcount,1)).and.(ztemp1.eq.elempair(jcount,2)))then
              icount=jcount
            endif
          enddo ! i2
        enddo ! i1
        if(layer.gt.num_layers_short_pair(icount))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_short_pair(layer,icount))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_short_pair(node,layer,icount)=actfunc
!!
     elseif(keyword.eq.'element_activation_ham')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2
        call checkelement(elementtemp1)
        call checkelement(elementtemp2)
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,&
          layer,node,actfunc
        icount=0
        jcount=0
        do i1=1,nelem
          do i2=i1,nelem
            jcount=jcount+1
            if((ztemp1.eq.elempair(jcount,1)).and.(ztemp2.eq.elempair(jcount,2)))then
              icount=jcount
            elseif((ztemp2.eq.elempair(jcount,1)).and.(ztemp1.eq.elempair(jcount,2)))then
              icount=jcount
            endif
          enddo ! i2
        enddo ! i1
        if(layer.gt.num_layers_ham(icount))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_ham(layer,icount))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_ham(node,layer,icount)=actfunc
!!
      endif
      goto 45
 46   continue
      close(nnunit)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read element-specific symmetry functions 
!! this must be done after the determination of elementindex  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(lshort.and.(nn_type_short.eq.1))then
        sym_short_atomic_count(:)=0
        num_funcvalues_short_atomic(:)=0
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        sym_short_pair_count(:)=0
        num_funcvalues_short_pair(:)=0
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        sym_elec_count(:)=0
        num_funcvalues_elec(:)=0
      endif
      if(lnntb)then
        sym_ham_count(:)=0
        num_funcvalues_ham(:)=0
      endif

      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 47   continue
      read(nnunit,*,END=48) keyword 
!!
      if(keyword.eq.'symfunction_short')then
        if(lshort.and.(nn_type_short.eq.1))then
          call readsymfunctionatomic(keyword,&
            maxnum_funcvalues_short_atomic,sym_short_atomic_count,function_type_short_atomic,symelement_short_atomic,&
            funccutoff_short_atomic,eta_short_atomic,zeta_short_atomic,rshift_short_atomic,lambda_short_atomic)
        endif
!!
      elseif(keyword.eq.'element_symfunction_short')then
        if(lshort.and.(nn_type_short.eq.1))then
          call readsymfunctionelementatomic(keyword,&
            maxnum_funcvalues_short_atomic,sym_short_atomic_count,function_type_short_atomic,symelement_short_atomic,&
            funccutoff_short_atomic,eta_short_atomic,zeta_short_atomic,rshift_short_atomic,lambda_short_atomic)
        endif
!!
      elseif(keyword.eq.'global_symfunction_short')then
        if(lshort.and.(nn_type_short.eq.1))then
          call readsymfunctionglobalatomic(keyword,&
            maxnum_funcvalues_short_atomic,sym_short_atomic_count,function_type_short_atomic,symelement_short_atomic,&
            funccutoff_short_atomic,eta_short_atomic,zeta_short_atomic,rshift_short_atomic,lambda_short_atomic)
        endif
!!
      elseif(keyword.eq.'symfunction_electrostatic')then
        if(lelec.and.(nn_type_elec.eq.1))then
          call readsymfunctionatomic(keyword,&
            maxnum_funcvalues_elec,sym_elec_count,function_type_elec,symelement_elec,&
            funccutoff_elec,eta_elec,zeta_elec,rshift_elec,lambda_elec)
        endif
!!
      elseif(keyword.eq.'element_symfunction_electrostatic')then
        if(lelec.and.(nn_type_elec.eq.1))then
          call readsymfunctionelementatomic(keyword,&
            maxnum_funcvalues_elec,sym_elec_count,function_type_elec,symelement_elec,&
            funccutoff_elec,eta_elec,zeta_elec,rshift_elec,lambda_elec)
        endif
!!
      elseif(keyword.eq.'global_symfunction_electrostatic')then
        if(lelec.and.(nn_type_elec.eq.1))then
          call readsymfunctionglobalatomic(keyword,&
            maxnum_funcvalues_elec,sym_elec_count,function_type_elec,symelement_elec,&
            funccutoff_elec,eta_elec,zeta_elec,rshift_elec,lambda_elec)
        endif
!!
      elseif(keyword.eq.'global_pairsymfunction_short')then
        if(lshort.and.(nn_type_short.eq.2))then
          call readsymfunctionglobalpair(keyword,&
            maxnum_funcvalues_short_pair,sym_short_pair_count,function_type_short_pair,symelement_short_pair,&
            funccutoff_short_pair,eta_short_pair,zeta_short_pair,rshift_short_pair,lambda_short_pair)
        endif
!!
      elseif(keyword.eq.'element_pairsymfunction_short')then
        if(lshort.and.(nn_type_short.eq.2))then
          call readsymfunctionelementpair(keyword,&
            maxnum_funcvalues_short_pair,sym_short_pair_count,function_type_short_pair,symelement_short_pair,&
            funccutoff_short_pair,eta_short_pair,zeta_short_pair,rshift_short_pair,lambda_short_pair)
        endif
!!
      elseif(keyword.eq.'pairsymfunction_short')then
        if(lshort.and.(nn_type_short.eq.2))then
          call readsymfunctionpair(keyword,&
            maxnum_funcvalues_ham,sym_ham_count,function_type_ham,symelement_ham,&
            funccutoff_ham,eta_ham,zeta_ham,rshift_ham,lambda_ham)
        endif
!!
      elseif(keyword.eq.'global_hamsymfunction')then
        if(lnntb)then
          call readsymfunctionglobalpair(keyword,&
            maxnum_funcvalues_short_pair,sym_short_pair_count,function_type_short_pair,symelement_short_pair,&
            funccutoff_short_pair,eta_short_pair,zeta_short_pair,rshift_short_pair,lambda_short_pair)
        endif
!!
      elseif(keyword.eq.'element_hamsymfunction')then
        if(lnntb)then
          call readsymfunctionelementpair(keyword,&
            maxnum_funcvalues_ham,sym_ham_count,function_type_ham,symelement_ham,&
            funccutoff_ham,eta_ham,zeta_ham,rshift_ham,lambda_ham)
        endif
!!
      elseif(keyword.eq.'hamsymfunction')then
        if(lnntb)then
          call readsymfunctionpair(keyword,&
            maxnum_funcvalues_ham,sym_ham_count,function_type_ham,symelement_ham,&
            funccutoff_ham,eta_ham,zeta_ham,rshift_ham,lambda_ham)
        endif
!!
      endif ! keyword
      goto 47
 48   continue
      close(nnunit)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! determination of the number of symmetry functions for each element and elementalpair
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          num_funcvalues_short_atomic(i1)=sym_short_atomic_count(i1)
          nodes_short_atomic(0,i1)=num_funcvalues_short_atomic(i1)
        endif
        if(lelec.and.(nn_type_elec.eq.1))then
          num_funcvalues_elec(i1)=sym_elec_count(i1)
          nodes_elec(0,i1)=num_funcvalues_elec(i1)
        endif
      enddo
      do i1=1,npairs
        if(lshort.and.(nn_type_short.eq.2))then
          num_funcvalues_short_pair(i1)=sym_short_pair_count(i1)
          nodes_short_pair(0,i1)=num_funcvalues_short_pair(i1)
        endif
        if(lnntb)then
          num_funcvalues_ham(i1)=sym_ham_count(i1)
          nodes_ham(0,i1)=num_funcvalues_ham(i1)
        endif
      enddo
!! do some checks
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_short_atomic(i1).eq.0)then
            write(ounit,*)'ERROR: No short range symfunctions specified for ',element(i1)
            stop
          endif
        enddo
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_elec(i1).eq.0)then
            write(ounit,*)'ERROR: No electrostatic symfunctions specified for ',element(i1)
            stop
          endif
        enddo
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        icount=0
        do i1=1,nelem
          do i2=i1,nelem
            icount=icount+1
            if(num_funcvalues_short_pair(icount).eq.0)then
              write(ounit,*)'WARNING: No short range pair symfunctions specified for ',element(i1),element(i2)
!!              stop !'
            endif
          enddo
        enddo
      endif
      if(nn_type_nntb.eq.1)then
        icount=0
        do i1=1,nelem
          do i2=i1,nelem
            icount=icount+1
            if(num_funcvalues_ham(icount).eq.0)then
              write(ounit,*)'WARNING: No hamiltonian symfunctions specified for ',element(i1),element(i2)
              stop !'
            endif
          enddo
        enddo
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read basis definitions for nn_type_nntb 1 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(nn_type_nntb.eq.1)then
        call readbasis()
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check for inconsistencies in the input.nn file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      call checkinputnn()
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate derived quantities 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! calculate number of weight parameters in short range NN type 1
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          wcount=0
          do i=1,num_layers_short_atomic(i1) ! loop over all hidden and output layers
            wcount=wcount+1
            windex_short_atomic(wcount,i1)=num_weights_short_atomic(i1)+1
            num_weights_short_atomic(i1)=num_weights_short_atomic(i1)&
              +nodes_short_atomic(i-1,i1)*nodes_short_atomic(i,i1)
            wcount=wcount+1
            windex_short_atomic(wcount,i1)=num_weights_short_atomic(i1)+1
            num_weights_short_atomic(i1)=num_weights_short_atomic(i1)&
              +nodes_short_atomic(i,i1) ! bias weights
          enddo
          maxnum_weights_short_atomic=max(maxnum_weights_short_atomic,num_weights_short_atomic(i1))
        enddo ! i1'
      endif ! lshort
!!
!! calculate number of weight parameters in electrostatic NN
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          wcount=0
          do i=1,num_layers_elec(i1)
            wcount=wcount+1
            windex_elec(wcount,i1)=num_weights_elec(i1)+1
            num_weights_elec(i1)=num_weights_elec(i1)+nodes_elec(i-1,i1)*nodes_elec(i,i1)
            wcount=wcount+1
            windex_elec(wcount,i1)=num_weights_elec(i1)+1
            num_weights_elec(i1)=num_weights_elec(i1)+nodes_elec(i,i1) ! bias weights
          enddo
          write(ounit,'(a,a3,i10)')' => electrostatic NN weights                     ',element(i1),num_weights_elec(i1)
          maxnum_weights_elec=max(maxnum_weights_elec,num_weights_elec(i1))
        enddo ! i1 '
      endif ! 
!!
!! calculate number of weight parameters in short range NN type 2
      if(lshort.and.(nn_type_short.eq.2))then
        do i1=1,npairs
          wcount=0
          do i=1,num_layers_short_pair(i1) ! loop over all hidden and output layers
            wcount=wcount+1
            windex_short_pair(wcount,i1)=num_weights_short_pair(i1)+1
            num_weights_short_pair(i1)=num_weights_short_pair(i1)+nodes_short_pair(i-1,i1)*nodes_short_pair(i,i1)
            wcount=wcount+1
            windex_short_pair(wcount,i1)=num_weights_short_pair(i1)+1
            num_weights_short_pair(i1)=num_weights_short_pair(i1)+nodes_short_pair(i,i1) ! bias weights
          enddo
          maxnum_weights_short_pair=max(maxnum_weights_short_pair,num_weights_short_pair(i1))
        enddo ! i1'
      endif ! lshort
!!
!! calculate number of weight parameters in short range NN type 2
      if(nn_type_nntb.eq.1)then
        do i1=1,npairs
          wcount=0
          do i=1,num_layers_ham(i1) ! loop over all hidden and output layers
            wcount=wcount+1
            windex_ham(wcount,i1)=num_weights_ham(i1)+1
            num_weights_ham(i1)=num_weights_ham(i1)+nodes_ham(i-1,i1)*nodes_ham(i,i1)
            wcount=wcount+1
            windex_ham(wcount,i1)=num_weights_ham(i1)+1
            num_weights_ham(i1)=num_weights_ham(i1)+nodes_ham(i,i1) ! bias weights
          enddo
          write(ounit,'(a,2a3,i10)')' => hamiltonian NN weights                       ',&
            element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),num_weights_ham(i1)
          maxnum_weights_ham=max(maxnum_weights_ham,num_weights_ham(i1))
        enddo ! i1'
      endif 
!!
!! set dummy dimensions for unused arrays
      if(nn_type_short.eq.1)then
        maxnum_weights_short_pair=1
        maxnum_weights_ham =1
      endif
      if(nn_type_short.eq.2)then
        maxnum_weights_short_atomic=1
        maxnum_weights_ham  =1
      endif
      if((.not.lelec).or.(lelec.and.(nn_type_elec.ne.1)))then
        maxnum_weights_elec=1
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read atom reference energies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(lremoveatomenergies)then
        call readatomenergies() 
      endif
!!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read activation functions of individual nodes 
!! WHY IS THIS COMMENTED???
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        open(nnunit,file='input.nn',form='formatted',status='old')
!!        rewind(nnunit)
!!        icount=0
!! 33     continue       
!!          read(nnunit,*,END=34)keyword
!!          if(keyword.eq.'node_activation_short')then
!!            backspace(nnunit)
!!            read(nnunit,*,ERR=99)dummy,elementtemp,layer,node
!!            call nuccharge(elementtemp,ztemp)
!!            if(layer.gt.num_layers_short_atomic(elementindex(ztemp)))then
!!              write(ounit,*)'Error: layer is too large in ',keyword
!!              stop
!!            endif
!!            if(node.gt.nodes_short_atomic(layer,elementindex(ztemp)))then
!!              write(ounit,*)'Error: node is too large in ',keyword
!!              stop
!!            endif
!!            backspace(nnunit)
!!            read(nnunit,*,ERR=99)dummy,element,layer,node,actfunc_short(node,layer,elementindex(ztemp))
!!          elseif(keyword.eq.'node_activation_electrostatic')then
!!            backspace(nnunit)
!!            read(nnunit,*,ERR=99)dummy,element,layer,node
!!            call nuccharge(elementtemp,ztemp)
!!            if(layer.gt.num_layersewald(elementindex(ztemp)))then
!!              write(ounit,*)'Error: layer is too large in ',keyword
!!              stop
!!            endif
!!            if(node.gt.nodes_ewald(layer,elementindex(ztemp)))then
!!              write(ounit,*)'Error: node is too large in ',keyword
!!              stop
!!            endif
!!            backspace(nnunit)
!!            read(nnunit,*,ERR=99)dummy,element,layer,node,actfunc_ewald(node,layer,elementindex(ztemp))
!!         endif
!!          goto 33
!! 34     continue      
!!        close(nnunit)
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! sort and analyze symmetry functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! sort short range atomic symmetry functions
      if((nn_type_short.eq.1).and.lshort)then
        call sortsymfunctions(&
          maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
          function_type_short_atomic,symelement_short_atomic,&
          eta_short_atomic,zeta_short_atomic,rshift_short_atomic,lambda_short_atomic,funccutoff_short_atomic)
      endif
!!
!! sort electrostatic symmetry functions
      if(lelec.and.(nn_type_elec.eq.1))then
        call sortsymfunctions(&
          maxnum_funcvalues_elec,num_funcvalues_elec,&
          function_type_elec,symelement_elec,&
          eta_elec,zeta_elec,rshift_elec,lambda_elec,funccutoff_elec)
      endif
!!
!! sort short range pair symmetry functions
      if((nn_type_short.eq.2).and.lshort)then
        call sortpairsymfunctions(npairs,&
          maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
          function_type_short_pair,symelement_short_pair,&
          eta_short_pair,zeta_short_pair,rshift_short_pair,lambda_short_pair,funccutoff_short_pair,&
          elempair)
      endif
!!
!! sort short range pair symmetry functions
      if(nn_type_nntb.eq.1)then
        call sortpairsymfunctions(npairs,&
          maxnum_funcvalues_ham,num_funcvalues_ham,&
          function_type_ham,symelement_ham,&
          eta_ham,zeta_ham,rshift_ham,lambda_ham,funccutoff_ham,&
          elempair)
      endif
!!
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
