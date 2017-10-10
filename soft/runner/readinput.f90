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
      integer tripletid                                             ! internal
      integer ndim                                                  ! internal
      integer nodes_short_atomic_temp(0:maxnum_layers_short_atomic) ! internal
      integer nodes_elec_temp(0:maxnum_layers_elec)                 ! internal
      integer nodes_short_pair_temp(0:maxnum_layers_short_pair)     ! internal
      integer nodes_ham_temp(0:maxnum_layers_ham)                   ! internal
      integer nodes_s_temp(0:maxnum_layers_s)                   ! internal
      integer nodes_hexton_temp(0:maxnum_layers_hexton)                   ! internal
      integer nodes_hextoff_temp(0:maxnum_layers_hextoff)                   ! internal
      integer nodes_dens_temp(0:maxnum_layers_dens)                   ! internal
      integer iseed                                     ! out, seed for weight initialization
      integer ielem                                     ! in (read also before from input.nn) 
      integer wcount                                    ! internal 
      integer i,j,k
      integer i0,i1,i2,i3                               ! internal
      integer ztemp                                     ! internal
      integer ztemp1                                    ! internal
      integer ztemp2                                    ! internal
      integer itemp                                     ! internal
      integer layer                                     ! internal
      integer node                                      ! internal
      integer maxnodes_hextoff_old                      ! internal
      integer sym_short_atomic_count(nelem)             ! internal
      integer sym_elec_count(nelem)                     ! internal
      integer sym_short_pair_count(npairs)              ! internal
      integer sym_ham_count(npairs)                     ! internal
      integer sym_s_count(npairs)                       ! internal
      integer sym_hexton_count(npairs)                  ! internal
      integer, allocatable :: sym_hextoff_count(:)              ! internal
      integer sym_dens_count(npairs)                    ! internal
      integer counter                                   ! internal
!!
      real*8 kalmanlambda_local                         ! internal
      real*8 kalmanlambdae_local                        ! internal
      real*8 kalmanlambda_ham_local                     ! internal CMH Kalman filter needs to be covered
      real*8 chargetemp                                 ! internal
!!
      character*40 dummy                                ! internal
      character*40 keyword                              ! internal
      character*1 actfunc_short_atomic_dummy(maxnum_layers_short_atomic)    ! internal
      character*1 actfunc_elec_dummy(maxnum_layers_elec)                    ! internal
      character*1 actfunc_short_pair_dummy(maxnum_layers_short_pair)        ! internal
      character*1 actfunc_ham_dummy(maxnum_layers_ham)                      ! internal
      character*1 actfunc_s_dummy(maxnum_layers_s)                          ! internal
      character*1 actfunc_hexton_dummy(maxnum_layers_hexton)                ! internal
      character*1 actfunc_hextoff_dummy(maxnum_layers_hextoff)              ! internal
      character*1, dimension (:,:,:), allocatable :: actfunc_hextoff_trans    ! internal
      character*1 actfunc_dens_dummy(maxnum_layers_dens)                    ! internal
      character*1 actfunc                               ! internal
      character*2 elementtemp                           ! internal
      character*2 elementtemp1                          ! internal
      character*2 elementtemp2                          ! internal
!!
      logical lelement(102)                             ! internal
      logical lprint                                    ! internal
!!
!! initializations
!! all inputs must be initialized here, except for the ones read before in getdimensions.f90
      write(ounit,*)'Reading control parameters from input.nn'
      write(ounit,*)'============================================================='
!!'
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
!        nodes_ham_temp(:)            =0
!        actfunc_ham_dummy(:)         =' '
        if(nntb_flag(1)) then
          nodes_s_temp(:)            =0
          actfunc_s_dummy(:)         =' '
        endif
        if(nntb_flag(2)) then
          nodes_hexton_temp(:)            =0
          actfunc_hexton_dummy(:)         =' '
        endif
        if(nntb_flag(3)) then
          nodes_hextoff_temp(:)            =0
          actfunc_hextoff_dummy(:)         =' '
        endif
        if(nntb_flag(4)) then
          nodes_dens_temp(:)            =0
          actfunc_dens_dummy(:)         =' '
        endif

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
!        windex_ham(:,:)             =0
!        num_weights_ham(:)          =0 
!        maxnum_weights_ham          =0
        if(nntb_flag(1)) then
          windex_s(:,:)             =0
          num_weights_s(:)          =0
          maxnum_weights_s          =0
        endif
        if(nntb_flag(2)) then
          windex_hexton(:,:)             =0
          num_weights_hexton(:)          =0
          maxnum_weights_hexton          =0
        endif
        if(nntb_flag(3)) then
          if(mode.ge.2) then
            windex_hextoff(:,:)             =0
            num_weights_hextoff(:)          =0
            maxnum_weights_hextoff          =0
          endif
        endif
        if(nntb_flag(4)) then
          windex_dens(:,:)             =0
          num_weights_dens(:)          =0
          maxnum_weights_dens          =0
        endif

      endif
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! start regular reading of input.nn '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!

      
      call readkeywords(iseed,&
        nodes_short_atomic_temp,nodes_elec_temp,nodes_short_pair_temp,nodes_ham_temp,&
        nodes_s_temp,nodes_hexton_temp,nodes_hextoff_temp,nodes_dens_temp,&
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
        if(nntb_flag(0)) then
          do i1=1,npairs
            nodes_ham(maxnum_layers_ham,i1)=1
          enddo
        endif
        if(nntb_flag(1))then
          do i1=1,npairs
            nodes_s(maxnum_layers_s,i1)=1
          enddo
        endif
        if(nntb_flag(2))then
          do i1=1,npairs
            nodes_hexton(maxnum_layers_hexton,i1)=1
          enddo
        endif
        if(nntb_flag(3))then
          if(mode.eq.1)then
            i1 = 1
            allocate(sym_hextoff_count(1))
          elseif(mode.eq.2)then
            i1 = 1
            
            nodes_hextoff(maxnum_layers_hextoff,i1)=1
            allocate(sym_hextoff_count(1))
          elseif(mode.eq.3)then
            allocate(sym_hextoff_count(ntriplets))
            do i1=1,ntriplets
              nodes_hextoff(maxnum_layers_hextoff,i1)=1
            enddo
          endif
        endif
        if(nntb_flag(4))then
          do i1=1,npairs
            nodes_dens(maxnum_layers_dens,i1)=1
          enddo
        endif

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
      elseif((keyword.eq.'global_activation_s').and.lnntb.and.nntb_flag(1))then
        if(count_global_activation_s.gt.0)then
          write(ounit,*)'Error: global_activation_s specified twice'
          stop
        endif
        call setglobalactivation(npairs,count_global_activation_s,&
          maxnum_layers_s,maxnodes_s,nodes_s,actfunc_s,&
          actfunc_s_dummy,keyword)
!!
      elseif((keyword.eq.'global_activation_hexton').and.lnntb.and.nntb_flag(2))then
        if(count_global_activation_hexton.gt.0)then
          write(ounit,*)'Error: global_activation_hexton specified twice'
          stop
        endif
        call setglobalactivation(npairs,count_global_activation_hexton,&
          maxnum_layers_hexton,maxnodes_hexton,nodes_hexton,actfunc_hexton,&
          actfunc_hexton_dummy,keyword)
!!
      elseif((keyword.eq.'global_activation_hextoff').and.lnntb.and.nntb_flag(3))then
        if(count_global_activation_hextoff.gt.0)then
          write(ounit,*)'Error: global_activation_hextoff specified twice'
          stop
        endif
        if(mode.eq.2)then
          call setglobalactivation(1,count_global_activation_hextoff,&
            maxnum_layers_hextoff,maxnodes_hextoff,nodes_hextoff,actfunc_hextoff,&
            actfunc_hextoff_dummy,keyword)
        elseif(mode.eq.3)then
          call setglobalactivation(ntriplets,count_global_activation_hextoff,&
            maxnum_layers_hextoff,maxnodes_hextoff,nodes_hextoff,actfunc_hextoff,&
            actfunc_hextoff_dummy,keyword)
        endif
!!
      elseif((keyword.eq.'global_activation_dens').and.lnntb.and.nntb_flag(4))then
        if(count_global_activation_dens.gt.0)then
          write(ounit,*)'Error: global_activation_dens specified twice'
          stop
        endif
        call setglobalactivation(npairs,count_global_activation_dens,&
          maxnum_layers_dens,maxnodes_dens,nodes_dens,actfunc_ham,&
          actfunc_dens_dummy,keyword)

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
      icount = 0
      do i=1,nelem
        do j=i,nelem
          do k=1,nelem
            icount=icount+1
            elemtriplet(icount,1)=nucelem(i)
            elemtriplet(icount,2)=nucelem(j)
            elemtriplet(icount,3)=nucelem(k)
          enddo
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
      elseif(keyword.eq.'element_hidden_layers_s')then
        call readelementlayerspair(maxnodes_s,&
          maxnum_layers_s,num_layers_s,nodes_ham,actfunc_s)
!!
      elseif(keyword.eq.'element_hidden_layers_hexton')then
        call readelementlayerspair(maxnodes_hexton,&
          maxnum_layers_hexton,num_layers_hexton,nodes_ham,actfunc_hexton)
!!
      elseif(keyword.eq.'element_hidden_layers_hextoff')then
        call readelementlayerspair(maxnodes_hextoff,&
          maxnum_layers_hextoff,num_layers_hextoff,nodes_hextoff,actfunc_hextoff)
!!
      elseif(keyword.eq.'element_hidden_layers_dens')then
        call readelementlayerspair(maxnodes_dens,&
          maxnum_layers_dens,num_layers_dens,nodes_dens,actfunc_dens)
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
!!
     elseif(keyword.eq.'element_activation_s')then
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
        if(layer.gt.num_layers_s(icount))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_s(layer,icount))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_s(node,layer,icount)=actfunc
!!
!!
     elseif(keyword.eq.'element_activation_hexton')then
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
        if(layer.gt.num_layers_hexton(icount))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_hexton(layer,icount))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_hexton(node,layer,icount)=actfunc
!!
!!
     elseif(keyword.eq.'element_activation_hextoff')then
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
        if(layer.gt.num_layers_hextoff(icount))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_hextoff(layer,icount))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_hextoff(node,layer,icount)=actfunc
!!
!!
     elseif(keyword.eq.'element_activation_dens')then
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
        if(layer.gt.num_layers_dens(icount))then
          write(ounit,*)'Error: layer is too large in ',keyword
          stop
        endif
        if(node.gt.nodes_dens(layer,icount))then
          write(ounit,*)'Error: node is too large in ',keyword
          stop
        endif
        actfunc_dens(node,layer,icount)=actfunc
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
        if(nntb_flag(0)) then
        sym_ham_count(:)=0
        num_funcvalues_ham(:)=0
        endif
        if(nntb_flag(1))then
        sym_s_count(:)=0
        num_funcvalues_s(:)=0
        endif
        if(nntb_flag(2))then
        sym_hexton_count(:)=0
        num_funcvalues_hexton(:)=0
        endif
        if(nntb_flag(3))then
          sym_hextoff_count(:)=0
          num_funcvalues_hextoff(:)=0
       endif
        if(nntb_flag(4))then
        sym_dens_count(:)=0
        num_funcvalues_dens(:)=0
        endif

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
      elseif(keyword.eq.'global_symfunction_s')then
        if(lnntb.and.nntb_flag(1))then
          call readsymfunctionglobalpair(keyword,&
            maxnum_funcvalues_s,sym_s_count,function_type_s,symelement_s,&
            funccutoff_s,eta_s,zeta_s,rshift_s,lambda_s)
        endif
!!
      elseif(keyword.eq.'global_symfunction_hexton')then
        if(lnntb.and.nntb_flag(2))then
          call readsymfunctionglobalpair(keyword,&
            maxnum_funcvalues_hexton,sym_hexton_count,function_type_hexton,symelement_hexton,&
            funccutoff_hexton,eta_hexton,zeta_hexton,rshift_hexton,lambda_hexton)
        endif
!!
      elseif(keyword.eq.'symfunction_hextoff')then
        if(lnntb.and.nntb_flag(3))then
          if(mode.le.2)then
            ndim=1
          else
            ndim=ntriplets
          endif
          call readsymfunctionhextoff(keyword,ndim,&
            maxnum_funcvalues_hextoff,sym_hextoff_count,function_type_hextoff,symelement_hextoff,&
            funccutoff_hextoff,eta_hextoff,zeta_hextoff,rshift_hextoff,lambda_hextoff)
        endif
!!
      elseif(keyword.eq.'global_symfunction_dens')then
        if(lnntb.and.nntb_flag(4))then
          call readsymfunctionglobalpair(keyword,&
            maxnum_funcvalues_dens,sym_dens_count,function_type_dens,symelement_dens,&
            funccutoff_dens,eta_dens,zeta_dens,rshift_dens,lambda_dens)
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
      enddo
      if(lnntb.and.nntb_flag(3))then
        if(mode.eq.3)then
          do i1=1,ntriplets
            num_funcvalues_hextoff(i1)=sym_hextoff_count(i1)
            nodes_hextoff(0,i1)=num_funcvalues_hextoff(i1)
          enddo
        else ! modes 1 or 2
          num_funcvalues_hextoff(1)=sym_hextoff_count(1)
          if(mode.eq.2)then
            nodes_hextoff(0,1)=num_funcvalues_hextoff(1)
          endif
        endif
      endif

      do i1=1,npairs
        if(lnntb)then
          if(nntb_flag(0)) then
            num_funcvalues_ham(i1)=sym_ham_count(i1)
            nodes_ham(0,i1)=num_funcvalues_ham(i1)
          endif
          if(nntb_flag(1))then
            num_funcvalues_s(i1)=sym_s_count(i1)
            nodes_s(0,i1)=num_funcvalues_s(i1)
          endif
          if(nntb_flag(2))then
            num_funcvalues_hexton(i1)=sym_hexton_count(i1)
            nodes_hexton(0,i1)=num_funcvalues_hexton(i1)
          endif
!          if(nntb_flag(3))then
!            num_funcvalues_hextoff(i1)=sym_hextoff_count(i1)
!            nodes_hextoff(0,i1)=num_funcvalues_hextoff(i1)
!          endif
          if(nntb_flag(4))then
            num_funcvalues_dens(i1)=sym_dens_count(i1)
            nodes_dens(0,i1)=num_funcvalues_dens(i1)
          endif
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
      if(lnntb)then
        if(nntb_flag(0)) then
        icount=0
        do i1=1,nelem
          do i2=i1,nelem
            icount=icount+1
            if(num_funcvalues_ham(icount).eq.0)then
              write(ounit,*)'WARNING: No hamiltonian symfunctions specified for ',element(i1),element(i2)
!              stop !'
            endif
          enddo
        enddo
        endif
        if(nntb_flag(1))then
        icount=0
        do i1=1,nelem
          do i2=i1,nelem
            icount=icount+1
            if(num_funcvalues_s(icount).eq.0)then
              write(ounit,*)'WARNING: No hamiltonian overlap symfunctions specified for ',element(i1),element(i2)
              stop !'
            endif
          enddo
        enddo
        endif
        if(nntb_flag(2))then
        icount=0
        do i1=1,nelem
          do i2=i1,nelem
            icount=icount+1
            if(num_funcvalues_hexton(icount).eq.0)then
              write(ounit,*)'WARNING: No hamiltonian hexton symfunctions specified for ',element(i1),element(i2)
              stop !'
            endif
          enddo
        enddo
        endif
        if(nntb_flag(4))then
          icount=0
          do i1=1,nelem
            do i2=i1,nelem
              icount=icount+1
              if(num_funcvalues_dens(icount).eq.0)then
                write(ounit,*)'WARNING: No hamiltonian dens symfunctions specified for ',element(i1),element(i2)
                stop !'
              endif
            enddo
          enddo
        endif

      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read basis definitions for nn_type_nntb 1 this is now done earlier in getdimensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(lnntb)then
         call readbasis()
         if(nntb_flag(3))then
           if((mode.eq.1).or.(mode.eq.2))then
             nodes_hextoff(maxnum_layers_hextoff,1) = &
               num_basis(elementindex(hextoff_training_triplet(1)))*&
               num_basis(elementindex(hextoff_training_triplet(2)))
             maxnodes_hextoff_old = maxnodes_hextoff
             maxnodes_hextoff=max(maxnodes_hextoff,nodes_hextoff(maxnum_layers_hextoff,1))
             allocate(actfunc_hextoff_trans(maxnodes_hextoff,maxnum_layers_hextoff,1))
             actfunc_hextoff_trans = ''
             actfunc_hextoff_trans(1:maxnodes_hextoff_old,1:maxnum_layers_hextoff,1)&
               =actfunc_hextoff(1:maxnodes_hextoff_old,1:maxnum_layers_hextoff,1)
             deallocate(actfunc_hextoff)
             allocate(actfunc_hextoff(maxnodes_hextoff,maxnum_layers_hextoff,1))
             actfunc_hextoff = actfunc_hextoff_trans
             actfunc_hextoff(1:maxnodes_hextoff,maxnum_layers_hextoff,1) = actfunc_hextoff(1,maxnum_layers_hextoff,1)    
             deallocate(actfunc_hextoff_trans)
           elseif(mode.eq.3)then
             maxnodes_hextoff_old = maxnodes_hextoff
             do i1=1,nelem
               do i2=1,nelem
                 if(nucelem(i2).ge.nucelem(i1))then
                   do i3=1,nelem
                     i = tripletindex(nucelem(i1),nucelem(i2),nucelem(i3))
                     nodes_hextoff(maxnum_layers_hextoff,i) = &
                       num_basis(elementindex(nucelem(i1)))*&
                       num_basis(elementindex(nucelem(i2)))
                     maxnodes_hextoff=max(maxnodes_hextoff,nodes_hextoff(maxnum_layers_hextoff,i))
                   enddo
                 endif
               enddo
             enddo
             allocate(actfunc_hextoff_trans(maxnodes_hextoff,maxnum_layers_hextoff,ntriplets))
             actfunc_hextoff_trans = ''
             actfunc_hextoff_trans(1:maxnodes_hextoff_old,1:maxnum_layers_hextoff,:)&
               =actfunc_hextoff(1:maxnodes_hextoff_old,1:maxnum_layers_hextoff,:)
             deallocate(actfunc_hextoff)
             allocate(actfunc_hextoff(maxnodes_hextoff,maxnum_layers_hextoff,ntriplets))
             actfunc_hextoff = actfunc_hextoff_trans
             do i1=1,ntriplets
               actfunc_hextoff(1:maxnodes_hextoff,maxnum_layers_hextoff,i1)& 
               =actfunc_hextoff(1,maxnum_layers_hextoff,i1)
             enddo
             deallocate(actfunc_hextoff_trans)
           endif
         else
           write(*,*) 'Error: nntb_flag modes 1,2,4 not implemented'
           stop
         endif
      endif
!      write(*,*) maxnum_basis
!      write(*,*) num_basis(1)
!      write(*,*) num_basis(2)
!      write(*,*) basis(1,:,:)
!      write(*,*) basis(2,:,:)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check for inconsistencies in the input.nn file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      call checkinputnn()
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write settings of input.nn/defaults to runner.out:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      call printinputnn(iseed,ielem,&
       nodes_short_atomic_temp,nodes_elec_temp,nodes_short_pair_temp,nodes_ham_temp,&
       nodes_s_temp,nodes_hexton_temp,nodes_hextoff_temp,nodes_dens_temp,&
       kalmanlambda_local,kalmanlambdae_local,& !! KALMAN FILTER will need expanding with Hamiltonian variants
       actfunc_short_atomic_dummy,actfunc_elec_dummy,actfunc_short_pair_dummy,actfunc_ham_dummy,&
       actfunc_s_dummy,actfunc_hexton_dummy,actfunc_hextoff_dummy,actfunc_dens_dummy)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Write the unique element pairs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(ounit,*)'Element pairs: ',npairs
      icount=0
      do i=1,nelem
        do j=i,nelem
          icount=icount+1
          write(ounit,'(a6,i4,2a3)')' pair ',icount,element(i),element(j)
        enddo
      enddo
      write(ounit,*)'============================================================='
!!'
!! CMH Write out the single triplet being used for Hextoff training. This would also be done for the single
!! pair being used for Overlap and Hexton training
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
          write(ounit,'(a,a3,i10)')' => short range NN weights type 1                ',element(i1),num_weights_short_atomic(i1)
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
          write(ounit,'(a,2a3,i10)')' => short range NN weights type 2                ',&
            element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),num_weights_short_pair(i1)
          maxnum_weights_short_pair=max(maxnum_weights_short_pair,num_weights_short_pair(i1))
        enddo ! i1'
      endif ! lshort
!!
!! calculate number of weight parameters in short range NN type 2
      
      if(lnntb.and.(mode.gt.1))then
        if(nntb_flag(0)) then
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
        if(nntb_flag(1)) then
        do i1=1,npairs
          wcount=0
          do i=1,num_layers_s(i1) ! loop over all hidden and output layers
            wcount=wcount+1
            windex_s(wcount,i1)=num_weights_s(i1)+1
            num_weights_s(i1)=num_weights_s(i1)+nodes_s(i-1,i1)*nodes_s(i,i1)
            wcount=wcount+1
            windex_s(wcount,i1)=num_weights_s(i1)+1
            num_weights_s(i1)=num_weights_s(i1)+nodes_s(i,i1) ! bias weights
          enddo
          write(ounit,'(a,2a3,i10)')' => hamiltonian Overlap  NN weights              ',&
            element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),num_weights_s(i1)
          maxnum_weights_s=max(maxnum_weights_s,num_weights_s(i1))
        enddo ! i1'
        endif
        if(nntb_flag(2)) then
        do i1=1,npairs
          wcount=0
          do i=1,num_layers_hexton(i1) ! loop over all hidden and output layers
            wcount=wcount+1
            windex_hexton(wcount,i1)=num_weights_hexton(i1)+1
            num_weights_hexton(i1)=num_weights_hexton(i1)+nodes_hexton(i-1,i1)*nodes_hexton(i,i1)
            wcount=wcount+1
            windex_hexton(wcount,i1)=num_weights_hexton(i1)+1
            num_weights_hexton(i1)=num_weights_hexton(i1)+nodes_hexton(i,i1) ! bias weights
          enddo
          write(ounit,'(a,2a3,i10)')' => hamiltonian Hextonsite NN weights            ',&
            element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),num_weights_hexton(i1)
          maxnum_weights_hexton=max(maxnum_weights_hexton,num_weights_hexton(i1))
        enddo ! i1'
        endif
        if(nntb_flag(3)) then
          if((mode.eq.1).or.(mode.eq.2)) then
            i1 = tripletindex(hextoff_training_triplet(1),hextoff_training_triplet(2),hextoff_training_triplet(3))
            wcount=0
            do i=1,num_layers_hextoff(1) ! loop over all hidden and output layers
              wcount=wcount+1
              windex_hextoff(wcount,1)=num_weights_hextoff(1)+1
              num_weights_hextoff(1)=num_weights_hextoff(1)+nodes_hextoff(i-1,1)*nodes_hextoff(i,1)
              wcount=wcount+1
              windex_hextoff(wcount,1)=num_weights_hextoff(1)+1
              num_weights_hextoff(1)=num_weights_hextoff(1)+nodes_hextoff(i,1) ! bias weights
            enddo
            
            write(ounit,'(a,3a3,i10)')' => hamiltonian Hextoff  NN weights            ',&
              element(elementindex(hextoff_training_triplet(1))),&
              element(elementindex(hextoff_training_triplet(2))),&
              element(elementindex(hextoff_training_triplet(3))),&
                      num_weights_hextoff(1)
            maxnum_weights_hextoff=max(maxnum_weights_hextoff,num_weights_hextoff(1))
          else
            do i1=1,ntriplets
              wcount=0
              do i=1,num_layers_hextoff(i1) ! loop over all hidden and output layers
                wcount=wcount+1
                windex_hextoff(wcount,i1)=num_weights_hextoff(i1)+1
                num_weights_hextoff(i1)=num_weights_hextoff(i1)+nodes_hextoff(i-1,i1)*nodes_hextoff(i,i1)
                wcount=wcount+1
                windex_hextoff(wcount,i1)=num_weights_hextoff(i1)+1
                num_weights_hextoff(i1)=num_weights_hextoff(i1)+nodes_hextoff(i,i1) ! bias weights
              enddo
              write(ounit,'(a,3a3,i10)')' => hamiltonian Hextoff  NN weights            ',&
              element(elementindex(elemtriplet(i1,1))),element(elementindex(elemtriplet(i1,2))),&
                element(elementindex(elemtriplet(i1,3))),num_weights_hextoff(i1)
              maxnum_weights_hextoff=max(maxnum_weights_hextoff,num_weights_hextoff(i1))
            enddo ! i1'
          endif
        endif
        if(nntb_flag(4)) then
        do i1=1,npairs
          wcount=0
          do i=1,num_layers_dens(i1) ! loop over all hidden and output layers
            wcount=wcount+1
            windex_dens(wcount,i1)=num_weights_dens(i1)+1
            num_weights_dens(i1)=num_weights_dens(i1)+nodes_dens(i-1,i1)*nodes_dens(i,i1)
            wcount=wcount+1
            windex_dens(wcount,i1)=num_weights_dens(i1)+1
            num_weights_dens(i1)=num_weights_dens(i1)+nodes_dens(i,i1) ! bias weights
          enddo
          write(ounit,'(a,2a3,i10)')' => hamiltonian Density  NN weights             ',&
            element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),num_weights_dens(i1)
          maxnum_weights_dens=max(maxnum_weights_dens,num_weights_dens(i1))
        enddo ! i1'
        endif
      endif 
!!
      write(ounit,*)'-------------------------------------------------------------'
!!'
!! set dummy dimensions for unused arrays
      if(nn_type_short.eq.1)then
        maxnum_weights_short_pair=1
        if(lnntb) then
          if(.not.nntb_flag(0)) then
            maxnum_weights_ham =1
          endif
          if(.not.nntb_flag(1)) then
            maxnum_weights_s =1
          endif
          if(.not.nntb_flag(2)) then
            maxnum_weights_hexton =1
          endif
          if(.not.nntb_flag(3)) then
            maxnum_weights_hextoff =1
          endif
          if(.not.nntb_flag(4)) then
            maxnum_weights_dens =1
          endif

        endif
      endif
      if(nn_type_short.eq.2)then
        maxnum_weights_short_atomic=1
        if(lnntb) then
          if(.not.nntb_flag(0)) then
            maxnum_weights_ham =1
          endif
          if(.not.nntb_flag(1)) then
            maxnum_weights_s =1
          endif
          if(.not.nntb_flag(2)) then
            maxnum_weights_hexton =1
          endif
          if(.not.nntb_flag(3)) then
            maxnum_weights_hextoff =1
          endif
          if(.not.nntb_flag(4)) then
            maxnum_weights_dens =1
          endif
        endif
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! print NN architectures 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(lshort.and.(nn_type_short.eq.1).and.(mode.ne.1))then
          do i3=1,nelem
            write(ounit,*)'-------------------------------------------------'
            write(ounit,*)'Atomic short range NN for element: ',element(i3)
            write(ounit,'(a,10i5)')' architecture    ',(nodes_short_atomic(i1,i3),i1=0,num_layers_short_atomic(i3))
            write(ounit,*)'-------------------------------------------------'
            itemp=0
            do i1=0,num_layers_short_atomic(i3)
              itemp=max(itemp,nodes_short_atomic(i1,i3))
            enddo ! i1
            do i1=1,itemp ! loop over all lines with hidden nodes
              if(i1.le.nodes_short_atomic(0,i3))then ! still input node to be printed
                if(i1.le.maxnodes_short_atomic)then ! still hidden nodes present
                  write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_short_atomic(i1,i2,i3),i2=1,num_layers_short_atomic(i3))
                else
                  write(ounit,'(i4,x,a3)')i1,'  G'
                endif
              else ! no input node in front of hidden nodes
                write(ounit,'(i4,4x,8a3)')i1,(actfunc_short_atomic(i1,i2,i3),i2=1,num_layers_short_atomic(i3))
              endif
            enddo
          enddo ! i3
        endif
!!
        if(lshort.and.(nn_type_short.eq.2).and.(mode.ne.1))then
          do i3=1,npairs
            write(ounit,*)'-------------------------------------------------'
            write(ounit,'(a33,2a3)')' Pair short range NN for element: ',&
              element(elementindex(elempair(i3,1))),element(elementindex(elempair(i3,2)))
            write(ounit,'(a,10i5)')' architecture    ',(nodes_short_pair(i1,i3),i1=0,num_layers_short_pair(i3))
            write(ounit,*)'-------------------------------------------------'
            itemp=0
            do i1=0,num_layers_short_pair(i3)
              itemp=max(itemp,nodes_short_pair(i1,i3))
            enddo ! i1
            do i1=1,itemp ! loop over all lines with hidden nodes
              if(i1.le.nodes_short_pair(0,i3))then ! still input node to be printed
                if(i1.le.maxnodes_short_pair)then ! still hidden nodes present
                  write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_short_pair(i1,i2,i3),i2=1,num_layers_short_pair(i3))
                else
                  write(ounit,'(i4,x,a3)')i1,'  G'
                endif
              else ! no input node in front of hidden nodes
                write(ounit,'(i4,4x,8a3)')i1,(actfunc_short_pair(i1,i2,i3),i2=1,num_layers_short_pair(i3))
              endif
            enddo
          enddo ! i3
        endif
!!
        if(lnntb)then
          if(nntb_flag(0).and.(mode.ne.1)) then
            do i3=1,npairs
              write(ounit,*)'-------------------------------------------------'
              write(ounit,'(a33,2a3)')' Hamiltonian NN for element: ',&
                element(elementindex(elempair(i3,1))),&
                element(elementindex(elempair(i3,2)))
              write(ounit,'(a,10i5)')' architecture    ',&
                (nodes_ham(i1,i3),i1=0,num_layers_ham(i3))
              write(ounit,*)'-------------------------------------------------'
              itemp=0
              do i1=0,num_layers_ham(i3)
                itemp=max(itemp,nodes_ham(i1,i3))
              enddo ! i1
              do i1=1,itemp ! loop over all lines with hidden nodes
                if(i1.le.nodes_ham(0,i3))then ! still input node to be printed
                  if(i1.le.maxnodes_ham)then ! still hidden nodes present
                    write(ounit,'(i4,x,9a3)')i1,'  G',&
                      (actfunc_ham(i1,i2,i3),i2=1,num_layers_ham(i3))
                  else
                    write(ounit,'(i4,x,a3)')i1,'  G'
                  endif
                else ! no input node in front of hidden nodes
                  write(ounit,'(i4,4x,8a3)')i1,&
                    (actfunc_ham(i1,i2,i3),i2=1,num_layers_ham(i3))
                endif
              enddo
            enddo ! i3
          endif
          if(nntb_flag(1).and.(mode.ne.1)) then
          do i3=1,npairs
            write(ounit,*)'-------------------------------------------------'
            write(ounit,'(a33,2a3)')' Hamiltonian Overlap NN for element: ',&
              element(elementindex(elempair(i3,1))),element(elementindex(elempair(i3,2)))
            write(ounit,'(a,10i5)')' architecture    ',&
              (nodes_s(i1,i3),i1=0,num_layers_s(i3))
            write(ounit,*)'-------------------------------------------------'
            itemp=0
            do i1=0,num_layers_s(i3)
              itemp=max(itemp,nodes_s(i1,i3))
            enddo ! i1
            do i1=1,itemp ! loop over all lines with hidden nodes
              if(i1.le.nodes_s(0,i3))then ! still input node to be printed
                if(i1.le.maxnodes_s)then ! still hidden nodes present
                  write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_s(i1,i2,i3),i2=1,num_layers_s(i3))
                else
                  write(ounit,'(i4,x,a3)')i1,'  G'
                endif
              else ! no input node in front of hidden nodes
                write(ounit,'(i4,4x,8a3)')i1,(actfunc_s(i1,i2,i3),i2=1,num_layers_s(i3))
              endif
            enddo
          enddo ! i3
          endif

          if(nntb_flag(2).and.(mode.ne.1)) then
          do i3=1,npairs
            write(ounit,*)'-------------------------------------------------'
            write(ounit,'(a33,2a3)')' Hamiltonian Hexton NN for element: ',&
              element(elementindex(elempair(i3,1))),element(elementindex(elempair(i3,2)))
            write(ounit,'(a,10i5)')' architecture    ',(nodes_hexton(i1,i3),i1=0,num_layers_hexton(i3))
            write(ounit,*)'-------------------------------------------------'
            itemp=0
            do i1=0,num_layers_hexton(i3)
              itemp=max(itemp,nodes_hexton(i1,i3))
            enddo ! i1
            do i1=1,itemp ! loop over all lines with hidden nodes
              if(i1.le.nodes_hexton(0,i3))then ! still input node to be printed
                if(i1.le.maxnodes_hexton)then ! still hidden nodes present
                  write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_hexton(i1,i2,i3),i2=1,num_layers_hexton(i3))
                else
                  write(ounit,'(i4,x,a3)')i1,'  G'
                endif
              else ! no input node in front of hidden nodes
                write(ounit,'(i4,4x,8a3)')i1,(actfunc_hexton(i1,i2,i3),i2=1,num_layers_hexton(i3))
              endif
            enddo
          enddo ! i3
          endif

!! print NN architecture for hextoff NNs
          if(nntb_flag(3).and.(mode.ne.1)) then
            lprint=.false.
            if(mode.eq.3)lprint=.true.
            if(mode.eq.3)then
              i0=ntriplets
            else
              i0 = 1
            endif
            do i3=1,i0
              if((mode.lt.3).and.(((elemtriplet(i3,1).eq.hextoff_training_triplet(1))&
                              .and.(elemtriplet(i3,2).eq.hextoff_training_triplet(2))&
                              .and.(elemtriplet(i3,3).eq.hextoff_training_triplet(3)))&
                              .or.&
                                  ((elemtriplet(i3,1).eq.hextoff_training_triplet(2))&
                              .and.(elemtriplet(i3,2).eq.hextoff_training_triplet(1))&
                              .and.(elemtriplet(i3,3).eq.hextoff_training_triplet(3)))))then
                lprint=.true.
              endif
              if(lprint)then
                write(ounit,*)'-------------------------------------------------'
                write(ounit,'(a38,3a3)')' Hamiltonian Hextoff NN for triplet: ',&
                  element(elementindex(hextoff_training_triplet(1))),&
                  element(elementindex(hextoff_training_triplet(2))),&
                  element(elementindex(hextoff_training_triplet(3)))
!! check if symmetry functions have been defined 
                if(nodes_hextoff(0,i3).eq.0)then
                  write(ounit,*)'ERROR: no symmetry functions have been specified for this atom triplet'
                  stop
                endif
                write(ounit,'(a,10i5)')' architecture    ',(nodes_hextoff(i1,i3),i1=0,num_layers_hextoff(i3))
                write(ounit,*)'-------------------------------------------------'
                itemp=0
                do i1=0,num_layers_hextoff(i3)
                  itemp=max(itemp,nodes_hextoff(i1,i3))
                enddo ! i1
                do i1=1,itemp ! loop over all lines with hidden nodes
                  if(i1.le.nodes_hextoff(0,i3))then ! still input node to be printed
                    if(i1.le.maxnodes_hextoff)then ! still hidden nodes present
                      write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_hextoff(i1,i2,i3),i2=1,num_layers_hextoff(i3))
                    else
                      write(ounit,'(i4,x,a3)')i1,'  G'
                    endif
                  else ! no input node in front of hidden nodes
                    write(ounit,'(i4,4x,8a3)')i1,(actfunc_hextoff(i1,i2,i3),i2=1,num_layers_hextoff(i3))
                  endif
                enddo
              endif
              if(mode.lt.3)then
                lprint=.false.
              endif
            enddo ! i3
          endif
 
          if(nntb_flag(4).and.(mode.ne.1)) then
          do i3=1,npairs
            write(ounit,*)'-------------------------------------------------'
            write(ounit,'(a33,2a3)')' Hamiltonian Density NN for element: ',&
              element(elementindex(elempair(i3,1))),element(elementindex(elempair(i3,2)))
            write(ounit,'(a,10i5)')' architecture    ',(nodes_dens(i1,i3),i1=0,num_layers_dens(i3))
            write(ounit,*)'-------------------------------------------------'
            itemp=0
            do i1=0,num_layers_dens(i3)
              itemp=max(itemp,nodes_dens(i1,i3))
            enddo ! i1
            do i1=1,itemp ! loop over all lines with hidden nodes
              if(i1.le.nodes_dens(0,i3))then ! still input node to be printed
                if(i1.le.maxnodes_dens)then ! still hidden nodes present
                  write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_dens(i1,i2,i3),i2=1,num_layers_dens(i3))
                else
                  write(ounit,'(i4,x,a3)')i1,'  G'
                endif
              else ! no input node in front of hidden nodes
                write(ounit,'(i4,4x,8a3)')i1,(actfunc_dens(i1,i2,i3),i2=1,num_layers_dens(i3))
              endif
            enddo
          enddo ! i3
          endif
        endif
!!
        if(lelec.and.(nn_type_elec.eq.1).and.(mode.ne.1))then
          do i3=1,nelem
            write(ounit,*)'---------------------------------------------------'
            write(ounit,*)'Electrostatic NN for element: ',element(i3)
            write(ounit,'(a,10i5)')' architecture    ',(nodes_elec(i1,i3),i1=0,num_layers_elec(i3))
            write(ounit,*)'---------------------------------------------------'
            itemp=0
            do i1=0,num_layers_elec(i3)
              itemp=max(itemp,nodes_elec(i1,i3))
            enddo ! i1
            do i1=1,itemp ! loop over all lines with hidden nodes
              if(i1.le.nodes_elec(0,i3))then ! still input node to be printed
                if(i1.le.maxnodes_elec)then ! still hidden nodes present
                  write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_elec(i1,i2,i3),i2=1,num_layers_elec(i3))
                else
                  write(ounit,'(i4,x,a3)')i1,'  G'
                endif
              else ! no input node in front of hidden nodes
                write(ounit,'(i4,4x,8a3)')i1,(actfunc_elec(i1,i2,i3),i2=1,num_layers_elec(i3))
              endif
            enddo
          enddo ! i3
        endif
      write(ounit,*)'-------------------------------------------------------------'
!!
!!'
!! write debugging print options
      write(debugunit,*)pstring
      if(pstring(1:1).eq.'1')then
        write(debugunit,*)'printing all short range weights'
      endif
      if(pstring(2:2).eq.'1')then
        write(debugunit,*)'printing all electrostatic weights'
      endif
      if(pstring(3:3).eq.'1')then
        write(debugunit,*)'printing all deshortdw derivatives'
      endif
      if(pstring(4:4).eq.'1')then
        write(debugunit,*)'printing all dfshortdw derivatives'
      endif
      if(pstring(5:5).eq.'1')then
        write(debugunit,*)'printing all ham weights'
      endif
      if(pstring(6:6).eq.'1')then
        write(debugunit,*)'printing all overlap weights'
      endif
      if(pstring(7:7).eq.'1')then
        write(debugunit,*)'printing all hexton range weights'
      endif
      if(pstring(8:8).eq.'1')then
        write(debugunit,*)'printing all hextoff range weights'
      endif
      if(pstring(9:9).eq.'1')then
        write(debugunit,*)'printing all density range weights'
      endif

      
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! sort and analyze symmetry functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! sort short range atomic symmetry functions
      if((nn_type_short.eq.1).and.lshort)then
        call sortsymfunctions(&
          maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
          function_type_short_atomic,symelement_short_atomic,&
          eta_short_atomic,zeta_short_atomic,rshift_short_atomic,&
          lambda_short_atomic,funccutoff_short_atomic)
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
      if(lnntb)then
        if(nntb_flag(0)) then
        call sortpairsymfunctions(npairs,&
          maxnum_funcvalues_ham,num_funcvalues_ham,&
          function_type_ham,symelement_ham,&
          eta_ham,zeta_ham,rshift_ham,lambda_ham,funccutoff_ham,&
          elempair)
        endif
        if(nntb_flag(1)) then
        call sortpairsymfunctions(npairs,&
          maxnum_funcvalues_s,num_funcvalues_s,&
          function_type_s,symelement_s,&
          eta_s,zeta_s,rshift_s,lambda_s,funccutoff_s,&
          elempair)
        endif
        if(nntb_flag(2)) then
        call sortpairsymfunctions(npairs,&
          maxnum_funcvalues_hexton,num_funcvalues_hexton,&
          function_type_hexton,symelement_hexton,&
          eta_hexton,zeta_hexton,rshift_hexton,lambda_hexton,funccutoff_hexton,&
          elempair)
        endif
        if(nntb_flag(3)) then
          if(mode.lt.3)then
            tripletid=tripletindex(hextoff_training_triplet(1),&
              hextoff_training_triplet(2),hextoff_training_triplet(3))
            call sorthextoffsymfunctions(1,ntriplets,tripletid,mode,&
              maxnum_funcvalues_hextoff,num_funcvalues_hextoff,&
              function_type_hextoff,symelement_hextoff,&
              eta_hextoff,zeta_hextoff,rshift_hextoff,lambda_hextoff,funccutoff_hextoff,&
              elemtriplet)
          else
            tripletid=0 ! dummy
            call sorthextoffsymfunctions(ntriplets,ntriplets,tripletid,mode,&
              maxnum_funcvalues_hextoff,num_funcvalues_hextoff,&
              function_type_hextoff,symelement_hextoff,&
              eta_hextoff,zeta_hextoff,rshift_hextoff,lambda_hextoff,funccutoff_hextoff,&
              elemtriplet)
          endif
        endif
        if(nntb_flag(4)) then
        call sortpairsymfunctions(npairs,&
          maxnum_funcvalues_dens,num_funcvalues_dens,&
          function_type_dens,symelement_dens,&
          eta_dens,zeta_dens,rshift_dens,lambda_dens,funccutoff_dens,&
          elempair)
        endif

      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write symmetry functions to output 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,*)' short range atomic symmetry &
                          &functions element ',element(i1),' :'
          write(ounit,*)'-------------------------------------------------------------'
          do i2=1,num_funcvalues_short_atomic(i1)
            if(function_type_short_atomic(i2,i1).eq.1)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.2)then
              write(ounit,'(i5,a3,i3,x,a3,3x,8x,3f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                eta_short_atomic(i2,i1),rshift_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.3)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                element(elementindex(symelement_short_atomic(i2,2,i1))),&
                eta_short_atomic(i2,i1),lambda_short_atomic(i2,i1),&
                zeta_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.4)then
              write(ounit,'(i5,a3,i3,x,a3,3x,16x,2f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                eta_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.5)then
              write(ounit,'(i5,a3,i3,4x,27x,f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),eta_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.6)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.8)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                element(elementindex(symelement_short_atomic(i2,2,i1))),&
                eta_short_atomic(i2,i1),rshift_short_atomic(i2,i1),&
                funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.9)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.4)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                element(elementindex(symelement_short_atomic(i2,2,i1))),&
                eta_short_atomic(i2,i1),lambda_short_atomic(i2,i1),&
                zeta_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            else
              write(ounit,*)'Error: printing unknown symfunction in readinput '
              stop
            endif
          enddo ! i2
        enddo ! i1=1,nelem
      endif ! lshort
!!
      if(lshort.and.(nn_type_short.eq.2))then
        do i1=1,npairs
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a47,2a3,a2)')' short range pair symmetry functions elements ',&
            element(elementindex(elempair(i1,1))),&
            &element(elementindex(elempair(i1,2))),' :'
          write(ounit,*)'-------------------------------------------------------------'
          do i2=1,num_funcvalues_short_pair(i1)
            if(function_type_short_pair(i2,i1).eq.1)then
              write(ounit,'(i5,2a3,i3,x,a3,3x,24x,f8.4)')&
                i2,element(elementindex(elempair(i1,1))),&
                element(elementindex(elempair(i1,2))),&
                function_type_short_pair(i2,i1),&
                element(elementindex(symelement_short_pair(i2,1,i1))),&
                funccutoff_short_pair(i2,i1)
            elseif(function_type_short_pair(i2,i1).eq.2)then
              write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                i2,element(elementindex(elempair(i1,1))),&
                element(elementindex(elempair(i1,2))),&
                function_type_short_pair(i2,i1),&
                funccutoff_short_pair(i2,i1)
            elseif(function_type_short_pair(i2,i1).eq.3)then  ! introduce by Jovan
              write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,17x,f8.3)')&
                i2,element(elementindex(elempair(i1,1))),&
                element(elementindex(elempair(i1,2))),&
                function_type_short_pair(i2,i1),&
                element(elementindex(symelement_short_pair(i2,1,i1))),&
                eta_short_pair(i2,i1),rshift_short_pair(i2,i1),funccutoff_short_pair(i2,i1)
            elseif(function_type_short_pair(i2,i1).eq.4)then  ! introduce by Jovan
              write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,4x,f8.3,8x,f8.3)')&
                i2,element(elementindex(elempair(i1,1))),&
                element(elementindex(elempair(i1,2))),&
                function_type_short_pair(i2,i1),&
                element(elementindex(symelement_short_pair(i2,1,i1))),&
                eta_short_pair(i2,i1),lambda_short_pair(i2,i1),&
                zeta_short_pair(i2,i1),funccutoff_short_pair(i2,i1)
            elseif(function_type_short_pair(i2,i1).eq.5)then  ! introduced by Jovan
              write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                i2,element(elementindex(elempair(i1,1))),&
                element(elementindex(elempair(i1,2))),&
                function_type_short_pair(i2,i1),&
                funccutoff_short_pair(i2,i1)
            elseif(function_type_short_pair(i2,i1).eq.6)then  ! introduced by Jovan
              write(ounit,'(i5,2a3,i3,4x,f8.4,3x,f8.3,17x,f8.3)')&
                i2,element(elementindex(elempair(i1,1))),&
                element(elementindex(elempair(i1,2))),&
                function_type_short_pair(i2,i1),eta_short_pair(i2,i1),rshift_short_pair(i2,i1),&
                funccutoff_short_pair(i2,i1)
            endif
          enddo ! i2
        enddo ! i1=1,npairs
      endif ! lshort
!!
      if(lnntb)then
        if(nntb_flag(0))then
          do i1=1,npairs
            write(ounit,*)'-------------------------------------------------------------'
            write(ounit,'(a47,2a3,a2)')' Hamiltonian symmetry functions elements ',&
              element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),' :'
            write(ounit,*)'-------------------------------------------------------------'
            do i2=1,num_funcvalues_ham(i1)
              if(function_type_ham(i2,i1).eq.1)then
                write(ounit,'(i5,2a3,i3,x,a3,3x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  element(elementindex(symelement_short_pair(i2,1,i1))),&
                  funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.2)then
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.3)then  ! introduce by Jovan
                write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  element(elementindex(symelement_ham(i2,1,i1))),&
                  eta_ham(i2,i1),rshift_ham(i2,i1),funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.4)then  ! introduce by Jovan
                write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,4x,f8.3,8x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  element(elementindex(symelement_ham(i2,1,i1))),&
                  eta_ham(i2,i1),lambda_ham(i2,i1),&
                  zeta_ham(i2,i1),funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.5)then  ! introduced by Jovan
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.6)then  ! introduced by Jovan
                write(ounit,'(i5,2a3,i3,4x,f8.4,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),eta_ham(i2,i1),rshift_ham(i2,i1),&
                  funccutoff_ham(i2,i1)
              endif
            enddo ! i2
          enddo ! i1=1,npairs
        endif

        if(nntb_flag(1))then
          do i1=1,npairs
            write(ounit,*)'-------------------------------------------------------------'
            write(ounit,'(a47,2a3,a2)')' Hamiltonian Overlap symmetry functions elements ',&
              element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),' :'
            write(ounit,*)'-------------------------------------------------------------'
            do i2=1,num_funcvalues_s(i1)
              if(function_type_s(i2,i1).eq.1)then
                write(ounit,'(i5,2a3,i3,x,a3,3x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_s(i2,i1),&
                  element(elementindex(symelement_short_pair(i2,1,i1))),&
                  funccutoff_s(i2,i1)
              elseif(function_type_s(i2,i1).eq.2)then
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_s(i2,i1),&
                  funccutoff_s(i2,i1)
              elseif(function_type_s(i2,i1).eq.3)then  
                write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_s(i2,i1),&
                  element(elementindex(symelement_s(i2,1,i1))),&
                  eta_s(i2,i1),rshift_s(i2,i1),funccutoff_s(i2,i1)
              elseif(function_type_s(i2,i1).eq.4)then  
                write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,4x,f8.3,8x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_s(i2,i1),&
                  element(elementindex(symelement_s(i2,1,i1))),&
                  eta_s(i2,i1),lambda_s(i2,i1),&
                  zeta_s(i2,i1),funccutoff_s(i2,i1)
              elseif(function_type_s(i2,i1).eq.5)then
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_s(i2,i1),&
                  funccutoff_s(i2,i1)
              elseif(function_type_s(i2,i1).eq.6)then
                write(ounit,'(i5,2a3,i3,4x,f8.4,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_s(i2,i1),eta_s(i2,i1),rshift_s(i2,i1),&
                  funccutoff_s(i2,i1)
              endif
            enddo ! i2
          enddo ! i1=1,npairs
        endif
        if(nntb_flag(2))then
          do i1=1,npairs
            write(ounit,*)'-------------------------------------------------------------'
            write(ounit,'(a47,2a3,a2)')' FIX THIS SECTION Hamiltonian Hexton symmetry functions elements ',&
              element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),' :'
            write(ounit,*)'-------------------------------------------------------------'
            do i2=1,num_funcvalues_hexton(i1)
              if(function_type_hexton(i2,i1).eq.1)then
                write(ounit,'(i5,2a3,i3,x,a3,3x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_hexton(i2,i1),&
                  element(elementindex(symelement_short_pair(i2,1,i1))),&
                  funccutoff_hexton(i2,i1)
              elseif(function_type_hexton(i2,i1).eq.2)then
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_hexton(i2,i1),&
                  funccutoff_hexton(i2,i1)
              elseif(function_type_hexton(i2,i1).eq.3)then
                write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_hexton(i2,i1),&
                  element(elementindex(symelement_hexton(i2,1,i1))),&
                  eta_hexton(i2,i1),rshift_hexton(i2,i1),funccutoff_hexton(i2,i1)
              elseif(function_type_hexton(i2,i1).eq.4)then
                write(ounit,'(i5,2a3,i3,x,a3,f8.4,3x,f8.3,4x,f8.3,8x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_hexton(i2,i1),&
                  element(elementindex(symelement_hexton(i2,1,i1))),&
                  eta_hexton(i2,i1),lambda_hexton(i2,i1),&
                  zeta_hexton(i2,i1),funccutoff_hexton(i2,i1)
              elseif(function_type_hexton(i2,i1).eq.5)then
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.4)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_hexton(i2,i1),&
                  funccutoff_hexton(i2,i1)
              elseif(function_type_hexton(i2,i1).eq.6)then
                write(ounit,'(i5,2a3,i3,4x,f8.4,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_hexton(i2,i1),eta_hexton(i2,i1),rshift_hexton(i2,i1),&
                  funccutoff_hexton(i2,i1)
              endif
            enddo ! i2
          enddo ! i1=1,npairs
        endif
        if(nntb_flag(3))then
          lprint=.false.
          if(mode.eq.3)lprint=.true.
          if(mode.eq.3)then
            ndim=ntriplets
          else ! mode 1 or 2
            ndim=1
          endif
          do i3=1,ndim
            if((mode.lt.3).and.(((elemtriplet(i3,1).eq.hextoff_training_triplet(1))&
                                .and.(elemtriplet(i3,2).eq.hextoff_training_triplet(2))&
                                .and.(elemtriplet(i3,3).eq.hextoff_training_triplet(3)))&
                                .or.&
                               ((elemtriplet(i3,1).eq.hextoff_training_triplet(2))&
                                .and.(elemtriplet(i3,2).eq.hextoff_training_triplet(1))&
                                .and.(elemtriplet(i3,3).eq.hextoff_training_triplet(3)))))then
              lprint=.true.
            endif
            if(lprint)then
              write(ounit,*)'-------------------------------------------------------------'
              if(mode.lt.3)then
                write(ounit,'(a47,3a3,a2)') &
                ' Hamiltonian Hextoff symmetry functions elements ',&
                element(elementindex(hextoff_training_triplet(1))),&
                element(elementindex(hextoff_training_triplet(2))),&
                element(elementindex(hextoff_training_triplet(3))),' :'
              else
                write(ounit,'(a47,3a3,a2)') &
                ' Hamiltonian Hextoff symmetry functions elements ',&
                element(elementindex(elemtriplet(i3,1))),&
                element(elementindex(elemtriplet(i3,2))),&
                element(elementindex(elemtriplet(i3,3)))
              endif
              write(ounit,*)'-------------------------------------------------------------'
              i1=i3
              do i2=1,num_funcvalues_hextoff(i3)
                if(function_type_hextoff(i2,i3).eq.1)then
                  write(ounit,'(i5,3a3,i3,x,f8.3,f8.3)')&
                  i2,element(elementindex(hextoff_training_triplet(1))),&
                  element(elementindex(hextoff_training_triplet(2))),&
                  element(elementindex(hextoff_training_triplet(3))),&
                  function_type_hextoff(i2,i3),eta_hextoff(i2,i3),&
                  funccutoff_hextoff(i2,i3)
                elseif(function_type_hextoff(i2,i3).eq.2)then
                  write(ounit,'(i5,3a3,i3,x,f8.3,f8.3)')&
                  i2,element(elementindex(hextoff_training_triplet(1))),&
                  element(elementindex(hextoff_training_triplet(2))),&
                  element(elementindex(hextoff_training_triplet(3))),&
                  function_type_hextoff(i2,i3),&
                  eta_hextoff(i2,i3),&
                  funccutoff_hextoff(i2,i3)
                elseif(function_type_hextoff(i2,i3).eq.3)then
                  write(ounit,'(i5,3a3,i3,x,f8.3,f8.3)')&
                  i2,element(elementindex(hextoff_training_triplet(1))),&
                  element(elementindex(hextoff_training_triplet(2))),&
                  element(elementindex(hextoff_training_triplet(3))),&
                  function_type_hextoff(i2,i3),&
                  eta_hextoff(i2,i3),&
                  funccutoff_hextoff(i2,i3)
                elseif(function_type_ham(i2,i3).eq.4)then  !unused
                  write(ounit,'(i5,2a3,i3,x,a3,f8.3,3x,f8.3,4x,f8.3,8x,f8.3)')&
                  i2,element(elementindex(elempair(i3,1))),&
                  element(elementindex(elempair(i3,2))),&
                  function_type_ham(i2,i3),&
                  element(elementindex(symelement_ham(i2,1,i3))),&
                  eta_ham(i2,i1),lambda_ham(i2,i3),&
                  zeta_ham(i2,i1),funccutoff_ham(i2,i3)
                elseif(function_type_ham(i2,i3).eq.5)then  !unused
                  write(ounit,'(i5,2a3,i3,x,6x,24x,f8.3)')&
                  i2,element(elementindex(elempair(i3,1))),&
                  element(elementindex(elempair(i3,2))),&
                  function_type_ham(i2,i3),&
                  funccutoff_ham(i2,i3)
                elseif(function_type_ham(i2,i3).eq.6)then  !unused
                  write(ounit,'(i5,2a3,i3,4x,f8.3,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i3,1))),&
                  element(elementindex(elempair(i3,2))),&
                  function_type_ham(i2,i3),eta_ham(i2,i3),rshift_ham(i2,i3),&
                  funccutoff_ham(i2,i3)
                endif
              enddo ! i2
            endif
            if(mode.lt.3)then
              lprint=.false.
            endif
          enddo ! i3
        endif
        if(nntb_flag(4))then
          do i1=1,npairs
            write(ounit,*)'-------------------------------------------------------------'
            write(ounit,'(a47,2a3,a2)')' ERROR FIX IT Hamiltonian Density symmetry functions elements ',&
              element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),' :'
            write(ounit,*)'-------------------------------------------------------------'
            do i2=1,num_funcvalues_ham(i1)
              if(function_type_ham(i2,i1).eq.1)then
                write(ounit,'(i5,2a3,i3,x,a3,3x,24x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  element(elementindex(symelement_short_pair(i2,1,i1))),&
                  funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.2)then
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.3)then  ! introduce by Jovan
                write(ounit,'(i5,2a3,i3,x,a3,f8.3,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  element(elementindex(symelement_ham(i2,1,i1))),&
                  eta_ham(i2,i1),rshift_ham(i2,i1),funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.4)then  ! introduce by Jovan
                write(ounit,'(i5,2a3,i3,x,a3,f8.3,3x,f8.3,4x,f8.3,8x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  element(elementindex(symelement_ham(i2,1,i1))),&
                  eta_ham(i2,i1),lambda_ham(i2,i1),&
                  zeta_ham(i2,i1),funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.5)then  ! introduced by Jovan
                write(ounit,'(i5,2a3,i3,x,6x,24x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),&
                  funccutoff_ham(i2,i1)
              elseif(function_type_ham(i2,i1).eq.6)then  ! introduced by Jovan
                write(ounit,'(i5,2a3,i3,4x,f8.3,3x,f8.3,17x,f8.3)')&
                  i2,element(elementindex(elempair(i1,1))),&
                  element(elementindex(elempair(i1,2))),&
                  function_type_ham(i2,i1),eta_ham(i2,i1),rshift_ham(i2,i1),&
                  funccutoff_ham(i2,i1)
              endif
            enddo ! i2
          enddo ! i1=1,npairs
        endif
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,*)' electrostatic symmetry functions element ',element(i1),' :'
          write(ounit,*)'-------------------------------------------------------------'
          do i2=1,num_funcvalues_elec(i1)
            if(function_type_elec(i2,i1).eq.1)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.2)then
              write(ounit,'(i5,a3,i3,x,a3,3x,8x,3f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                eta_elec(i2,i1),rshift_elec(i2,i1),funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.3)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                element(elementindex(symelement_elec(i2,2,i1))),&
                eta_elec(i2,i1),lambda_elec(i2,i1),&
                zeta_elec(i2,i1),funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.4)then
              write(ounit,'(i5,a3,i3,x,a3,3x,16x,2f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                eta_elec(i2,i1),funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.5)then
              write(ounit,'(i5,a3,i3,4x,27x,f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),eta_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.6)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.8)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                element(elementindex(symelement_elec(i2,2,i1))),&
                eta_elec(i2,i1),rshift_elec(i2,i1),&
                funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.9)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                element(elementindex(symelement_elec(i2,2,i1))),&
                eta_elec(i2,i1),lambda_elec(i2,i1),&
                zeta_elec(i2,i1),funccutoff_elec(i2,i1)
            else
              write(ounit,*)'Error: printing unknown symfunctione in readinput '
              stop
            endif
          enddo ! i2
        enddo ! i1=1,nelem
      endif ! lelec 
      write(ounit,*)'-------------------------------------------------------------'
!!'
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
