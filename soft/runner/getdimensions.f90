!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - initialization.f90
!!
      subroutine getdimensions()
!!
      use fileunits
      use nnflags
      use globaloptions
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use nnham
!!
      implicit none
!!
      integer i,j,k,i1,i2,i3              ! internal
      integer ztemp                ! internal
      integer ztemp1               ! internal
      integer ztemp2               ! internal
      integer ztemp3               ! internal
      integer function_type_local  ! internal
      integer, dimension(:)  , allocatable :: nodes_short_local     ! internal
      integer, dimension(:)  , allocatable :: nodes_ewald_local     ! internal
      integer, dimension(:)  , allocatable :: nodes_pair_local      ! internal
      integer, dimension(:)  , allocatable :: nodes_s_local      ! internal
      integer, dimension(:)  , allocatable :: nodes_hexton_local      ! internal
      integer, dimension(:)  , allocatable :: nodes_hextoff_local      ! internal
      integer, dimension(:)  , allocatable :: nodes_dens_local      ! internal
      integer, dimension(:)  , allocatable :: nodes_ham_local      ! internal
      integer, dimension(:)  , allocatable :: num_funcvalues_local  ! internal
      integer, dimension(:)  , allocatable :: num_funcvaluese_local ! internal
      integer, dimension(:,:), allocatable :: num_funcvaluesp_local ! internal
      integer, dimension(:,:), allocatable :: num_funcvalues_s_local ! internal
      integer, dimension(:,:), allocatable :: num_funcvalues_hexton_local ! internal
      integer, dimension(:,:,:), allocatable :: num_funcvalues_hextoff_local ! internal
      integer, dimension(:,:), allocatable :: num_funcvalues_dens_local ! internal
      integer, dimension(:,:), allocatable :: num_funcvalues_ham_local ! internal

!!
      character*40 dummy       ! internal
      character*40 keyword     ! internal
      character*2 elementtemp  ! internal
      character*2 elementtemp1 ! internal
      character*2 elementtemp2 ! internal
      character*2 elementtemp3 ! internal
!!
      logical lfounddebug
      logical lfound_num_layersshort
      logical lfound_num_layersewald
      logical lfound_num_layerspair
      logical lfound_num_layers_s
      logical lfound_num_layers_hexton
      logical lfound_num_layers_hextoff
      logical lfound_num_layers_dens
      logical lfound_num_layers_ham
      logical lfound_luseatomenergies
      logical lfound_luseatomcharges
      logical lfound_nelem
!!
      lfounddebug            =.false.
      lfound_num_layersshort =.false.
      lfound_num_layersewald =.false.
      lfound_num_layerspair  =.false.
      lfound_num_layers_ham     =.false.

      lfound_num_layers_s    =.false.
      lfound_num_layers_hexton    =.false.
      lfound_num_layers_hextoff    =.false.
      lfound_num_layers_dens    =.false.
      lfound_luseatomenergies=.false.
      lfound_luseatomcharges =.false.
      lfound_nelem           =.false.
      maxnodes_short_atomic  =0
      maxnodes_elec          =0
      maxnodes_short_pair    =0
      maxnodes_ham           =0
      maxnodes_s             =0
      maxnodes_hexton        =0
      maxnodes_hextoff       =0
      maxnodes_dens          =0
      
      
      maxnum_funcvalues_short_atomic      =0
      maxnum_funcvalues_elec              =0
      maxnum_funcvalues_short_pair        =0
      maxnum_funcvalues_ham               =0
      maxnum_funcvalues_s                 =0
      maxnum_funcvalues_hexton            =0
      maxnum_funcvalues_hextoff           =0
      maxnum_funcvalues_dens              =0
      maxnum_layers_short_atomic          =0
      maxnum_layers_elec                  =0
      maxnum_layers_short_pair            =0
      maxnum_layers_ham                   =0
      maxnum_layers_s                     =0
      maxnum_layers_hexton                =0
      maxnum_layers_hextoff               =0
      maxnum_layers_dens                  =0

      nn_type_short          =1 ! default value consistent with readinput.f90
      nn_type_elec           =1 !
      nn_type_nntb           =0 ! default value consistent with readinput.f90
      nntb_flag              =.false.
      lshort                 =.false.
      lelec                  =.false.
      lnntb                  =.false.
      mode                   =0
!!
!! determine which NNs to expect in input.nn
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 13   continue
      read(nnunit,*,END=14)keyword
      if(keyword.eq.'nn_type_short')then
        backspace(nnunit)
        read(nnunit,*)dummy,nn_type_short
      endif
      if(keyword.eq.'nn_type_nntb')then
        backspace(nnunit)
        read(nnunit,*)dummy,nn_type_nntb
      endif
      if(keyword.eq.'runner_mode')then
        backspace(nnunit)
        read(nnunit,*)dummy,mode
      endif
      if(keyword.eq.'use_short_nn')then
        lshort=.true.
      endif
      if(keyword.eq.'use_electrostatics')then
        lelec=.true.
      endif
      if(keyword.eq.'use_hamiltonian')then
        lnntb=.true.
      endif
      if(keyword.eq.'use_ham')then
        nntb_flag(0)=.true.
      endif
      if(keyword.eq.'use_overlap')then
        nntb_flag(1)=.true.
      endif
      if(keyword.eq.'use_hexton')then
        nntb_flag(2)=.true.
      endif
      if(keyword.eq.'use_hextoff')then
        nntb_flag(3)=.true.
      endif
      if(keyword.eq.'use_dens')then
        nntb_flag(4)=.true.
      endif

      if((keyword.eq.'electrostatic_type').or.(keyword.eq.'nn_type_elec'))then
        backspace(nnunit)
        read(nnunit,*)dummy,nn_type_elec
      endif
      goto 13
 14   continue
      close(nnunit)



       

!!
!! determine maxnum_layers_short_atomic, maxnum_layersewald, maxnum_layerspair, maxnum_layersham 
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 10   continue
      read(nnunit,*,END=20)keyword
      if(keyword.eq.'debug_mode')then
        lfounddebug=.true.
        ldebug=.true.
      elseif(keyword.eq.'global_hidden_layers_short')then
        lfound_num_layersshort=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_short_atomic
        maxnum_layers_short_atomic=maxnum_layers_short_atomic + 1
      elseif(keyword.eq.'global_hidden_layers_electrostatic')then
        lfound_num_layersewald=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_elec
        maxnum_layers_elec=maxnum_layers_elec + 1
      elseif(keyword.eq.'global_hidden_layers_pair')then
        lfound_num_layerspair=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_short_pair
        maxnum_layers_short_pair=maxnum_layers_short_pair + 1
      elseif(keyword.eq.'global_hidden_layers_ham')then
        lfound_num_layers_ham=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_ham
        maxnum_layers_ham=maxnum_layers_ham + 1
      elseif(keyword.eq.'global_hidden_layers_s')then
        lfound_num_layers_s=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_s
        maxnum_layers_s=maxnum_layers_s + 1
      elseif(keyword.eq.'global_hidden_layers_hexton')then
        lfound_num_layers_hexton=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_hexton
        maxnum_layers_hexton=maxnum_layers_hexton + 1
      elseif(keyword.eq.'global_hidden_layers_hextoff')then
        lfound_num_layers_hextoff=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_hextoff
        maxnum_layers_hextoff=maxnum_layers_hextoff + 1
      elseif(keyword.eq.'global_hidden_layers_dens')then
        lfound_num_layers_dens=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,maxnum_layers_dens
        maxnum_layers_dens=maxnum_layers_dens + 1


      elseif(keyword.eq.'use_atom_energies')then
        lfound_luseatomenergies=.true.
        luseatomenergies=.true.
      elseif(keyword.eq.'use_atom_charges')then
        lfound_luseatomcharges=.true.
        luseatomcharges=.true.
      elseif(keyword.eq.'number_of_elements')then
        lfound_nelem=.true.
        backspace(nnunit)
        read(nnunit,*)dummy,nelem
        npairs=0
        do i1=1,nelem
          do i2=i1,nelem
            npairs=npairs+1
          enddo
        enddo
        ntriplets=0
!! loop over unique pairs of elements
        do i1=1,nelem
          do i2=i1,nelem
!! loop over all possible neighboring atoms
            do i3=1,nelem
              ntriplets = ntriplets+1
            enddo
          enddo
        enddo
      endif
      goto 10
 20   continue
      close(nnunit)
!!
!! check if mandatory keywords have been found
      if(lshort.and.(nn_type_short.eq.1).and.(maxnum_layers_short_atomic.eq.0))then
        write(ounit,*)'Error: please specify global_hidden_layers_short'
        stop
      endif
      if(lshort.and.(nn_type_short.eq.2).and.(maxnum_layers_short_pair.eq.0))then
        write(ounit,*)'Error: please specify global_hidden_layers_pair'
        stop
      endif

      if((nn_type_nntb.eq.1).and.(maxnum_layers_ham.eq.0) .and. &
        (maxnum_layers_s.eq.0) .and. (maxnum_layers_hexton.eq.0) .and.&
        (maxnum_layers_hextoff.eq.0) .and. (maxnum_layers_dens.eq.0))then         ! Right now this is the same for all parts - we are unlikely to train all of them at once
        write(ounit,*)'Error: please specify global_hidden_layers_ham/s/hexton/hextoff/dens'
        stop
      endif
      if(lelec.and.(nn_type_elec.eq.1).and.(maxnum_layers_elec.eq.0))then
        write(ounit,*)'Error: please specify global_hidden_layers_electrostatic'
        stop   !'
      endif
!!
      allocate(nucelem(nelem))
      allocate(element(nelem))
!!
!! get nucelem
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 11   continue
      read(nnunit,*,END=12)keyword
      if(keyword.eq.'elements')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(element(i),i=1,nelem)
      endif
      goto 11
 12   continue
      close(nnunit)
      do i=1,nelem
        call nuccharge(element(i),nucelem(i))
      enddo
      call sortelements()
!!
!! check if everything has been found and put defaults
      if(.not.lfounddebug)then
        ldebug=.false.
      endif
      if(.not.lfound_num_layersshort)then
        maxnum_layers_short_atomic=0 
      endif
      if(.not.lfound_num_layersewald)then
        maxnum_layers_elec=0 
      endif
      if(.not.lfound_num_layerspair)then
        maxnum_layers_short_pair=0 
      endif
      if(.not.lfound_num_layers_ham)then
        maxnum_layers_ham=0 
      endif
      if(.not.lfound_num_layers_s)then
        maxnum_layers_s=0
      endif
      if(.not.lfound_num_layers_hexton)then
        maxnum_layers_hexton=0
      endif
      if(.not.lfound_num_layers_hextoff)then
        maxnum_layers_hextoff=0
      endif
      if(.not.lfound_num_layers_dens)then
        maxnum_layers_dens=0
      endif

      if(.not.lfound_luseatomenergies)then
         luseatomenergies=.false. 
      endif
      if(.not.lfound_luseatomcharges)then
         luseatomcharges=.false. 
      endif
      if(.not.lfound_nelem)then
        write(ounit,*)'Error: nelem not found in input.nn'
        stop
      endif
!!
!! determine maxnodes_short_atomic, maxnodes_ewald, maxnum_funcvalues, maxnum_funcvaluese
      allocate(num_funcvalues_local(102))
      allocate(num_funcvaluese_local(102))
      allocate(num_funcvaluesp_local(102,102))
      allocate(num_funcvalues_ham_local(102,102))
      allocate(num_funcvalues_s_local(102,102))
      allocate(num_funcvalues_hexton_local(102,102))
      allocate(num_funcvalues_hextoff_local(102,102,102))
      allocate(num_funcvalues_dens_local(102,102))
      num_funcvalues_local(:)=0
      num_funcvaluese_local(:)=0
      num_funcvaluesp_local(:,:)=0
      num_funcvalues_ham_local(:,:)=0
      num_funcvalues_s_local(:,:)=0
      num_funcvalues_hexton_local(:,:)=0
      num_funcvalues_hextoff_local(:,:,:)=0
      num_funcvalues_dens_local(:,:)=0

!!
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
!! 
      if(maxnum_layers_short_atomic.gt.0)then
        allocate (nodes_short_local(0:maxnum_layers_short_atomic))
        nodes_short_local(:)=0
      endif
      if(maxnum_layers_elec.gt.0)then
        allocate (nodes_ewald_local(0:maxnum_layers_elec))
        nodes_ewald_local(:)=0
      endif
      if(maxnum_layers_short_pair.gt.0)then
        allocate (nodes_pair_local(0:maxnum_layers_short_pair))
        nodes_pair_local(:)=0
      endif
      if(maxnum_layers_ham.gt.0)then
        allocate (nodes_ham_local(0:maxnum_layers_ham))
        nodes_ham_local(:)=0
      endif
      if(maxnum_layers_s.gt.0)then
        allocate (nodes_s_local(0:maxnum_layers_s))
        nodes_s_local(:)=0
      endif
      if(maxnum_layers_hexton.gt.0)then
        allocate (nodes_hexton_local(0:maxnum_layers_hexton))
        nodes_hexton_local(:)=0
      endif
      if(maxnum_layers_hextoff.gt.0)then
        allocate (nodes_hextoff_local(0:maxnum_layers_hextoff))
        nodes_hextoff_local(:)=0
      endif
      if(maxnum_layers_dens.gt.0)then
        allocate (nodes_dens_local(0:maxnum_layers_dens))
        nodes_dens_local(:)=0
      endif

!!
 30   continue
      read(nnunit,*,END=40)keyword
      if(keyword.eq.'global_nodes_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_short_local(i),i=1,maxnum_layers_short_atomic-1)
!!
      elseif(keyword.eq.'global_nodes_electrostatic')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_ewald_local(i),i=1,maxnum_layers_elec-1)
!!
      elseif(keyword.eq.'global_nodes_pair')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_pair_local(i),i=1,maxnum_layers_short_pair-1)
!!
      elseif(keyword.eq.'global_nodes_ham')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_ham_local(i),i=1,maxnum_layers_ham-1)
      elseif(keyword.eq.'global_nodes_s')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_s_local(i),i=1,maxnum_layers_s-1)
      elseif(keyword.eq.'global_nodes_hexton')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_hexton_local(i),i=1,maxnum_layers_hexton-1)
      elseif(keyword.eq.'global_nodes_hextoff')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_hextoff_local(i),i=1,maxnum_layers_hextoff-1)
      elseif(keyword.eq.'global_nodes_dens')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(nodes_dens_local(i),i=1,maxnum_layers_dens-1)

!!
      elseif(keyword.eq.'element_symfunction_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,function_type_local
!!        write(ounit,*)dummy,elementtemp
        call nuccharge(elementtemp,ztemp)
        if((function_type_local.eq.1).or.&
           (function_type_local.eq.2).or.&
           (function_type_local.eq.4))then
          num_funcvalues_local(ztemp)=num_funcvalues_local(ztemp)+nelem 
        elseif((function_type_local.eq.3).or.&
           (function_type_local.eq.8).or.&
           (function_type_local.eq.9))then
          num_funcvalues_local(ztemp)=num_funcvalues_local(ztemp)+nelem ! 2 neighbors of the same element
!! add cross terms
          if(nelem.gt.1)then
            do i1=1,nelem-1
              num_funcvalues_local(ztemp)=num_funcvalues_local(ztemp)+i1    ! 2 neighbors of different element
            enddo
          endif
        elseif(function_type_local.eq.5)then
          num_funcvalues_local(ztemp)=num_funcvalues_local(ztemp)+1 
        elseif(function_type_local.eq.6)then
          num_funcvalues_local(ztemp)=num_funcvalues_local(ztemp)+1 
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions ',function_type_local
          stop !'
        endif
!!
      elseif(keyword.eq.'element_symfunction_electrostatic')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp,function_type_local
!!        write(ounit,*)dummy,elementtemp
        call nuccharge(elementtemp,ztemp)
        if((function_type_local.eq.1).or.&
           (function_type_local.eq.2).or.&
           (function_type_local.eq.4))then
          num_funcvaluese_local(ztemp)=num_funcvaluese_local(ztemp)+nelem 
        elseif((function_type_local.eq.3).or.&
               (function_type_local.eq.8).or.&
               (function_type_local.eq.9))then
          num_funcvaluese_local(ztemp)=num_funcvaluese_local(ztemp)+nelem ! 2 neighbors of the same element
!! add cross terms
          if(nelem.gt.1)then
            do i1=1,nelem-1
              num_funcvaluese_local(ztemp)=num_funcvaluese_local(ztemp)+i1    ! 2 neighbors of different element
            enddo
          endif
        elseif(function_type_local.eq.5)then
          num_funcvaluese_local(ztemp)=num_funcvaluese_local(ztemp)+1 
        elseif(function_type_local.eq.6)then
          num_funcvaluese_local(ztemp)=num_funcvaluese_local(ztemp)+1 
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
!!
      elseif(keyword.eq.'global_symfunction_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_local
!!        write(ounit,*)dummy,function_type_local
        if((function_type_local.eq.1).or.&
          (function_type_local.eq.2).or.&
          (function_type_local.eq.4))then
          do i1=1,nelem
            num_funcvalues_local(nucelem(i1))=num_funcvalues_local(nucelem(i1))+nelem 
          enddo
        elseif((function_type_local.eq.3).or.&
               (function_type_local.eq.8).or.&
               (function_type_local.eq.9))then
          do i1=1,nelem
            num_funcvalues_local(nucelem(i1))=num_funcvalues_local(nucelem(i1))+nelem ! 2 neighbors of the same element
          enddo
          do i1=1,nelem
            if(nelem.gt.1)then
              do i2=1,nelem-1
                num_funcvalues_local(nucelem(i1))=num_funcvalues_local(nucelem(i1))+i2    ! 2 neighbors of different element
              enddo ! i2
            endif
          enddo ! i1
        elseif(function_type_local.eq.5)then
          do i1=1,nelem
            num_funcvalues_local(nucelem(i1))=num_funcvalues_local(nucelem(i1))+1 
          enddo
        elseif(function_type_local.eq.6)then
          do i1=1,nelem
            num_funcvalues_local(nucelem(i1))=num_funcvalues_local(nucelem(i1))+1 
          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
!!
      elseif(keyword.eq.'global_symfunction_electrostatic')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_local
!!        write(ounit,*)dummy
        if((function_type_local.eq.1).or.&
          (function_type_local.eq.2).or.&
          (function_type_local.eq.4))then
          do i1=1,nelem
            num_funcvaluese_local(nucelem(i1))=num_funcvaluese_local(nucelem(i1))+nelem
          enddo
        elseif((function_type_local.eq.3).or.&
               (function_type_local.eq.8).or.&
               (function_type_local.eq.9))then
          do i1=1,nelem
            num_funcvaluese_local(nucelem(i1))=num_funcvaluese_local(nucelem(i1))+nelem ! 2 neighbors of the same element
          enddo
          do i1=1,nelem
            if(nelem.gt.1)then
              do i2=1,nelem-1
                num_funcvaluese_local(nucelem(i1))=num_funcvaluese_local(nucelem(i1))+i2    ! 2 neighbors of different element
              enddo ! i2
            endif
          enddo ! i1
        elseif(function_type_local.eq.5)then
          do i1=1,nelem
            num_funcvaluese_local(nucelem(i1))=num_funcvaluese_local(nucelem(i1))+1
          enddo
        elseif(function_type_local.eq.6)then
          do i1=1,nelem
            num_funcvaluese_local(nucelem(i1))=num_funcvaluese_local(nucelem(i1))+1
          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
!!
      elseif(keyword.eq.'symfunction_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp
!!        write(ounit,*)dummy,elementtemp
        call nuccharge(elementtemp,ztemp)
        num_funcvalues_local(ztemp)=num_funcvalues_local(ztemp)+1 
!!
      elseif(keyword.eq.'symfunction_electrostatic')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp
!!        write(ounit,*)dummy,elementtemp
        call nuccharge(elementtemp,ztemp)
        num_funcvaluese_local(ztemp)=num_funcvaluese_local(ztemp)+1 
!!
      elseif(keyword.eq.'pairsymfunction_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,function_type_local
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        if(function_type_local.eq.1)then
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        elseif(function_type_local.eq.2)then ! just bond length
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        elseif(function_type_local.eq.3)then  ! Introduced By Jovan
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        elseif(function_type_local.eq.4)then  ! Introduced By Jovan
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        elseif(function_type_local.eq.5)then ! just bond length
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        elseif(function_type_local.eq.6)then ! Introduced Jovan exp*fc
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
!!
      elseif(keyword.eq.'element_pairsymfunction_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,function_type_local
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        if(function_type_local.eq.1)then
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+nelem
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+nelem
          endif
        elseif(function_type_local.eq.2)then ! just bond length
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        elseif(function_type_local.eq.3)then     ! introduce by Jovan
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+nelem
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+nelem
          endif
        elseif(function_type_local.eq.4)then     ! introduce by Jovan
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+nelem
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+nelem
          endif
        elseif(function_type_local.eq.5)then ! just bond length
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        elseif(function_type_local.eq.6)then ! Introduced Jovan exp*fc
          if(ztemp1.le.ztemp2)then
            num_funcvaluesp_local(ztemp1,ztemp2)=num_funcvaluesp_local(ztemp1,ztemp2)+1
          else
            num_funcvaluesp_local(ztemp2,ztemp1)=num_funcvaluesp_local(ztemp2,ztemp1)+1
          endif
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
!!
      elseif(keyword.eq.'global_pairsymfunction_short')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_local
        if(function_type_local.eq.1)then
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvaluesp_local(nucelem(i1),nucelem(i2))&
                =num_funcvaluesp_local(nucelem(i1),nucelem(i2))+nelem
            enddo
          enddo
        elseif(function_type_local.eq.2)then ! just bond length
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvaluesp_local(nucelem(i1),nucelem(i2))&
                =num_funcvaluesp_local(nucelem(i1),nucelem(i2))+1
            enddo
          enddo
        elseif(function_type_local.eq.3)then    ! introduced by Jovan
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvaluesp_local(nucelem(i1),nucelem(i2))&
                =num_funcvaluesp_local(nucelem(i1),nucelem(i2))+nelem
            enddo
          enddo
        elseif(function_type_local.eq.4)then    ! introduced by Jovan
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvaluesp_local(nucelem(i1),nucelem(i2))&
                =num_funcvaluesp_local(nucelem(i1),nucelem(i2))+nelem
            enddo
          enddo
        elseif(function_type_local.eq.5)then ! just bond length
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvaluesp_local(nucelem(i1),nucelem(i2))&
                =num_funcvaluesp_local(nucelem(i1),nucelem(i2))+1
            enddo
          enddo
        elseif(function_type_local.eq.6)then ! Introduced by Jovan exp*fc
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvaluesp_local(nucelem(i1),nucelem(i2))&
                =num_funcvaluesp_local(nucelem(i1),nucelem(i2))+1
            enddo
          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
!!
!      elseif(keyword.eq.'hamsymfunction')then
!        backspace(nnunit)
!        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,function_type_local
!        call nuccharge(elementtemp1,ztemp1)
!        call nuccharge(elementtemp2,ztemp2)
!        if(function_type_local.eq.1)then
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        elseif(function_type_local.eq.2)then 
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        elseif(function_type_local.eq.3)then  
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        elseif(function_type_local.eq.4)then  
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        elseif(function_type_local.eq.5)then 
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        elseif(function_type_local.eq.6)then
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        else
!          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
!          stop !'
!        endif
!!
!      elseif(keyword.eq.'element_hamsymfunction')then
!        backspace(nnunit)
!        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,function_type_local
!        call nuccharge(elementtemp1,ztemp1)
!        call nuccharge(elementtemp2,ztemp2)
!        if(function_type_local.eq.1)then
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+nelem
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+nelem
!          endif
!        elseif(function_type_local.eq.2)then 
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        elseif(function_type_local.eq.3)then    
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+nelem
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+nelem
!          endif
!        elseif(function_type_local.eq.4)then     
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+nelem
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+nelem
!          endif
!        elseif(function_type_local.eq.5)then 
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        elseif(function_type_local.eq.6)then
!          if(ztemp1.le.ztemp2)then
!            num_funcvaluesh_local(ztemp1,ztemp2)=num_funcvaluesh_local(ztemp1,ztemp2)+1
!          else
!            num_funcvaluesh_local(ztemp2,ztemp1)=num_funcvaluesh_local(ztemp2,ztemp1)+1
!          endif
!        else
!          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
!          stop !'
!        endif
!!
      elseif(keyword.eq.'symfunction_hextoff')then
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elementtemp1,elementtemp2,elementtemp3,function_type_local
        call nuccharge(elementtemp1,ztemp1)
        call nuccharge(elementtemp2,ztemp2)
        call nuccharge(elementtemp3,ztemp3)
        if(function_type_local.eq.1)then         ! this is the distance only between atom A and B
              num_funcvalues_hextoff_local(ztemp1,ztemp2,ztemp3)&
                =num_funcvalues_hextoff_local(ztemp1,ztemp2,ztemp3)+1
        elseif(function_type_local.eq.2)then ! this is the distance between AC + B
              num_funcvalues_hextoff_local(ztemp1,ztemp2,ztemp3)&
                =num_funcvalues_hextoff_local(ztemp1,ztemp2,ztemp3)+1
        elseif(function_type_local.eq.3)then    ! this is the distance between AC - BC
              num_funcvalues_hextoff_local(ztemp1,ztemp2,ztemp3)&
                =num_funcvalues_hextoff_local(ztemp1,ztemp2,ztemp3)+1
!        elseif(function_type_local.eq.4)then    ! introduced by Jovan
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+nelem
!            enddo
!          enddo
!        elseif(function_type_local.eq.5)then ! just bond length
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
!        elseif(function_type_local.eq.6)then ! Introduced by Jovan exp*fc
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
      elseif(keyword.eq.'global_symfunction_hexton')then
        write(*,*) 'dimensions for symfunction_hexton'
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_local
        if(function_type_local.eq.1)then         ! this is the distance only between atom A and B
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_hexton_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_hexton_local(nucelem(i1),nucelem(i2))+ 1
            enddo
          enddo
        elseif(function_type_local.eq.2)then ! this is the distance between AC + BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_hexton_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_hexton_local(nucelem(i1),nucelem(i2))+ nelem
            enddo
          enddo
        elseif(function_type_local.eq.3)then    ! this is the distance between AC - BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_hexton_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_hexton_local(nucelem(i1),nucelem(i2))+nelem
            enddo
          enddo
!        elseif(function_type_local.eq.4)then    ! introduced by Jovan
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+nelem
!            enddo
!          enddo
!        elseif(function_type_local.eq.5)then ! just bond length
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
!        elseif(function_type_local.eq.6)then ! Introduced by Jovan exp*fc
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
      elseif(keyword.eq.'global_symfunction_s')then
        write(*,*) 'dimensions for symfunction_s'
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_local
        if(function_type_local.eq.1)then         ! this is the distance only between atom A and B
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_s_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_s_local(nucelem(i1),nucelem(i2))+ 1
            enddo
          enddo
        elseif(function_type_local.eq.2)then ! this is the distance between AC + BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_s_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_s_local(nucelem(i1),nucelem(i2))+ nelem
            enddo
          enddo
        elseif(function_type_local.eq.3)then    ! this is the distance between AC - BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_s_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_s_local(nucelem(i1),nucelem(i2))+nelem
            enddo
          enddo
!        elseif(function_type_local.eq.4)then    ! introduced by Jovan
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+nelem
!            enddo
!          enddo
!        elseif(function_type_local.eq.5)then ! just bond length
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
!        elseif(function_type_local.eq.6)then ! Introduced by Jovan exp*fc
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
      elseif(keyword.eq.'global_symfunction_dens')then
        write(*,*) 'dimensions for symfunction_dens'
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_local
        if(function_type_local.eq.1)then         ! this is the distance only between atom A and B
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_dens_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_dens_local(nucelem(i1),nucelem(i2))+ 1
            enddo
          enddo
        elseif(function_type_local.eq.2)then ! this is the distance between AC + BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_dens_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_dens_local(nucelem(i1),nucelem(i2))+ nelem
            enddo
          enddo
        elseif(function_type_local.eq.3)then    ! this is the distance between AC - BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_dens_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_dens_local(nucelem(i1),nucelem(i2))+nelem
            enddo
          enddo
!        elseif(function_type_local.eq.4)then    ! introduced by Jovan
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+nelem
!            enddo
!          enddo
!        elseif(function_type_local.eq.5)then ! just bond length
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
!        elseif(function_type_local.eq.6)then ! Introduced by Jovan exp*fc
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif
      elseif(keyword.eq.'global_symfunction_ham')then
        write(*,*) 'dimensions for symfunction_ham'
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,function_type_local
        if(function_type_local.eq.1)then         ! this is the distance only between atom A and B
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_ham_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_ham_local(nucelem(i1),nucelem(i2))+ 1
            enddo
          enddo
        elseif(function_type_local.eq.2)then ! this is the distance between AC + BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_ham_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_ham_local(nucelem(i1),nucelem(i2))+ nelem
            enddo
          enddo
        elseif(function_type_local.eq.3)then    ! this is the distance between AC - BC
          do i1=1,nelem
            do i2=i1,nelem
              num_funcvalues_ham_local(nucelem(i1),nucelem(i2))&
                =num_funcvalues_ham_local(nucelem(i1),nucelem(i2))+nelem
            enddo
          enddo
!        elseif(function_type_local.eq.4)then    ! introduced by Jovan
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+nelem
!            enddo
!          enddo
!        elseif(function_type_local.eq.5)then ! just bond length
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
!        elseif(function_type_local.eq.6)then ! Introduced by Jovan exp*fc
!          do i1=1,nelem
!            do i2=i1,nelem
!              num_funcvaluesh_local(nucelem(i1),nucelem(i2))&
!                =num_funcvaluesh_local(nucelem(i1),nucelem(i2))+1
!            enddo
!          enddo
        else
          write(ounit,*)'ERROR: symfunction_type not implemented in getdimensions'
          stop !'
        endif


      endif

      goto 30
 40   continue
      close(nnunit)
!!
      if(maxnum_layers_short_atomic.gt.0)then
        do i=1,maxnum_layers_short_atomic ! loop over all hidden and output layers
          maxnodes_short_atomic=max(maxnodes_short_atomic,nodes_short_local(i))
        enddo
      endif
      if(maxnum_layers_elec.gt.0)then
        do i=1,maxnum_layers_elec
          maxnodes_elec=max(maxnodes_elec,nodes_ewald_local(i))
        enddo
      endif
      if(maxnum_layers_short_pair.gt.0)then
        do i=1,maxnum_layers_short_pair ! loop over all hidden and output layers
          maxnodes_short_pair=max(maxnodes_short_pair,nodes_pair_local(i))
        enddo
      endif
      if(maxnum_layers_ham.gt.0)then
        do i=1,maxnum_layers_ham ! loop over all hidden and output layers
          maxnodes_ham=max(maxnodes_ham,nodes_ham_local(i))
        enddo
      endif
      if(maxnum_layers_s.gt.0)then
        do i=1,maxnum_layers_s ! loop over all hidden and output layers
          maxnodes_s=max(maxnodes_s,nodes_s_local(i))
        enddo
      endif
      if(maxnum_layers_hexton.gt.0)then
        do i=1,maxnum_layers_hexton ! loop over all hidden and output layers
          maxnodes_hexton=max(maxnodes_hexton,nodes_hexton_local(i))
        enddo
      endif
      if(maxnum_layers_hextoff.gt.0)then
        do i=1,maxnum_layers_hextoff ! loop over all hidden and output layers          
          maxnodes_hextoff=max(maxnodes_hextoff,nodes_hextoff_local(i))
        enddo
      endif
      if(maxnum_layers_dens.gt.0)then
        do i=1,maxnum_layers_dens ! loop over all hidden and output layers
          maxnodes_dens=max(maxnodes_dens,nodes_dens_local(i))
        enddo
      endif

!!
      if(allocated(nodes_short_local))deallocate(nodes_short_local)
      if(allocated(nodes_ewald_local))deallocate(nodes_ewald_local)
      if(allocated(nodes_pair_local))deallocate(nodes_pair_local)
      if(allocated(nodes_ham_local))deallocate(nodes_ham_local)
      if(allocated(nodes_s_local))deallocate(nodes_s_local)
      if(allocated(nodes_hexton_local))deallocate(nodes_hexton_local)
      if(allocated(nodes_hextoff_local))deallocate(nodes_hextoff_local)
      if(allocated(nodes_dens_local))deallocate(nodes_dens_local)

!!
      do i=1,102
        maxnum_funcvalues_short_atomic =max(maxnum_funcvalues_short_atomic ,num_funcvalues_local(i))
        maxnum_funcvalues_elec=max(maxnum_funcvalues_elec,num_funcvaluese_local(i))
      enddo
      do i=1,102
        do j=i,102
          maxnum_funcvalues_short_pair=max(maxnum_funcvalues_short_pair,num_funcvaluesp_local(i,j))
          maxnum_funcvalues_ham=max(maxnum_funcvalues_ham,num_funcvalues_hexton_local(i,j))       !!! CMH this will need to be turned off eventually
          maxnum_funcvalues_s=max(maxnum_funcvalues_s,num_funcvalues_s_local(i,j))
          maxnum_funcvalues_hexton=max(maxnum_funcvalues_hexton,num_funcvalues_hexton_local(i,j))
!          maxnum_funcvalues_hextoff=max(maxnum_funcvalues_hextoff,num_funcvalues_hextoff_local(i,j))
          maxnum_funcvalues_dens=max(maxnum_funcvalues_dens,num_funcvalues_dens_local(i,j))

        enddo
      enddo
      do i=1,102
        do j=1,102
          do k=1,102
            maxnum_funcvalues_hextoff=max(maxnum_funcvalues_hextoff,num_funcvalues_hextoff_local(i,j,k))
          enddo
        enddo
      enddo
!!
      deallocate(num_funcvalues_local)
      deallocate(num_funcvaluese_local)
      deallocate(num_funcvaluesp_local)
      deallocate(num_funcvalues_s_local)
      deallocate(num_funcvalues_hexton_local)
      deallocate(num_funcvalues_hextoff_local)
      deallocate(num_funcvalues_dens_local)

      deallocate(nucelem)
      deallocate(element)
      
      return
!!
 99   continue
      write(ounit,*)'Error in getdimensions ',keyword
      stop

      end
