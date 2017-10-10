!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

      module symfunctions 
!!
      implicit none
!!
      integer, dimension(:,:), allocatable :: function_type_short_atomic
      integer, dimension(:,:), allocatable :: function_type_elec
      integer, dimension(:,:), allocatable :: function_type_short_pair
      integer, dimension(:,:), allocatable :: function_type_ham
      integer, dimension(:,:), allocatable :: function_type_s
      integer, dimension(:,:), allocatable :: function_type_hexton
      integer, dimension(:,:), allocatable :: function_type_hextoff
      integer, dimension(:,:), allocatable :: function_type_dens

      integer, dimension(:,:,:), allocatable :: symelement_short_atomic
      integer, dimension(:,:,:), allocatable :: symelement_elec
      integer, dimension(:,:,:), allocatable :: symelement_short_pair
      integer, dimension(:,:,:), allocatable :: symelement_ham
      integer, dimension(:,:,:), allocatable :: symelement_s
      integer, dimension(:,:,:), allocatable :: symelement_hextoff
      integer, dimension(:,:,:), allocatable :: symelement_hexton
      integer, dimension(:,:,:), allocatable :: symelement_dens

      real*8, dimension(:,:) , allocatable :: funccutoff_short_atomic
      real*8, dimension(:,:) , allocatable :: funccutoff_elec
      real*8, dimension(:,:) , allocatable :: funccutoff_short_pair
      real*8, dimension(:,:) , allocatable :: funccutoff_ham
      real*8, dimension(:,:) , allocatable :: funccutoff_s
      real*8, dimension(:,:) , allocatable :: funccutoff_hexton
      real*8, dimension(:,:) , allocatable :: funccutoff_hextoff
      real*8, dimension(:,:) , allocatable :: funccutoff_dens

      real*8, dimension(:,:) , allocatable :: eta_short_atomic
      real*8, dimension(:,:) , allocatable :: eta_elec
      real*8, dimension(:,:) , allocatable :: eta_short_pair
      real*8, dimension(:,:) , allocatable :: eta_ham
      real*8, dimension(:,:) , allocatable :: eta_s
      real*8, dimension(:,:) , allocatable :: eta_hexton
      real*8, dimension(:,:) , allocatable :: eta_hextoff
      real*8, dimension(:,:) , allocatable :: eta_dens

      real*8, dimension(:,:) , allocatable :: zeta_short_atomic
      real*8, dimension(:,:) , allocatable :: zeta_elec
      real*8, dimension(:,:) , allocatable :: zeta_short_pair
      real*8, dimension(:,:) , allocatable :: zeta_ham
      real*8, dimension(:,:) , allocatable :: zeta_s
      real*8, dimension(:,:) , allocatable :: zeta_hexton
      real*8, dimension(:,:) , allocatable :: zeta_hextoff
      real*8, dimension(:,:) , allocatable :: zeta_dens

      real*8, dimension(:,:) , allocatable :: lambda_short_atomic
      real*8, dimension(:,:) , allocatable :: lambda_elec
      real*8, dimension(:,:) , allocatable :: lambda_short_pair
      real*8, dimension(:,:) , allocatable :: lambda_ham
      real*8, dimension(:,:) , allocatable :: lambda_s
      real*8, dimension(:,:) , allocatable :: lambda_hexton
      real*8, dimension(:,:) , allocatable :: lambda_hextoff
      real*8, dimension(:,:) , allocatable :: lambda_dens

      real*8, dimension(:,:) , allocatable :: rshift_short_atomic
      real*8, dimension(:,:) , allocatable :: rshift_elec
      real*8, dimension(:,:) , allocatable :: rshift_short_pair
      real*8, dimension(:,:) , allocatable :: rshift_ham
      real*8, dimension(:,:) , allocatable :: rshift_s
      real*8, dimension(:,:) , allocatable :: rshift_hexton
      real*8, dimension(:,:) , allocatable :: rshift_hextoff
      real*8, dimension(:,:) , allocatable :: rshift_dens

      real*8 maxcutoff_short_atomic 
      real*8 maxcutoff_elec   
      real*8 maxcutoff_short_pair
      real*8 maxcutoff_ham
      real*8 maxcutoff_s
      real*8 maxcutoff_hexton
      real*8 maxcutoff_hextoff
      real*8 maxcutoff_dens

      contains 

      subroutine allocatesymfunctions()
!!
      use nnflags 
      use globaloptions 
!!
      implicit none
!!
      if(lshort.and.(nn_type_short.eq.1))then
        allocate(function_type_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        function_type_short_atomic(:,:)=0
        allocate(symelement_short_atomic(maxnum_funcvalues_short_atomic,2,nelem))
        symelement_short_atomic(:,:,:)=0
        allocate(funccutoff_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        funccutoff_short_atomic(:,:)=0.0d0
        allocate(eta_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        eta_short_atomic(:,:)=0.0d0
        allocate(zeta_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        zeta_short_atomic(:,:)=0.0d0
        allocate(lambda_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        lambda_short_atomic(:,:)=0.0d0
        allocate(rshift_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        rshift_short_atomic(:,:)=0.0d0
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        allocate(function_type_short_pair(maxnum_funcvalues_short_pair,npairs))
        function_type_short_pair(:,:)=0
        allocate(symelement_short_pair(maxnum_funcvalues_short_pair,2,npairs))
        symelement_short_pair(:,:,:)=0
        allocate(funccutoff_short_pair(maxnum_funcvalues_short_pair,npairs))
        funccutoff_short_pair(:,:)=0.0d0
        allocate(eta_short_pair(maxnum_funcvalues_short_pair,npairs))
        eta_short_pair(:,:)=0.0d0
        allocate(zeta_short_pair(maxnum_funcvalues_short_pair,npairs))
        zeta_short_pair(:,:)=0.0d0
        allocate(lambda_short_pair(maxnum_funcvalues_short_pair,npairs))
        lambda_short_pair(:,:)=0.0d0
        allocate(rshift_short_pair(maxnum_funcvalues_short_pair,npairs))
        rshift_short_pair(:,:)=0.0d0
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        allocate(function_type_elec(maxnum_funcvalues_elec,nelem))
        function_type_elec(:,:)=0
        allocate(symelement_elec(maxnum_funcvalues_elec,2,nelem))
        symelement_elec(:,:,:)=0
        allocate(funccutoff_elec(maxnum_funcvalues_elec,nelem))
        funccutoff_elec(:,:)=0.0d0
        allocate(eta_elec(maxnum_funcvalues_elec,nelem))
        eta_elec(:,:)=0.0d0
        allocate(zeta_elec(maxnum_funcvalues_elec,nelem))
        zeta_elec(:,:)=0.0d0
        allocate(lambda_elec(maxnum_funcvalues_elec,nelem))
        lambda_elec(:,:)=0.0d0
        allocate(rshift_elec(maxnum_funcvalues_elec,nelem))
        rshift_elec(:,:)=0.0d0
      endif
!!
      if(lnntb)then
        if(nntb_flag(1))then
          allocate(function_type_s(maxnum_funcvalues_s,npairs))
          function_type_s(:,:)=0
          allocate(symelement_s(maxnum_funcvalues_s,2,npairs))
          symelement_s(:,:,:)=0
          allocate(funccutoff_s(maxnum_funcvalues_s,npairs))
          funccutoff_s(:,:)=0.0d0
          allocate(eta_s(maxnum_funcvalues_s,npairs))
          eta_s(:,:)=0.0d0
          allocate(zeta_s(maxnum_funcvalues_s,npairs))
          zeta_s(:,:)=0.0d0
          allocate(lambda_s(maxnum_funcvalues_s,npairs))
          lambda_s(:,:)=0.0d0
          allocate(rshift_s(maxnum_funcvalues_s,npairs))
          rshift_s(:,:)=0.0d0

        endif
        if(nntb_flag(2))then
          allocate(function_type_hexton(maxnum_funcvalues_hexton,npairs))
          function_type_hexton(:,:)=0
          allocate(symelement_hexton(maxnum_funcvalues_hexton,2,npairs))
          symelement_hexton(:,:,:)=0
          allocate(funccutoff_hexton(maxnum_funcvalues_hexton,npairs))
          funccutoff_hexton(:,:)=0.0d0
          allocate(eta_hexton(maxnum_funcvalues_hexton,npairs))
          eta_hexton(:,:)=0.0d0
          allocate(zeta_hexton(maxnum_funcvalues_hexton,npairs))
          zeta_hexton(:,:)=0.0d0
          allocate(lambda_hexton(maxnum_funcvalues_hexton,npairs))
          lambda_hexton(:,:)=0.0d0
          allocate(rshift_hexton(maxnum_funcvalues_hexton,npairs))
          rshift_hexton(:,:)=0.0d0

        endif
        if(nntb_flag(3))then
          if(mode.eq.3)then
            allocate(function_type_hextoff(maxnum_funcvalues_hextoff,ntriplets))
            allocate(symelement_hextoff(maxnum_funcvalues_hextoff,2,ntriplets))
            allocate(funccutoff_hextoff(maxnum_funcvalues_hextoff,ntriplets))
            allocate(eta_hextoff(maxnum_funcvalues_hextoff,ntriplets))
            allocate(zeta_hextoff(maxnum_funcvalues_hextoff,ntriplets))
            allocate(lambda_hextoff(maxnum_funcvalues_hextoff,ntriplets))
            allocate(rshift_hextoff(maxnum_funcvalues_hextoff,ntriplets))
          else
            allocate(function_type_hextoff(maxnum_funcvalues_hextoff,1))
            allocate(symelement_hextoff(maxnum_funcvalues_hextoff,2,1))
            allocate(funccutoff_hextoff(maxnum_funcvalues_hextoff,1))
            allocate(eta_hextoff(maxnum_funcvalues_hextoff,1))
            allocate(zeta_hextoff(maxnum_funcvalues_hextoff,1))
            allocate(lambda_hextoff(maxnum_funcvalues_hextoff,1))
            allocate(rshift_hextoff(maxnum_funcvalues_hextoff,1))
          endif
          function_type_hextoff(:,:)=0
          symelement_hextoff(:,:,:)=0
          funccutoff_hextoff(:,:)=0.0d0
          eta_hextoff(:,:)=0.0d0
          zeta_hextoff(:,:)=0.0d0
          lambda_hextoff(:,:)=0.0d0
          rshift_hextoff(:,:)=0.0d0

        endif
        if(nntb_flag(4))then
          allocate(function_type_dens(maxnum_funcvalues_dens,npairs))
          function_type_dens(:,:)=0
          allocate(symelement_dens(maxnum_funcvalues_dens,2,npairs))
          symelement_dens(:,:,:)=0
          allocate(funccutoff_dens(maxnum_funcvalues_dens,npairs))
          funccutoff_dens(:,:)=0.0d0
          allocate(eta_dens(maxnum_funcvalues_dens,npairs))
          eta_dens(:,:)=0.0d0
          allocate(zeta_dens(maxnum_funcvalues_dens,npairs))
          zeta_dens(:,:)=0.0d0
          allocate(lambda_dens(maxnum_funcvalues_dens,npairs))
          lambda_dens(:,:)=0.0d0
          allocate(rshift_dens(maxnum_funcvalues_dens,npairs))
          rshift_dens(:,:)=0.0d0

        endif

      endif
      end subroutine allocatesymfunctions
 
      subroutine deallocatesymfunctions()
!!
      use nnflags 
      use globaloptions
!!
      implicit none

      if(lshort.and.(nn_type_short.eq.1))then
        deallocate(function_type_short_atomic)
        deallocate(symelement_short_atomic)
        deallocate(funccutoff_short_atomic)
        deallocate(eta_short_atomic)
        deallocate(zeta_short_atomic)
        deallocate(lambda_short_atomic)
        deallocate(rshift_short_atomic)
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        deallocate(function_type_short_pair)
        deallocate(symelement_short_pair)
        deallocate(funccutoff_short_pair)
        deallocate(eta_short_pair)
        deallocate(zeta_short_pair)
        deallocate(lambda_short_pair)
        deallocate(rshift_short_pair)
      endif

      if(lelec.and.(nn_type_elec.eq.1))then
        deallocate(function_type_elec)
        deallocate(symelement_elec)
        deallocate(funccutoff_elec)
        deallocate(eta_elec)
        deallocate(zeta_elec)
        deallocate(lambda_elec)
        deallocate(rshift_elec)
      endif
!!
      if(lnntb)then
        if(nntb_flag(1))then
          deallocate(function_type_s)
          deallocate(symelement_s)
          deallocate(funccutoff_s)
          deallocate(eta_s)
          deallocate(zeta_s)
          deallocate(lambda_s)
          deallocate(rshift_s)
        endif
        if(nntb_flag(2))then
          deallocate(function_type_hexton)
          deallocate(symelement_hexton)
          deallocate(funccutoff_hexton)
          deallocate(eta_hexton)
          deallocate(zeta_hexton)
          deallocate(lambda_hexton)
          deallocate(rshift_hexton)
        endif
        if(nntb_flag(3))then
          deallocate(function_type_hextoff)
          deallocate(symelement_hextoff)
          deallocate(funccutoff_hextoff)
          deallocate(eta_hextoff)
          deallocate(zeta_hextoff)
          deallocate(lambda_hextoff)
          deallocate(rshift_hextoff)
        endif
        if(nntb_flag(4))then
          deallocate(function_type_dens)
          deallocate(symelement_dens)
          deallocate(funccutoff_dens)
          deallocate(eta_dens)
          deallocate(zeta_dens)
          deallocate(lambda_dens)
          deallocate(rshift_dens)
        endif

      endif
      end subroutine deallocatesymfunctions

      subroutine distribute_symfunctions()
!!
      use mpi_mod
      use nnflags 
      use nnham
      use globaloptions
!!
      implicit none
      integer ndim
!!
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(function_type_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(function_type_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_short_pair,&
          maxnum_funcvalues_short_pair*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(function_type_elec,&
          maxnum_funcvalues_elec*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_elec,&
          maxnum_funcvalues_elec*nelem*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lnntb)then
        if(nntb_flag(0))then
          call mpi_bcast(function_type_ham,&
            maxnum_funcvalues_ham*npairs,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(symelement_ham,&
            maxnum_funcvalues_ham*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(funccutoff_ham,&
            maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(eta_ham,&
            maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(zeta_ham,&
            maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(lambda_ham,&
            maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(rshift_ham,&
            maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        endif

        if(nntb_flag(1))then
          call mpi_bcast(function_type_s,&
            maxnum_funcvalues_s*npairs,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(symelement_s,&
            maxnum_funcvalues_s*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(funccutoff_ham,&
            maxnum_funcvalues_s*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(eta_s,&
            maxnum_funcvalues_s*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(zeta_s,&
            maxnum_funcvalues_s*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(lambda_s,&
            maxnum_funcvalues_s*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(rshift_s,&
            maxnum_funcvalues_s*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        endif
        if(nntb_flag(2))then
          call mpi_bcast(function_type_hexton,&
            maxnum_funcvalues_hexton*npairs,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(symelement_hexton,&
            maxnum_funcvalues_hexton*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(funccutoff_hexton,&
            maxnum_funcvalues_hexton*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(eta_hexton,&
            maxnum_funcvalues_hexton*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(zeta_hexton,&
            maxnum_funcvalues_hexton*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(lambda_hexton,&
            maxnum_funcvalues_hexton*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(rshift_hexton,&
            maxnum_funcvalues_hexton*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        endif
        if(nntb_flag(3))then
          if(mode.eq.3)then
            ndim=ntriplets
          else
            ndim=1
          endif
          call mpi_bcast(function_type_hextoff,&
            maxnum_funcvalues_hextoff*ndim,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(symelement_hextoff,&
            maxnum_funcvalues_hextoff*ndim,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(funccutoff_hextoff,&
            maxnum_funcvalues_hextoff*ndim,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(eta_hextoff,&
            maxnum_funcvalues_hextoff*ndim,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(zeta_hextoff,&
            maxnum_funcvalues_hextoff*ndim,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(lambda_hextoff,&
            maxnum_funcvalues_hextoff*ndim,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(rshift_hextoff,&
            maxnum_funcvalues_hextoff*ndim,mpi_real8,0,mpi_comm_world,mpierror)
        endif
        if(nntb_flag(4))then
          call mpi_bcast(function_type_dens,&
            maxnum_funcvalues_dens*npairs,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(symelement_dens,&
            maxnum_funcvalues_dens*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
          call mpi_bcast(funccutoff_dens,&
            maxnum_funcvalues_dens*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(eta_dens,&
            maxnum_funcvalues_dens*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(zeta_dens,&
            maxnum_funcvalues_dens*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(lambda_dens,&
            maxnum_funcvalues_dens*npairs,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(rshift_dens,&
            maxnum_funcvalues_dens*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        endif

      endif
      end subroutine distribute_symfunctions

      end module symfunctions 
