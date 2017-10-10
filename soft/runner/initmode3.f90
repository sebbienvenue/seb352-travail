!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! Purpose: read and distribute scaling data and weights, determine maxcutoffs and num_pairs
!!
      subroutine initmode3(num_atoms,num_pairs,zelem,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        eshortmin,eshortmax,chargemin,chargemax,&
        lattice,xyzstruct,lperiodic&
        )
!!
      use fileunits
      use mpi_mod
      use nnflags
      use globaloptions
      use symfunctions
      use nnewald
      use nnshort_atomic
      use nnshort_pair
      use predictionoptions
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5                                              ! internal
      integer num_atoms                                                   ! in 
      integer num_pairs                                                   ! out 
      integer zelem(max_num_atoms)                                        ! in 
      integer jcount                                                      ! internal
      integer paircount                                                   ! internal
!!
      real*8 lattice(3,3)                                                 ! in 
      real*8 xyzstruct(3,max_num_atoms)                                   ! in 
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)  ! out 
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)  ! out 
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)   ! out 
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)                  ! out 
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)                  ! out 
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)                   ! out 
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)     ! out 
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)     ! out 
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)      ! out 
      real*8 eshortmin                                                    ! out 
      real*8 eshortmax                                                    ! out 
      real*8 chargemin(nelem)                                             ! out 
      real*8 chargemax(nelem)                                             ! out 
      real*8 rdummy(nelem)                                                ! internal
      real*8 dummy                                                        ! internal
!!
      logical lperiodic                                                   ! in 
!!
!!=====================================================================
!! determine maxcutoffs 
!!=====================================================================
      maxcutoff_short_atomic       =0.0d0
      maxcutoff_short_pair         =0.0d0
      maxcutoff_elec               =0.0d0
!! get maxcutoff_short_atomic
      if(lshort.and.(nn_type_short.eq.1))then
        do i2=1,nelem
          do i1=1,num_funcvalues_short_atomic(i2)
            maxcutoff_short_atomic=max(maxcutoff_short_atomic,funccutoff_short_atomic(i1,i2))
          enddo ! i1
        enddo
      endif ! lshort
!! get maxcutoff_short_pair
      if(lshort.and.(nn_type_short.eq.2))then
        do i2=1,npairs
          do i1=1,num_funcvalues_short_pair(i2)
           maxcutoff_short_pair=max(maxcutoff_short_pair,funccutoff_short_pair(i1,i2))
          enddo ! i1
        enddo
      endif ! lshort
!! get maxcutoff_elec
      if(lelec.and.(nn_type_elec.eq.1))then
        do i2=1,nelem
          do i1=1,num_funcvalues_elec(i2)
            maxcutoff_elec=max(maxcutoff_elec,funccutoff_elec(i1,i2))
          enddo ! i1
        enddo ! i2
      endif
!!
!!=====================================================================
!! For nn_type_short 2 we need to determine num_pairs
!!=====================================================================
      if((lshort.and.(nn_type_short.eq.2)).or.(lnntb))then
        if(mpirank.eq.0)then
          call getnumpairs(num_atoms,num_pairs,zelem,&
            maxcutoff_short_pair,lattice,xyzstruct,lperiodic)
        endif ! 'mpirank.eq.0
        call mpi_bcast(num_pairs,1,mpi_integer,0,mpi_comm_world,mpierror)
      endif
!!
!!=====================================================================
!! read and distribute the scaling data from file scaling.data
!!=====================================================================
      minvalue_short_atomic(:,:)   =0.0d0
      maxvalue_short_atomic(:,:)   =0.0d0
      avvalue_short_atomic(:,:)    =0.0d0
      minvalue_short_pair(:,:)     =0.0d0
      maxvalue_short_pair(:,:)     =0.0d0
      avvalue_short_pair(:,:)      =0.0d0
      minvalue_elec(:,:) =0.0d0
      maxvalue_elec(:,:) =0.0d0
      avvalue_elec(:,:)  =0.0d0
      chargemin(:)       =0.0d0
      chargemax(:)       =0.0d0
      if(mpirank.eq.0)then
        if(lshort.and.(nn_type_short.eq.1))then
          call readscale(nelem,1,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
            eshortmin,eshortmax,rdummy,rdummy)
        endif ! lshort
        if(lshort.and.(nn_type_short.eq.2))then
          call readscale(npairs,2,&
            maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
            minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
            eshortmin,eshortmax,rdummy,rdummy)
        endif ! lshort
        if(lelec.and.(nn_type_elec.eq.1))then
          call readscale(nelem,3,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            minvalue_elec,maxvalue_elec,avvalue_elec,&
            dummy,dummy,chargemin,chargemax)
        endif ! lelec
      endif ! mpirank.eq.0
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(minvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eshortmin,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eshortmax,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(minvalue_short_pair,npairs*maxnum_funcvalues_short_pair,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvalue_short_pair,npairs*maxnum_funcvalues_short_pair,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvalue_short_pair,npairs*maxnum_funcvalues_short_pair,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eshortmin,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eshortmax,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(minvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(chargemin,nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(chargemax,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!!=====================================================================
!! read and distribute the weight parameters 
!!=====================================================================
      if(mpirank.eq.0)then
        if(lshort.and.(nn_type_short.eq.1))then
          call readweights(0,nelem,&
            maxnum_weights_short_atomic,num_weights_short_atomic,&
            weights_short_atomic)
        endif ! lshort
        if(lshort.and.(nn_type_short.eq.2))then
          call readweights(2,npairs,&
            maxnum_weights_short_pair,num_weights_short_pair,&
            weights_short_pair)
        endif ! lshort
        if(lelec.and.(nn_type_elec.eq.1))then
          call readweights(1,nelem,&
            maxnum_weights_elec,num_weights_elec,&
            weights_elec)
        endif
      endif ! mpirank.eq.0
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(weights_short_atomic,nelem*maxnum_weights_short_atomic,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif ! lshort
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(weights_short_pair,npairs*maxnum_weights_short_pair,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif ! lshort
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(weights_elec,nelem*maxnum_weights_elec,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!!=====================================================================
!! write weights for debugging
!!=====================================================================
      if(ldebug.and.(mpisize.eq.1))then
        if(lshort.and.(nn_type_short.eq.1))then
          write(ounit,'(a82)')'                    element   from layer node    to layer node              weight'
          do i1=1,nelem !'
            jcount=0
            do i2=1,num_layers_short_atomic(i1)
              do i3=1,nodes_short_atomic(i2-1,i1)
                do i4=1,nodes_short_atomic(i2,i1)
                  jcount=jcount+1
                  write(ounit,'(a20,i7,8x,2i5,7x,2i5,f20.14)')'WEIGHTS SHORT      ',&
                    nucelem(i1),i2-1,i3,i2,i4,weights_short_atomic(jcount,i1)
                enddo ! i4'
              enddo ! i3
              do i3=1,nodes_short_atomic(i2,i1)
                jcount=jcount+1
                write(ounit,'(a20,i7,25x,2i5,f20.14)')'BIAS WEIGHTS SHORT ',&
                  nucelem(i1),i2,i3,weights_short_atomic(jcount,i1)
              enddo
            enddo ! i2
          enddo ! i1
          write(ounit,*)'-------------------------------------------------------------'
        endif ! lshort
        if(lshort.and.(nn_type_short.eq.2))then
          write(ounit,'(a82)')'                    element   from layer node    to layer node              weight'
          paircount = 0 
          do i1=1,nelem
           do i2=i1,nelem
            paircount = paircount +1
            jcount=0
            do i3=1,num_layers_short_pair(paircount)
              do i4=1,nodes_short_pair(i3-1,paircount)
                do i5=1,nodes_short_pair(i3,paircount)
                  jcount=jcount+1
                  write(ounit,'(a20,i7,i7,8x,2i5,7x,2i5,f20.14)')'WEIGHTS SHORT      ',&
                    nucelem(i1),nucelem(i2),i3-1,i4,i3,i5,weights_short_pair(jcount,paircount)
                enddo ! i5'
              enddo ! i4
              do i4=1,nodes_short_pair(i3,paircount)
                jcount=jcount+1
                write(ounit,'(a20,i7,i7,25x,2i5,f20.14)')'BIAS WEIGHTS SHORT ',&
                  nucelem(i1),nucelem(i2),i3,i4,weights_short_pair(jcount,paircount)
              enddo
             enddo ! i3
            enddo ! i2
          enddo ! i1
         write(ounit,*)'-------------------------------------------------------------'
        endif ! lshort
        if(lelec.and.(nn_type_elec.eq.1))then
          write(ounit,'(a82)')'                    element   from layer node    to layer node              weight'
          do i1=1,nelem !'
            jcount=0
            do i2=1,num_layers_elec(i1)
              do i3=1,nodes_elec(i2-1,i1)
                do i4=1,nodes_elec(i2,i1)
                  jcount=jcount+1
                  write(ounit,'(a21,i6,8x,2i5,7x,2i5,f20.14)')'WEIGHTS CHARGE      ',&
                    nucelem(i1),i2-1,i3,i2,i4,weights_elec(jcount,i1)
                enddo ! i4'
              enddo ! i3
              do i3=1,nodes_elec(i2,i1)
                jcount=jcount+1
                write(ounit,'(a21,i6,25x,2i5,f20.14)')'BIAS WEIGHTS CHARGE ',&
                  nucelem(i1),i2,i3,weights_elec(jcount,i1)
              enddo
            enddo ! i2
          enddo ! i1
          write(ounit,*)'-------------------------------------------------------------'
        endif ! lelec
      endif ! ldebug
!!
!!=====================================================================
!! dump full configuration of RuNNer to file for external MD program (not completed yet) 
!!=====================================================================
      if(lpreparemd)then
        call preparemd()
      endif
!!
      return
      end
