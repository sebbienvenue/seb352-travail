!!=====================================================================
!!SK: this is a checkup if the i-pi wrapper is running 
!!SK: like in FHIaims-ipi-wrapper (from MR) 
!!SK: still needs to be tested
!!=====================================================================
!!
function drv_check(info)
    implicit none
    logical drv_check
    integer, intent(in) :: info
!!
    call mpi_bcast(info) 
    drv_check = .false.
    if (info.le.0) drv_check = .true.
end function
!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - main.f90 if ipi is used
!!
      subroutine predict_ipi()
!!
      use mpi_mod
      use fileunits
      use predictionoptions 
      use nnflags 
      use globaloptions
      use symfunctions 
      use timings
      use nnewald
      use nnshort_atomic
      use nnshort_pair
      use runneripi
!!
      implicit none
!!=====================================================================
!! Variables for RuNNer
!!=====================================================================
      integer zelem(max_num_atoms)                           ! internal
      integer num_atoms                                      ! internal
      integer num_pairs                                      ! internal
      integer num_atoms_element(nelem)                       ! internal 
      integer i1,i2,i3                                       ! internal
      integer counter                                        ! internal
      integer pairs_charge(2,max_num_pairs)                  ! internal 
!!
      real*8 lattice(3,3)                                    ! internal
      real*8 xyzstruct(3,max_num_atoms)                      ! internal
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)  ! internal
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)  ! internal
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)   ! internal
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)     ! internal
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)     ! internal
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)      ! internal
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)  ! internal
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)  ! internal
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)   ! internal
      real*8, dimension(:,:)  , allocatable :: sens          ! internal
      real*8, dimension(:,:)  , allocatable :: sense         ! internal
      real*8 eshort                                          ! internal
      real*8 volume                                          ! internal
!! DFT data (not necessarily provided in predicition mode)
      real*8 totalcharge                                     ! internal
      real*8 totalenergy                                     ! internal
      real*8 totalforce(3,max_num_atoms)                     ! internal
      real*8 atomcharge(max_num_atoms)                       ! internal
      real*8 atomenergy(max_num_atoms)                       ! internal
!! NN data
      real*8 nntotalenergy                                   ! internal 
      real*8 nnshortenergy                                   ! internal 
      real*8 nntotalcharge                                   ! internal 
      real*8 nnshortforce(3,max_num_atoms)                   ! internal 
      real*8 nntotalforce(3,max_num_atoms)                   ! internal 
      real*8 nnatomcharge(max_num_atoms)                     ! internal 
      real*8 nnatomenergy(max_num_atoms)                     ! internal
      real*8 nnelecforce(3,max_num_atoms)                    ! internal 
      real*8 nnpairenergy(max_num_pairs)                     ! internal
      real*8 nnstress(3,3)                                   ! internal 
      real*8 nnstress_short(3,3)                             ! internal 
      real*8 nnstress_elec(3,3)                              ! internal 
      real*8 nnelecenergy                                    ! internal 
      real*8 atomenergysum                                   ! internal
      real*8 forcevdw(3,max_num_atoms)                       ! internal 
      real*8 energyvdw                                       ! internal
      real*8 eshortmin                                       ! internal
      real*8 eshortmax                                       ! internal
      real*8 chargemin(nelem)                                ! internal
      real*8 chargemax(nelem)                                ! internal
      real*8 pressure                                        ! internal
      real*8 forcesum(3)                                     ! internal
      real*8 au2gpa                                          ! internal
!!
      character*2 elementsymbol(max_num_atoms)               ! internal
!!
      logical lperiodic                                      ! internal
!!
!!=====================================================================
!!SK: variables i-pi-wrapper 
!!=====================================================================
!!
    integer, parameter :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
    integer socket, inet, port        ! socket ID & address of the server
    character*1024 :: host
!!
    character*12 :: header
    logical :: isinit=.true., hasdata=.false., drv_check
    character*1024 :: parbuffer
!!
    integer nat, siz, info, thisrep  !nat=number of atoms 
    integer ccmd, i, replicaid 
!!
!!SK: system specific cell parameters i-pi internal 
!!SK: callh= lattice vectors
!!SK: callih= reciprocal lattice vectors
!!SK: vir= stress 
!!SK: pot= Epot
    real*8 :: cellh(3,3), cellih(3,3), vir(3,3), pot
!!
    real*8 combuf(3,max_num_atoms) !!SK: this we need for xyz from i-pi 
    real*8 forces(3,max_num_atoms) !!SK: this we need to pass forces to i-pi
    character*1024 srvaddress      !!SK: i-pi internal server adress
    character*100000 extra_string  !!SK: to pass other information e.g. charges to i-pi
!!=====================================================================
      call zerotime(daymode3,timemode3start,timemode3end) 
      call abstime(timemode3start,daymode3)
!!=====================================================================
!! initializations RuNNer 
!!=====================================================================
      if((mpirank.eq.0).and.(.not.lmd))then
        write(ounit,*)'Prediction Mode i-pi'
      endif
      timeio            = 0.0d0
      timeoutput        = 0.0d0
      eshort            = 0.0d0
      nnshortenergy     = 0.0d0
      nntotalenergy     = 0.0d0
      nnelecenergy      = 0.0d0
      nnatomenergy(:)   = 0.0d0
      nntotalforce(:,:) = 0.0d0
      nnshortforce(:,:) = 0.0d0
      nnelecforce(:,:)  = 0.0d0
      nntotalcharge     = 0.0d0
      atomenergysum     = 0.0d0
      au2gpa=29419.844d0 ! 1 Ha/Bohr3 = 29419.844 GPa
      nnstress_short(:,:)=0.0d0
      nnstress_elec(:,:) =0.0d0
!!
      if(lshort.and.(nn_type_short.eq.1))then
        allocate(sens(nelem,maxnum_funcvalues_short_atomic))
      elseif(lshort.and.(nn_type_short.eq.2))then
        allocate(sens(npairs,maxnum_funcvalues_short_pair))
      endif
!! we need to allocate sense also for nn_type_elec 3 and 4 because of subroutine call below
      if(lelec.and.(nn_type_elec.eq.1).or.(nn_type_elec.eq.3).or.(nn_type_elec.eq.4))then
        allocate(sense(nelem,maxnum_funcvalues_elec))
      endif
!!
      if(ldebug.and.(mpisize.gt.1).and.(.not.lmd))then
        ldebug=.false.
        write(ounit,*)'### WARNING ### ldebug is switched off for parallel case'
      endif !'
!!
!!=====================================================================
!! timining initialization of mode 3 
!!=====================================================================
      if(lfinetime)then
        dayio=0
        call abstime(timeiostart,dayio)
      endif ! lfinetime
!!
!!=====================================================================
!! read structure and derive all structure related data, also mpi distribution
!!SK: i-pi just provides xyz coordinates and number of atoms, but no
!!SK: information about the elements - so we need to initialize the structure
!!SK: USERS: make sure that you use the same order of elements in input.data
!!SK: and you xyz structure for i-pi
!!=====================================================================
      call getstructure_mode3(num_atoms,zelem,&
        num_atoms_element,lattice,xyzstruct,&
        totalenergy,totalcharge,totalforce,atomenergy,atomcharge,&
        elementsymbol,lperiodic&
        )
!!
!!=====================================================================
!! read and distribute scaling data and weights, determine num_pairs if needed 
!!=====================================================================
      call initmode3(num_atoms,num_pairs,zelem,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        eshortmin,eshortmax,chargemin,chargemax,&
        lattice,xyzstruct,lperiodic& 
        )
!!
!!=====================================================================
!! timining initialization of mode 3 
!!=====================================================================
      if(lfinetime)then
        call abstime(timeioend,dayio)
        timeio=timeioend-timeiostart
      endif ! lfinetime
!!
!!=====================================================================
!!SK: at this point the i-pi wrapper starts - socket specific stuff 
!!SK: we don't need to change anything
!!SK: I still need to test inet connection - but it should be no problem
!!=====================================================================
    replicaid=0
    thisrep=0
    srvaddress=trim(desc_str)//":"//trim(ipistring)//":"//trim(ipisocket)
    inet=1    !!SK: 1 if we want a internet socket (remote name or localhost or ip)   0 if we want UNIX socket
    host=srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)//achar(0)  !!SK: host name (or unix socket name)
    read(srvaddress(INDEX(srvaddress,':',back=.true.)+1:),*) port     !!SK: port number (kind of ignored for unix sockets)
    
    IF (mpirank.eq.0) write(*,*) " @ DRIVER MODE: Connecting to host:port", &
             trim(srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)), port
    
    IF (srvaddress(1:INDEX(srvaddress,':')-1).eq.('UNIX')) THEN
      inet=0
      host=srvaddress(6:INDEX(srvaddress,':',back=.true.)-1)//achar(0)    
    ENDIF

IF (mpirank.eq.0) call open_socket(socket, inet, port, host)          

!!=====================================================================
!!SK: driver loop starts
!!SK: TO DO: drv_check not working yet, thus commented out
!!=====================================================================
    
    driver_loop: DO
      if (mpirank.eq.0) call readbuffer(socket, header, MSGLEN, info)
!      if (drv_check(info)) exit  
      
      if (mpirank.eq.0) write(*,*) " @ DRIVER MODE: Message from server: ", header
      if (trim(header) == "STATUS") then
         if (mpirank.eq.0) then  
            if (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN, info)
            else if (isinit) then
               call writebuffer(socket,"READY       ",MSGLEN, info)
            else 
               call writebuffer(socket,"NEEDINIT       ",MSGLEN, info)
            endif
         endif
!         if (drv_check(info)) exit  
      else if (trim(header)=="INIT") then
         if (mpirank.eq.0) call readbuffer(socket, thisrep, 4, info) !!SK: reading replica info 
!         if (drv_check(info)) exit 
         if (mpirank.eq.0) write(*,*) " @ DRIVER MODE: Receiving replica",thisrep  
         replicaid = thisrep
         if (mpirank.eq.0) then
            call readbuffer(socket, nat, 4, info) 
            call readbuffer(socket, parbuffer, nat, info)
        endif
!        if (drv_check(info)) exit  

!!==================================================
!!SK: here we get informtion from ipi
!!SK: USERS: make sure that you use the same order
!!SK: in xyz and input.data with i-pi and check the 
!!SK: units in the i-pi input xml file
!!=================================================
        
        isinit=.true.
      else if (trim(header) == "POSDATA") then !!SK: driver is sending the positions of the atoms. 
         if (mpirank.eq.0) then        
            call readbuffer(socket, cellh, 9*8, info)    !!SK: lattice vectors
            call readbuffer(socket, cellih, 9*8, info)   !!SK: reziprocal lattice vectors
            call readbuffer(socket, nat, 4, info)        !!SK: number of atoms
         endif         
!         if (drv_check(info)) exit  

!!SK: DO WE NEED THAT FOR MPI? 
!!         call mp_bcast_mat(cellh) ! 2d matrix 
!!         call mp_bcast_mat(cellih) ! 2d matrix
!!         call mp_bcast_int(nat) ! integer

!!SK: Here we get the coordinates from i-pi with combuf
         if (mpirank.eq.0) call readbuffer(socket, combuf, nat*3*8, info)  

!!         if (drv_check(info)) exit
      
!!         call mp_bcast_arr(combuf) ! real array

!!====================================================
!!SK: Here the work begins! 
!!SK: pass lattice, number of atoms, xyz coord from 
!!SK: i-pi to RuNNer
!!====================================================
!!SK: Be careful with the cell parameter storage!
!!SK: RuNNer: rows (ax, ay, az; bx, by, bz; cx, cy, cz)
!!SK: i-pi:  rows (ax, ay, az; bx, by, bz; cx, cy, cz)


!SK: TO DO: checkup for nat = num_atoms
!SK: TO DO: checkup for periodic structure

           lattice=cellh 
           num_atoms = nat
           xyzstruct = combuf 
         if (mpirank.eq.0) write(*,*) " @ DRIVER MODE: Received positions "

!!====================================================
!!SK: Here we calculate everything we need (energies, 
!!SK  forces with RuNNer
!!SK: The subroutines are from predict.f90
!!SK: 'Things' we don't need are just commented out (e.g. output) 
!!====================================================
        
        if(lshort)then
          if(nn_type_short.eq.1)then 
            call predictionshortatomic(&
              num_atoms,num_atoms_element,zelem,&
              lattice,xyzstruct,&
              minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
              eshortmin,eshortmax,&
              nntotalenergy,nnshortforce,&
              nnatomenergy,nnshortenergy,nnstress_short,&
              atomenergysum,sens,lperiodic)
          elseif(nn_type_short.eq.2)then
            call predictionshortpair(&
              num_atoms,num_atoms_element,zelem,&
              lattice,xyzstruct,&
              minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
              eshortmin,eshortmax,&
              nntotalenergy,nnshortforce,&
              nnatomenergy,nnpairenergy,nnshortenergy,&
              nnstress_short,pairs_charge,&
              atomenergysum,sens,lperiodic)
          endif
        endif ! lshort
!! electrostatic part, only if charges are not used as second output of atomic NN
        if(lelec.and.((nn_type_elec.eq.1).or.(nn_type_elec.eq.3)&
           .or.(nn_type_elec.eq.4)))then
          call predictionelectrostatic(&
            num_atoms,zelem,&
            minvalue_elec,maxvalue_elec,avvalue_elec,&
            lattice,xyzstruct,&
            nntotalcharge,nnatomcharge,&
            chargemin,chargemax,nnelecenergy,&
            nnelecforce,nnstress_elec,sense,lperiodic)
        else
!! if no electrostatics is used, initialize charge array as zero for proper output
          nnatomcharge(:)=0.0d0 
        endif ! lelec
!!
!!======================================================================
!! combine the short range and electrostatic energies
!!======================================================================
      nntotalenergy=nnshortenergy+nnelecenergy
!!
!!======================================================================
!! add energies of free atoms
!!======================================================================
      if(lremoveatomenergies.and.lshort)then
        call addatoms(num_atoms,&
          zelem,num_atoms_element,&
          atomenergysum,nnatomenergy)
        nntotalenergy=nntotalenergy+atomenergysum
      endif ! lremoveatomenergies
!!
!!======================================================================
!! print short range and electrostatic forces separately for debugging'
!!======================================================================
      if((mpirank.eq.0).and.(.not.lmd))then
        if(ldoforces)then
          write(ounit,*)'-------------------------------------------------------------'
          if(ldebug)then
            if(lshort)then
              write(ounit,*)'NN short range forces'
              do i1=1,num_atoms
                write(ounit,'(i6,3f18.8)')i1,nnshortforce(1,i1),&
                nnshortforce(2,i1),nnshortforce(3,i1)
              enddo
            endif
            if(lelec)then
              write(ounit,*)'NN electrostatic forces'
              do i1=1,num_atoms
                write(ounit,'(i6,3f18.8)')i1,nnelecforce(1,i1),&
                nnelecforce(2,i1),nnelecforce(3,i1)
              enddo
            endif
            write(ounit,*)'-------------------------------------------------------------'
          endif ! ' ldebug
        endif ! ldoforces
      endif !' mpirank.eq.0
!!
!!======================================================================
!! combination of short-range and electrostatic forces
!! must be done after the separate output of short and electrostatic forces
!!======================================================================
      if(ldoforces)then
        nntotalforce(:,:)=nnshortforce(:,:)+nnelecforce(:,:)
      endif ! ldoforces
!!
!!======================================================================
!! calculate the volume, needed also for stress
!!======================================================================
      if(lperiodic)then
        volume=0.0d0
        call getvolume(lattice,volume)
        if((mpirank.eq.0).and.(.not.lmd))then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a8,f20.8,a6)')' volume ',volume,' Bohr3'
        endif ! mpirank.eq.0
      endif !' lperiodic
!!
!!======================================================================
!! combination of short-range and electrostatic stress 
!!======================================================================
      if(ldostress.and.lperiodic)then
        nnstress(:,:)=nnstress_short(:,:)+nnstress_elec(:,:)
        nnstress(:,:)=nnstress(:,:)/volume
      endif ! ldostress
!!
!!======================================================================
!! add vdW if requested FIXME: This is not working yet
!!======================================================================
        energyvdw=0.0d0
        if(nn_type_vdw.eq.1)then
          call addvdw(num_atoms,zelem,&
            xyzstruct,lattice,energyvdw,forcevdw,&
            lperiodic) 
          nntotalenergy=nntotalenergy+energyvdw
          do i1=1,num_atoms
            do i2=1,3
              nntotalforce(i2,i1)=nntotalforce(i2,i1)+forcevdw(i2,i1)
            enddo
          enddo
        endif
!!
!!======================================================================
!! now write results to files 
!!======================================================================
!!      if(lfinetime)then
!!        dayoutput=0
!!        call abstime(timeoutputstart,dayoutput)
!!      endif ! lfinetime
!!!!
!!      if((mpirank.eq.0).and.(.not.lmd))then
!!        open(nneunit,file='energy.out',form='formatted',status='replace')
!!          write(nneunit,'(f18.8)')nntotalenergy
!!        close(nneunit)
!!!!
!!!! delete nnforces just in case an old file is present
!!        open(nnfunit,file='nnforces.out',form='formatted',status='unknown')
!!        close(nnfunit,status='delete')
!!        if(ldoforces)then
!!          open(nnfunit,file='nnforces.out',form='formatted',status='replace')
!!          do i1=1,num_atoms !'
!!            write(nnfunit,'(3f18.8)')(nntotalforce(i2,i1),i2=1,3)
!!          enddo
!!          close(nnfunit)
!!        endif ! ldoforces
!!
!!======================================================================
!! check sum of forces if requested 
!!======================================================================
        if(lcheckf)then
          forcesum(:)=0.0d0
          do i1=1,num_atoms
            do i2=1,3
              forcesum(i2)=forcesum(i2)+nntotalforce(i2,i1)
            enddo ! i2
          enddo ! i1
          write(ounit,*)'Sum of forces in x, y and z:'
          write(ounit,'(3f14.8)')forcesum(1),forcesum(2),forcesum(3)
          do i1=1,3
            if(abs(forcesum(i1)).gt.0.000001d0)then
              write(ounit,*)'Error in forces: ',forcesum(i1)
              stop
            endif
          enddo ! i1
        endif ! lcheckf
!!
!! delete nnstress just in case an old file is present
!!        open(nnsunit,file='nnstress.out',form='formatted',status='unknown')
!!        close(nnsunit,status='delete')
!!        if(ldostress.and.lperiodic)then
!!          open(nnsunit,file='nnstress.out',form='formatted',status='replace')
!!          do i1=1,3 !'
!!            write(nnsunit,'(3f18.8)')(nnstress(i1,i2),i2=1,3)
!!          enddo ! i1
!!          close(nnsunit)
!!        endif ! ldostress
!!
!! delete nnatoms just in case an old file is present
!!        open(nnaunit,file='nnatoms.out',form='formatted',status='unknown')
!!        close(nnaunit,status='delete')
!!        open(nnaunit,file='nnatoms.out',form='formatted',status='replace')
!!          do i1=1,num_atoms
!!            write(nnaunit,'(i6,x,a2,x,2f14.6)')i1,elementsymbol(i1),&
!!            nnatomcharge(i1),nnatomenergy(i1)
!!          enddo ! i1
!!        close(nnaunit)
!!
!!        if(nn_type_short.eq.2)then
!!          open(nnpunit,file='nnpairs.out',form='formatted',status='unknown')
!!          close(nnpunit,status='delete')
!!          open(nnpunit,file='nnpairs.out',form='formatted',status='replace')
!!            do i1=1,num_pairs 
!!              write(nnpunit,'(i6,x,a2,x,a2,x,2f14.6)')i1,&
!!                element(elementindex(pairs_charge(1,i1))),&
!!                element(elementindex(pairs_charge(2,i1))),&
!!                0.d0,nnpairenergy(i1)
!!            enddo ! i1
!!          close(nnpunit)
!!        endif
!!
!!======================================================================
!! write output data in RuNNer format and put in all results we have
!!======================================================================
!!        open(outunit,file='output.data',form='formatted',status='unknown')
!!        close(outunit,status='delete')
!!        open(outunit,file='output.data',form='formatted',status='replace')
!!          write(outunit,*)'begin'
!!          write(outunit,*)'comment File written by RuNNer in prediction mode'
!!          if(lperiodic)then !'
!!            write(outunit,'(a8,3f14.6)')'lattice ',(lattice(1,i1),i1=1,3)
!!            write(outunit,'(a8,3f14.6)')'lattice ',(lattice(2,i1),i1=1,3)
!!            write(outunit,'(a8,3f14.6)')'lattice ',(lattice(3,i1),i1=1,3)
!!          endif
!!          do i1=1,num_atoms
!!            write(outunit,'(a5,3f14.6,x,a2,x,5f14.6)')&
!!            'atom ',(xyzstruct(i2,i1),i2=1,3),&
!!            elementsymbol(i1),nnatomcharge(i1),nnatomenergy(i1),&
!!            (nntotalforce(i2,i1),i2=1,3)
!!          enddo
!!          write(outunit,'(a7,f20.6)')'energy ',nntotalenergy
!!          write(outunit,'(a7,f14.6)')'charge ',nntotalcharge
!!          write(outunit,*)'end'
!!        close(outunit)
!!      endif ! mpirank.eq.0
!!
!!======================================================================
!! write results to standard out for process 0
!!======================================================================
!!      if((mpirank.eq.0).and.(.not.lmd))then
!!        write(ounit,*)'-------------------------------------------------------------'
!!!!        if((.not.lshort).and.(lremoveatomenergies))then
!!!!          write(ounit,*)'WARNING: short range energy is just sum of atomic energies because lshort=F'
!!!!        endif
!!        write(ounit,'(a29,f20.8,a3)')' NN sum of free atom energies ',atomenergysum,' Ha'
!!        write(ounit,'(a29,f20.8,a3)')' NN short range energy        ',nntotalenergy-nnelecenergy-atomenergysum,' Ha'
!!        write(ounit,'(a29,f20.8,a3)')' NN electrostatic energy      ',nnelecenergy,' Ha'
!!        write(ounit,'(a29,f20.8,a3)')' van der Waals energy         ',energyvdw,' Ha'
!!        write(ounit,'(a29,f20.8,a3)')' NNenergy                     ',nntotalenergy,' Ha'
!!        write(ounit,*)'-------------------------------------------------------------'
!!        if(nn_type_short.eq.1)then
!!          do i1=1,num_atoms
!!            write(ounit,'(a14,i6,x,a2,x,f14.6)')' NNatomenergy ',i1,elementsymbol(i1),&
!!              nnatomenergy(i1)
!!          enddo ! i1
!!        elseif(nn_type_short.eq.2)then
!!          do i1=1,num_pairs
!!            write(ounit,'(a14,i6,x,a2,x,a2,x,f14.6)')'NNpairenergy ',i1,&
!!              element(elementindex(pairs_charge(1,i1))),&
!!              element(elementindex(pairs_charge(2,i1))),&
!!              nnpairenergy(i1)
!!          enddo ! i1
!!        endif
!!        write(ounit,*)'-------------------------------------------------------------'
!!        if(lelec)then
!!          do i1=1,num_atoms
!!            write(ounit,'(a,i6,f14.8)')' NNcharge ',i1,nnatomcharge(i1)
!!          enddo
!!        endif ! 
!!
!!        write(ounit,*)'-------------------------------------------------------------'
!!        write(ounit,*)'Neural Network Forces:'
!!        if(ldoforces)then
!!          do i1=1,num_atoms
!!            write(ounit,'(a9,i6,3f16.8,a8)')' NNforces ',i1,nntotalforce(1,i1),&
!!            nntotalforce(2,i1),nntotalforce(3,i1),' Ha/Bohr'
!!          enddo ! i1
!!        write(ounit,*)'-------------------------------------------------------------'
!!        endif !ldoforces'
!!        if(ldostress.and.lperiodic)then
!!          do i1=1,3
!!            write(ounit,'(a9,3f16.8,a10)')' NNstress ',nnstress(i1,1),&
!!            nnstress(i1,2),nnstress(i1,3),' Ha/Bohr^3'
!!          enddo ! i1
!!          write(ounit,*)'-------------------------------------------------------------'
!!          pressure=(nnstress(1,1)+nnstress(2,2)+nnstress(3,3))/3.0d0
!!          pressure=pressure*au2gpa
!!          write(ounit,'(a11,f16.8,a4)')' NNpressure ',pressure,' GPa'
!!          write(ounit,'(a11,f16.8,a5)')' NNpressure ',pressure*10.0d0,' kbar'
!!          write(ounit,*)'-------------------------------------------------------------'
!!        endif !ldostress'
!!!! sensitivity:
!!        if(lsens.and.lshort.and.(nn_type_short.eq.1))then
!!          do i1=1,nelem
!!            do i2=1,num_funcvalues_short_atomic(i1)
!!              write(ounit,'(a15,x,a2,x,i5,f16.8)')' NNsensitivity ',&
!!!! CHANGE ANDI: Mean square average is probably a better sensitivity value
!!               !element(i1),i2,sens(i1,i2)/dble(num_atoms)
!!                element(i1),i2,sqrt(sens(i1,i2)/dble(num_atoms_element(i1)))
!!!! END CHANGE
!!            enddo ! i2
!!            write(ounit,*)'-------------------------------------------------------------'
!!          enddo ! i1
!!        endif ! lsens
!!        if(lsens.and.(nn_type_short.eq.2))then
!!          counter = 0
!!          do i1=1,nelem
!!           do i2=i1,nelem
!!            counter = counter + 1
!!            do i3=1,num_funcvalues_short_pair(counter)
!!              write(ounit,'(a15,x,a2,x,a2,x,i5,f16.8)')' NNsensitivity ',&
!!!! CHANGE ANDI: Mean square average is probably a better sensitivity value.
!!!!              What is the correct normalization for nn_type_short 2 (num_atoms vs. num_atoms_element(i1))?
!!               !element(i1),element(i2),i3,sens(counter,i3)/dble(num_atoms)
!!                element(i1),element(i2),i3,sqrt(sens(counter,i3)/dble(num_atoms))
!!!! END CHANGE
!!            enddo ! i3
!!            write(ounit,*)'-------------------------------------------------------------'
!!          enddo ! i2
!!         enddo ! i1
!!        endif ! lsens
!!      endif ! mpirank.eq.0
!!
!!======================================================================
!! write structures in various formats to files for visualization for process 0
!!======================================================================
!!!! xyz file 
!!      if((mpirank.eq.0).and.(.not.lmd))then
!!        if(lwritexyz)then
!!          call writexyz(num_atoms,&
!!            xyzstruct,elementsymbol)
!!        endif ! lwritexyz
!!!! pdb file
!!        if(lwritepdb)then
!!          call writepdb(num_atoms,&
!!            lattice,xyzstruct,elementsymbol,lperiodic)
!!        endif ! lwritepdb
!!
!!!! povray file
!!        if(lwritepov)then 
!!          call writepov(num_atoms,zelem,&
!!            lattice,xyzstruct,lperiodic)
!!        endif ! lwritepov
!!!! pwscf input file
!!        if(lwritepw)then
!!          call writepwscf(num_atoms,&
!!            lattice,xyzstruct,elementsymbol,lperiodic)
!!        endif ! lwritepw
!!      endif ! mpirank.eq.0
!!!!
!!!! calculate final output timing
!!      if(lfinetime)then
!!        call abstime(timeoutputend,dayoutput)
!!        timeoutput=timeoutputend-timeoutputstart
!!      endif ! lfinetime
!!!!

!!======================================================================
!!SK: Now pass all results to i-pi
!!SK: here the driver loop ends 
!!======================================================================

         pot=nntotalenergy 
!!         forces=RESHAPE(nntotalforce, (/ 3 * nat /) )   ! return force in atomic units
         forces=nntotalforce   ! return force in atomic units           
         vir= nnstress

!!SK: additional information can be passed here
!!         extra_string = trim(comm_string)   
!!         comm_string=' ' ! resets comm_string after it has been passed              

         hasdata=.true.
         isinit = .false. !SK: resets init so that we will get replicaindex again at next step!

      else if (trim(header)=="GETFORCE") then
         if (mpirank.eq.0) write(*,*) " @ DRIVER MODE: Returning v,forces,stress "
         if (mpirank.eq.0) then      
            call writebuffer(socket,"FORCEREADY  ",MSGLEN, info)
            call writebuffer(socket,pot,8, info)
            call writebuffer(socket,nat,4, info)            
            call writebuffer(socket,forces,3*nat*8, info)  !!SK: do we need combuf?

            call writebuffer(socket,vir,9*8, info)
!!            siz=len(trim(extra_string))
            call writebuffer(socket,siz,4, info)
!!            call writebuffer(socket,trim(extra_string),siz, info)
         endif
!        if (drv_check(info)) exit 
         hasdata=.false.
      endif
    ENDDO driver_loop    
    
!!======================================================================
!! deallocate arrays
!!======================================================================
      if(lshort.and.(nn_type_short.eq.1))then
        deallocate(sens) 
      elseif(lshort.and.(nn_type_short.eq.2))then
        deallocate(sens) 
      endif
      if(lelec.and.(nn_type_elec.eq.1).or.(nn_type_elec.eq.3).or.(nn_type_elec.eq.4))then
        deallocate(sense) 
      endif
!!
      call abstime(timemode3end,daymode3)
!!
      return
      end
