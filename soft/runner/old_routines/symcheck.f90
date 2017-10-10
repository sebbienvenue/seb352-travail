!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine symcheck(&
        nblock,nelem,nucelem,max_num_atoms,&
        elementindex,&
        maxnum_weightsshort,maxnum_weightsewald,&
        num_weightsshort,num_weightsewald,&
        maxnum_funcvalues,maxnum_funcvaluese,&
        num_funcvalues,num_funcvaluese,&
        symthres,chargethres,forcethres,energythres,&
        lshort,lewald,luseforces,lscalesym,lcentersym,&
        ldebug)
!!
      use mpi_mod
      use fileunits
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5   ! internal
      integer j1,j2,j3         ! internal
      integer nelem            ! in
      integer elementindex(102) ! in
      integer maxnum_weightsshort  ! in
      integer maxnum_weightsewald  ! in
      integer num_weightsshort(nelem)  ! in
      integer num_weightsewald(nelem)  ! in
      integer maxnum_funcvalues    ! in
      integer maxnum_funcvaluese   ! in
      integer num_funcvalues(nelem)    ! in
      integer num_funcvaluese(nelem)   ! in
      integer refstruct         ! internal
      integer refatom           ! internal
!!      integer, dimension(:) , allocatable :: zelemref
      integer, dimension(:) , allocatable :: natoms
      integer zelemref
      integer zelemtest
      integer nattemp           ! internal
      integer idummy            ! internal
      integer numvalues         ! internal
      integer nstruct           ! internal
      integer nucelem(nelem)    ! in
      integer nblock            ! in
      integer max_num_atoms     ! in
!!
      real*8, dimension(:) , allocatable :: symfunctionref 
      real*8, dimension(:) , allocatable :: symfunctiontest 
      real*8 rdummy        ! internal
      real*8 refcharge     ! internal
      real*8 refenergy     ! internal
      real*8 refforce(3)   ! internal
      real*8 refxyz(3)     ! internal
      real*8 testcharge    ! internal
      real*8 testforce(3)  ! internal
      real*8 deltasym      ! internal
      real*8 symthres      ! in
      real*8 chargethres   ! in
      real*8 forcethres    ! in
      real*8 energythres   ! in
      real*8 minvalue(nelem,maxnum_funcvalues)                 ! internal 
      real*8 maxvalue(nelem,maxnum_funcvalues)                 ! internal 
      real*8 avvalue(nelem,maxnum_funcvalues)                  ! internal 
      real*8 minvaluee(nelem,maxnum_funcvaluese)                 ! internal 
      real*8 maxvaluee(nelem,maxnum_funcvaluese)                 ! internal 
      real*8 avvaluee(nelem,maxnum_funcvaluese)                  ! internal 
!!
      logical lshort       ! in
      logical lewald       ! in
      logical luseforces   ! in
      logical ldebug       ! in
      logical lscalesym    ! in
      logical lcentersym   ! in
      logical lperiodic    ! internal
!!
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
      if(mpisize.gt.1)then
        write(ounit,*)'Error: mode 4 is not parallelized'
        stop
      endif
!!
!! initial checks of the input
      if(lshort.and.lewald)then
        write(ounit,*)'Error: both lshort and lewald are set'
        stop
      elseif((.not.lshort).and.(.not.lewald))then
        write(ounit,*)'Error: both lshort and lewald are NOT set'
        stop
      endif
      if(luseforces.and.(.not.lshort))then
        write(ounit,*)'Error: luseforces is set, but lshort is not set'
        stop
      endif
!!
      if(lewald)then
        write(ounit,*)'Checking symmetry functions for atomic charges'
      elseif(luseforces)then
        write(ounit,*)'Checking symmetry functions for atomic forces'
      elseif(lshort)then
        write(ounit,*)'Checking symmetry functions for energies'
      endif
!!
!! count the number of structures in the training set
!!-------------------------------------------------------
      nstruct=0
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)
 30   continue
      read(symunit,*,end=31) nattemp
      do i1=1,nattemp
        read(symunit,*)idummy
      enddo ! i1
      read(symunit,*)rdummy
      nstruct=nstruct+1
      goto 30
 31   continue
      close(symunit)
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Number of structures in training set: ',nstruct
!!'
!! count number of atoms in each structure
!!-------------------------------------------------------
      allocate(natoms(nstruct))
      natoms(:)=0
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)
      do i1=1,nstruct
        read(symunit,*)natoms(i1)
        do i2=1,natoms(i1)
          read(symunit,*)idummy
        enddo ! i2
        read(symunit,*)rdummy
      enddo ! i1
      close(symunit)
!!
!! initializations
!!-------------------------------------------------------
      if(lshort)then
        numvalues=num_funcvalues
      elseif(lewald)then
        numvalues=num_funcvaluese
      else
        write(ounit,*)'Error in symcheck'
        stop
      endif
      allocate(symfunctionref(numvalues))
      allocate(symfunctiontest(numvalues))
!!     symthres    = 0.01d0
!!      chargethres = 0.0001d0
!!      forcethres  = 0.0001d0
!!  
!! 
!! identify the reference structure and count the atoms
!!      refstruct = 1
!!      refatom   = 2
!!      write(ounit,*)'Please enter the structure and the atom '
!!      read(*,*)refstruct,refatom
!!      write(ounit,*)'structure ',refstruct,' atom ',refatom
!!
!!
!! write headers
!!-------------------------------------------------------
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,'(a,f14.6)')' symmetry function threshold: ',symthres
      if(lewald)then !'
        write(ounit,'(a,f14.6)')' atomic charge criterion:     ',chargethres
      elseif(luseforces)then
        write(ounit,'(a,f14.6)')' atomic force criterion:      ',forcethres
      else
        write(ounit,*)'Error in symcheck:'
        write(ounit,*)'Currently the symmetry function check is only '
        write(ounit,*)'implemented for atomic charges and atomic forces'
        stop
      endif
      write(ounit,*)'-------------------------------------------------------------'
      if(lewald)then
        write(ounit,'(3a)')'       refstruct   refatom   element',&
          '  structure   atom',&
          '      deltasym     refcharge    testcharge         Delta'
      elseif(luseforces)then
        write(ounit,'(3a)')'       refstruct   refatom   element',&
          '  structure   atom',&
          '      deltasym         refforce     testforce         Delta'
      else

      endif
      write(ounit,'(3a)')' --------------------------------------------',&
        '--------------------------------------------------',&
        '-------------------'
!!'
!! loop over all reference atoms
!!-------------------------------------------------------
      do j1=1,nstruct
        refstruct=j1
        do j2=1,natoms(j1)
          refatom=j2
!!          write(ounit,*)'Checking structure, atom ',refstruct,refatom
!!
!! open files for reading the reference data
!!-------------------------------------------------------
          if(lshort)then
            open(symunit,file='function.data',form='formatted',status='old')
            rewind(symunit)
            open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
            rewind(trainstructunit)
          endif ! lshort
          if(lewald)then
            open(symunit,file='functione.data',form='formatted',status='old')
            rewind(symunit)
            open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
            rewind(trainstructunit) !'
          endif ! lewald
!!
!! identify reference structure and read reference data
!!-------------------------------------------------------
            i1=0  ! counter for structures
 10         continue
            i1=i1+1
            read(symunit,*,END=11)nattemp 
            read(trainstructunit,*)idummy,lperiodic
            if(lperiodic)then
              do i3=1,3
                read(trainstructunit,*)rdummy
              enddo
            endif
            i2=0  ! counter for atoms
            if(refstruct.eq.i1)then
              do i2=1,nattemp
                if(i2.eq.refatom)then 
                  read(symunit,*,END=12)zelemref,(symfunctionref(i3),i3=1,numvalues)
                  if(lewald)then
                    read(trainstructunit,*)idummy,refxyz(1),refxyz(2),refxyz(3),refcharge
                  elseif(luseforces)then
                    read(trainstructunit,*)idummy,refxyz(1),refxyz(2),refxyz(3),refcharge,&
                      rdummy,refforce(1),refforce(2),refforce(3)
                  endif
                  goto 13
                else
                  read(symunit,*)idummy
                  read(trainstructunit,*)idummy
                endif
              enddo
            else
              do i2=1,nattemp
                read(symunit,*)idummy
                read(trainstructunit,*)idummy
              enddo
              read(symunit,*)rdummy
              goto 10
            endif
 11         continue
            write(ounit,*)'Error: number of requested reference structure is too large ',refstruct
            stop
 12         continue
            write(ounit,*)'Error: number of requested reference atom is too large ',refatom
            stop
 13         continue
!!
!! close the files for determining the reference atoms
!!--------------------------------------------------------
            if(lshort)then
              close(symunit)
              close(trainstructunit)
            endif ! lshort
            if(lewald)then
              close(symunit)
              close(trainstructunit)
            endif ! lewald
!!
!! write summary of the reference
!!--------------------------------------------------------
!!        write(ounit,*)'reference atom found:'
!!        write(ounit,'(i3,x,200f14.10)')zelemref,&
!!          (symfunctionref(i3),i3=1,numvalues)
!!        if(lewald)then
!!          write(ounit,'(a,f20.8)')' reference charge : ',refcharge
!!          write(ounit,'(a,3f20.8)')' reference xyz    : ',(refxyz(i3),i3=1,3)
!!        elseif(luseforces)then
!!          write(ounit,'(a,3f20.8)')' reference forces : ',(refforce(i3),i3=1,3)
!!          write(ounit,'(a,3f20.8)')' reference xyz    : ',(refxyz(i3),i3=1,3)
!!        else
!!          write(ounit,*)'reference energy : ',refenergy
!!        endif
!!        write(ounit,*)'-------------------------------------------------------------'
!!
!!            if(lewald)then
!!              write(ounit,'(a,3i5,f14.8)')&
!!                ' Checking ',refstruct,refatom,zelemref,refcharge
!!            elseif(luseforces)then
!!              write(ounit,'(a,3i5,3f14.8)')&
!!                ' Checking ',refstruct,refatom,zelemref,(refforce(i3),i3=1,3)
!!            else
!!              write(ounit,*)'not finished for lshort'
!!            endif
!!'
!!
!!
!! open files for symmetry function checks
!!--------------------------------------------------------
            if(lshort)then
              open(symunit,file='function.data',form='formatted',status='old')
              rewind(symunit)
              open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
              rewind(trainstructunit)
            endif ! lshort
            if(lewald)then
              open(symunit,file='functione.data',form='formatted',status='old')
              rewind(symunit)
              open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
              rewind(trainstructunit)
            endif ! lewald
!!
!! read information for comparisons 
!!--------------------------------------------------------
            i1=0  ! counter for structures
 21         continue
            i1=i1+1
            read(symunit,*,END=22)nattemp
            read(trainstructunit,*)idummy,lperiodic
            if(lperiodic)then
              do i3=1,3
                read(trainstructunit,*)rdummy
              enddo
            endif
            do i2=1,nattemp
              read(symunit,*)zelemtest,(symfunctiontest(i3),i3=1,numvalues)
              read(trainstructunit,*)idummy,rdummy,rdummy,rdummy,testcharge,&
                rdummy,testforce(1),testforce(2),testforce(3)
              if(lewald)then
                if((zelemref.eq.zelemtest).and.(.not.((refstruct.eq.i1).and.(refatom.eq.i2))))then
!! the deltasym criterion is probably too basic
                  deltasym=0.0d0
                  do i3=1,numvalues
                    deltasym=deltasym+(symfunctionref(i3)-symfunctiontest(i3))**2
                  enddo
                  deltasym=dsqrt(deltasym)
                  if((deltasym.lt.symthres).and.(abs(refcharge-testcharge).gt.chargethres))then
                    write(ounit,'(a,3i10,i11,i7,f14.8,3f14.8)')&
                      'FOUND ',refstruct,refatom,zelemref,i1,i2,deltasym,&
                      refcharge,testcharge,refcharge-testcharge
!!                    write(ounit,*)'Problematic structure found:'
!!                    do i3=1,numvalues
!!                      write(ounit,'(i6,2f20.8)')i3,symfunctionref(i3),symfunctiontest(i3)
!!                    enddo ! i3
                  endif
                endif ! zelemref.eq.zelemtest
              elseif(luseforces)then
                if((zelemref.eq.zelemtest).and.(.not.((refstruct.eq.i1).and.(refatom.eq.i2))))then
!! the deltasym criterion is probably too basic
                  deltasym=0.0d0
                  do i3=1,numvalues
                    deltasym=deltasym+(symfunctionref(i3)-symfunctiontest(i3))**2
                  enddo
                  deltasym=dsqrt(deltasym)
                  if((deltasym.lt.symthres).and.&
                    (abs(refforce(1)-testforce(1)).gt.forcethres))then !.or.&
!!                     (abs(refforce(2)-testforce(2)).gt.forcethres).or.&
!!                     (abs(refforce(3)-testforce(3)).gt.forcethres)))then    
                    write(ounit,'(a,3i10,i11,i7,f14.8,a3,3f14.8)')&
                      'FOUND ',refstruct,refatom,zelemref,i1,i2,deltasym,&
                      ' X ',refforce(1),testforce(1),refforce(1)-testforce(1)
                  endif
                  if((deltasym.lt.symthres).and.&
                    (abs(refforce(2)-testforce(2)).gt.forcethres))then
                    write(ounit,'(a,3i10,i11,i7,f14.8,a3,3f14.8)')&
                      'FOUND ',refstruct,refatom,zelemref,i1,i2,deltasym,&
                      ' Y ',refforce(2),testforce(2),refforce(2)-testforce(2)
                  endif
                  if((deltasym.lt.symthres).and.&
                    (abs(refforce(3)-testforce(3)).gt.forcethres))then
                    write(ounit,'(a,3i10,i11,i7,f14.8,a3,3f14.8)')&
                      'FOUND ',refstruct,refatom,zelemref,i1,i2,deltasym,&
                      ' Z ',refforce(3),testforce(3),refforce(3)-testforce(3)
                  endif
                endif ! zelemref.eq.zelemtest
              else

              endif
            enddo ! i2
            read(symunit,*)rdummy
            goto 21
 22         continue
!!
!! close files
!!--------------------------------------------------------
            if(lshort)then
              close(symunit)
              close(trainstructunit)
            endif ! lshort
            if(lewald)then
              close(symunit)
              close(trainstructunit)
            endif ! lewald
!!
        enddo ! j2
      enddo ! j1
!!
!! final cleanup 
!!--------------------------------------------------------
      deallocate(symfunctionref)
      deallocate(symfunctiontest)
      deallocate(natoms)
!!
      return
      end
