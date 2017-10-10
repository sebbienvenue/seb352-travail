!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine environmentanalysis(iswitch,maxcutoff_local)
!!
      use fileunits
      use nnflags
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer i1,i2,i3,i4,i5                                      ! internal
      integer dim1,dim2,dim3,dim4                                 ! internal
      integer, dimension(:,:,:,:,:), allocatable :: temparray     ! internal 
      integer iswitch                                             ! in
      integer memsize                                             ! internal
      integer num_atoms                                           ! internal
      integer memthres                                            ! internal
      integer unit1                                               ! internal
      integer unit2                                               ! internal
      integer nstruct                                             ! internal
      integer, dimension(:), allocatable :: zelem                 ! internal 
      integer num_atomselement(nelem)                             ! internal
      integer idummy                                              ! internal
      integer jdummy                                              ! internal
      integer, allocatable :: lsta(:,:)                           ! numbers of neighbors
      integer, allocatable :: lstc(:)                             ! identification of atom
      integer, allocatable :: lste(:)                             ! nuclear charge of atom
      integer counter(4)                                          ! counter for number of atoms of each element, dimensioned for nelem<=4
      integer max_num_neighbors_local                             ! internal
      integer max_num_neighbors_temp                              ! internal
!!
      real*8 maxcutoff_local                                      ! in
      real*8 zdummy                                               ! internal
      real*8 zdummy1                                              ! internal
      real*8 zdummy2                                              ! internal
      real*8 lattice_local(3,3)                                   ! internal
      real*8 xyz_local(3,max_num_atoms)                           ! internal
      real*8 fxyz_local(3,max_num_atoms)                           ! internal
      real*8, allocatable :: lstb(:,:)                            ! xyz and r_ij
      real*8 symfunction_local(maxnum_funcvalues_short_atomic,max_num_atoms)   ! internal
      real*8 force                                                ! internal
!!
      logical lperiodic                                           ! internal
!!
      nstruct = 0
      allocate(lsta(2,max_num_atoms),lstc(listdim),lste(listdim),lstb(listdim,4))
!!
!! for dimensioning temparray we need to know max_num_neighbors_local over ALL structures
!! loop over all structures, determine for each max_num_neighbors_local and take the largest one
      if(.not.lenforcemaxnumneighborsatomic)then
        max_num_neighbors_local=0
        if(iswitch.eq.0)then
          open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
          rewind(trainstructunit)
          unit1=trainstructunit
          if(lshort)then
            unit2=symunit
            open(unit2,file='function.data',form='formatted',status='old')
            rewind(unit2)
          elseif(lelec.and.(nn_type_elec.eq.1))then
            unit2=symeunit
            open(unit2,file='functione.data',form='formatted',status='old')
            rewind(unit2)
          endif
        elseif(iswitch.eq.1)then
          open(teststructunit,file='teststruct.data',form='formatted',status='old')
          rewind(teststructunit)
          unit1=teststructunit
          if(lshort)then
            unit2=tymunit
            open(unit2,file='testing.data',form='formatted',status='old')
            rewind(unit2)
          elseif(lelec.and.(nn_type_elec.eq.1))then
            unit2=tymeunit
            open(unit2,file='testinge.data',form='formatted',status='old')
            rewind(unit2)
          endif
        endif
 10     continue
        if((lshort.and.(nn_type_short.eq.1)).or.(lelec.and.(nn_type_elec.eq.1)))then
          read(unit2,*,END=20)num_atoms
          nstruct=nstruct+1
          num_atomselement(:)=0
          do i1=1,num_atoms
            read(unit2,*)idummy
            backspace(unit2)
            read(unit2,*)jdummy,(symfunction_local(i2,i1),i2=1,num_funcvalues_short_atomic(elementindex(idummy)))
          enddo
        elseif(nn_type_short.eq.2)then
          write(ounit,*)'ERROR: environment_analysis is not implemented for nn_type_short 2'
          stop !'
        endif
        read(unit2,*)zdummy             ! energy line in function.data
        read(unit1,*)idummy,lperiodic
        allocate(zelem(num_atoms))
        if(lperiodic)then
          do i1=1,3
            read(unit1,*)lattice_local(i1,1),lattice_local(i1,2),lattice_local(i1,3)
          enddo
        endif
        do i1=1,num_atoms
          read(unit1,*)zelem(i1),xyz_local(1,i1),xyz_local(2,i1),xyz_local(3,i1),&
            zdummy1,zdummy2,fxyz_local(1,i1),fxyz_local(2,i1),fxyz_local(3,i1)
        enddo
!!
!! get max_num_neighbors_local for this structure
        call getmaxnumneighbors(1,num_atoms,&
          num_atoms,max_num_atoms,max_num_neighbors_temp,&
          maxcutoff_local,lattice_local,xyz_local,lperiodic)
        max_num_neighbors_local=max(max_num_neighbors_local,max_num_neighbors_temp)
!!
        deallocate(zelem)
        goto 10
 20     continue
        close(unit2)
        close(unit1)
      else ! lenforcemaxnumneighborsatomic
        max_num_neighbors_local=max_num_neighbors_atomic_input
      endif
!!
!! FIXME: We need a more clever way to dimension the array temparray
!! It needs too much memory, and not all combinations of element numbers will be present
!! (maybe introduce pointer array?)
!!
!! JB: this is wrong because in periodic structures an atom can have more neighbors than there are atoms in the structure
!!      dim1=max_num_atoms
!!      dim2=max_num_atoms
!!      dim3=max_num_atoms
!!      dim4=max_num_atoms
      dim1=max_num_neighbors_local
      dim2=max_num_neighbors_local
      dim3=max_num_neighbors_local
      dim4=max_num_neighbors_local
      if(nelem.lt.4)dim4=0
      if(nelem.lt.3)dim3=0
      if(nelem.lt.2)dim2=0
!!
!! check if array is too large
      memsize=(dim1+1)*(dim2+1)*(dim3+1)*(dim4+1)*2
      memthres=1024*1024*1024 ! roughly 1 GB
      if(memsize.gt.memthres)then
        write(ounit,'(2a,i20)')'### WARNING ### : environmentanalysis.f90:',&
          ' Analysis of input.data cannot be done (insufficient memory) ',memsize
        return !'
      endif
!!
!! allocate array for analysis
      allocate(temparray(0:dim1,0:dim2,0:dim3,0:dim4,nelem))
      temparray(:,:,:,:,:)=0
!!
!! now analyze all environments in all structures
      if(iswitch.eq.0)then
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
        unit1=trainstructunit
        if(lshort)then
          unit2=symunit
          open(unit2,file='function.data',form='formatted',status='old')
          rewind(unit2)
        elseif(lelec.and.(nn_type_elec.eq.1))then
          unit2=symeunit
          open(unit2,file='functione.data',form='formatted',status='old')
          rewind(unit2)
        endif
      elseif(iswitch.eq.1)then
        open(teststructunit,file='teststruct.data',form='formatted',status='old')
        rewind(teststructunit)
        unit1=teststructunit
        if(lshort)then
          unit2=tymunit
          open(unit2,file='testing.data',form='formatted',status='old')
          rewind(unit2)
        elseif(lelec.and.(nn_type_elec.eq.1))then
          unit2=tymeunit
          open(unit2,file='testinge.data',form='formatted',status='old')
          rewind(unit2)
        endif
      endif
 30   continue
      if((lshort.and.(nn_type_short.eq.1)).or.(lelec.and.(nn_type_elec.eq.1)))then
        read(unit2,*,END=40)num_atoms
        nstruct=nstruct+1
        num_atomselement(:)=0
        do i1=1,num_atoms
          read(unit2,*)idummy
          backspace(unit2)
          read(unit2,*)jdummy,(symfunction_local(i2,i1),i2=1,num_funcvalues_short_atomic(elementindex(idummy)))
        enddo
      elseif(nn_type_short.eq.2)then
        write(ounit,*)'ERROR: environment_analysis is not implemented for nn_type_short 2'
        stop !'
      endif
      read(unit2,*)zdummy             ! energy line in function.data
      read(unit1,*)idummy,lperiodic
      allocate(zelem(num_atoms))
      if(lperiodic)then
        do i1=1,3
          read(unit1,*)lattice_local(i1,1),lattice_local(i1,2),lattice_local(i1,3)
        enddo
      endif
      do i1=1,num_atoms
        read(unit1,*)zelem(i1),xyz_local(1,i1),xyz_local(2,i1),xyz_local(3,i1),&
          zdummy1,zdummy2,fxyz_local(1,i1),fxyz_local(2,i1),fxyz_local(3,i1)
      enddo

!! analyze environments
!! get neighbor list here for atoms 1 to num_atoms
      call neighbor_para(1,num_atoms,&
        num_atoms,zelem,lsta,lstb,lstc,lste,&
        maxcutoff_local,lattice_local,xyz_local,lperiodic)
!!
      do i1=1,num_atoms
        if(num_atoms.gt.2)then
          if(lsta(1,i1).eq.lsta(2,i1))then
            if(iswitch.eq.0)then
              write(ounit,'(a,i8,a,i8,a)')'### WARNING ### training structure ',&
                nstruct,' atom ',i1,' has just one neighbor but is no dimer'
              force=0.0d0
              force=fxyz_local(1,i1)**2 + fxyz_local(2,i1)**2 + fxyz_local(3,i1)**2
              force=dsqrt(force)
              write(ounit,'(a,2i8,f14.8,100f14.6)')'WARNING ',nstruct,i1,force,&
                (symfunction_local(i2,i1),i2=1,num_funcvalues_short_atomic(elementindex(zelem(i1))))
            endif
          endif
        endif
        if(lsta(1,i1).gt.lsta(2,i1))then
          if(iswitch.eq.0)then
            write(ounit,'(a,i8,a,i8,a)')'### WARNING ### training structure ',&
              nstruct,' has a non-bonded atom ',i1
            force=0.0d0
            force=fxyz_local(1,i1)**2 + fxyz_local(2,i1)**2 + fxyz_local(3,i1)**2
            force=dsqrt(force)
            write(ounit,'(a,2i8,f14.8,100f14.6)')'WARNING ',nstruct,i1,force,&
            (symfunction_local(i2,i1),i2=1,num_funcvalues_short_atomic(elementindex(zelem(i1))))
          endif !'
        endif
      enddo ! i1



      do i1=1,num_atoms
        counter(:)=0
        do i2=lsta(1,i1),lsta(2,i1) ! loop over neighbors 
          counter(elementindex(lste(i2)))=counter(elementindex(lste(i2)))+1
        enddo ! i2
        temparray(counter(1),counter(2),counter(3),counter(4),elementindex(zelem(i1)))&
          =temparray(counter(1),counter(2),counter(3),counter(4),elementindex(zelem(i1)))+1
      enddo ! i1
!!
!! next structure
      deallocate(zelem)
      goto 30
 40   continue
      close(unit2)
      close(unit1)
!!
!! output
      do i5=1,nelem
        write(ounit,*)'-------------------------------------------------------------'
        if(iswitch.eq.0)then
          write(ounit,'(a,a2)')' Training atomic environments of ',element(i5) !'
        elseif(iswitch.eq.1)then
          write(ounit,'(a,a2)')' Testing atomic environments of ',element(i5) !'
        endif
        write(ounit,'(18x,4(a2,7x))')(element(i1),i1=1,nelem)
        do i1=0,dim1
          do i2=0,dim2
            do i3=0,dim3
              do i4=0,dim4
                if((i1+i2+i3+i4).gt.0)then ! avoid output for structures with no elements
                  if(temparray(i1,i2,i3,i4,i5).gt.0)then ! print only existing clusters
                    if(nelem.eq.1)then
                      write(ounit,'(a,a2,i8,16x,i20)')&
                        ' CLUSTER ',element(i5),i1,temparray(i1,i2,i3,i4,i5)
                    elseif(nelem.eq.2)then
                      write(ounit,'(a,a2,i8,x,i8,11x,i20)')&
                        ' CLUSTER ',element(i5),i1,i2,temparray(i1,i2,i3,i4,i5)
                    elseif(nelem.eq.3)then
                      write(ounit,'(a,a2,i8,x,i8,x,i8,6x,i20)')&
                        ' CLUSTER ',element(i5),i1,i2,i3,temparray(i1,i2,i3,i4,i5)
                    elseif(nelem.eq.4)then
                      write(ounit,'(a,a2,i8,x,i8,x,i8,x,i8,x,i20)')&
                        ' CLUSTER ',element(i5),i1,i2,i3,i4,temparray(i1,i2,i3,i4,i5)
                    else
                      write(ounit,*)'ERROR in analyzeinput.f90, valid only up to 4 elements'
                      stop
                    endif
                  endif
                endif
              enddo
            enddo
          enddo
        enddo 
      enddo
!!
      deallocate(temparray)
      deallocate(lsta,lstc,lste,lstb)
!!
      return
      end
