!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine for the atomic case 

!! called by:
!! - fitting_atomic.f90
!! - fitting_short_atomic.f90
!!
      subroutine findcontradictions(ntrain,maxnum_funcvalues_local,&
        num_funcvalues_local,&
        minvalue_local,maxvalue_local,avvalue_local,&
        scmin_local,scmax_local)
!!
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer ntrain                                      ! in
      integer npoints1,npoints2                           ! internal
      integer ncount1,ncount2                             ! internal
      integer nstart1,nstart2                             ! internal
      integer nend1,nend2                                 ! internal
      integer nref                                        ! internal
      integer nstruct1,nstruct2                           ! internal
      integer idummy                                      ! internal
      integer i1,i2,i3,i4,i5                              ! internal
      integer num_atoms_local                             ! internal      
      integer ncontra(nelem)                              ! internal
      integer maxnum_funcvalues_local                     ! in
      integer zelem_local(max_num_atoms)                  ! internal
      integer num_funcvalues_local(nelem)                 ! in

      real*8 zdummy                                       ! internal
      real*8 zdummy1                                      ! internal
      real*8 zdummy2                                      ! internal
      real*8 deltaf                                       ! internal
      real*8 deltaf1                                      ! internal
      real*8 deltaf2                                      ! internal
      real*8 deltag                                       ! internal
      real*8 deltagmax(nelem)                             ! internal
      real*8 scmin_local                                  ! in
      real*8 scmax_local                                  ! in
      real*8 symfunction_local(maxnum_funcvalues_local,max_num_atoms)     ! internal
      real*8 symfunction1(maxnum_funcvalues_local,max_num_atoms,nblock)   ! internal
      real*8 symfunction2(maxnum_funcvalues_local,max_num_atoms,nblock)   ! internal
      real*8 minvalue_local(nelem,maxnum_funcvalues_local)                ! in
      real*8 maxvalue_local(nelem,maxnum_funcvalues_local)                ! in
      real*8 avvalue_local(nelem,maxnum_funcvalues_local)                 ! in


      logical lperiodic_local                             ! internal

!! arrays for first block of points
      integer     num_atoms_list1(nblock)
      integer     zelem_list1(nblock,max_num_atoms)
      real*8      lattice_list1(3,3,nblock)
      real*8      xyzstruct_list1(3,max_num_atoms,nblock)
      real*8      totalforce_list1(3,max_num_atoms,nblock)
      logical     lperiodic_list1(nblock)
      integer     index1(nblock)
!!
!! arrays for second block of points
      integer     num_atoms_list2(nblock)
      integer     zelem_list2(nblock,max_num_atoms)
      real*8      lattice_list2(3,3,nblock)
      real*8      xyzstruct_list2(3,max_num_atoms,nblock)
      real*8      totalforce_list2(3,max_num_atoms,nblock)
      logical     lperiodic_list2(nblock)
      integer     index2(nblock)

!!
!! initializations
      nstart1   =1
      nstruct1  =0
      nstruct2  =0
      ncontra(:)=0
      deltagmax(:)=0.0d0
!!      deltagthres=0.002d0
!!      deltafthres=0.0001d0
!!
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,'(a)')' Searching for contradictory forces: '
!'
      write(ounit,'(a,f14.6,a,f14.6)')' Criteria: Delta F > ',deltafthres,' and Delta G < ',deltagthres
      write(ounit,'(a)')'                    structure    atom  structure    atom        deltaG        deltaF'
!!'
!!================================================================
!! read reference block
!!================================================================
      ncount1=ntrain
 11   continue
      nstruct1=0
      nstart2=1
      if(ncount1.gt.nblock)then
        npoints1=nblock
        ncount1=ncount1-nblock
        nend1 =nstart1+npoints1-1
      else
        npoints1=ncount1
        ncount1=ncount1-npoints1
        nend1 =nstart1+npoints1-1
      endif
!!      write(ounit,*)'reading reference structures ',nstart1,nend1,npoints1

      nref=0
      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      rewind(trainstructunit)
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)

 20   continue
      nstruct1=nstruct1+1
!! read previous dummy structures which are not used
      if(nstruct1.gt.ntrain)goto 21
      if(nstruct1.lt.nstart1)then ! do not use these structures 
        read(symunit,*)num_atoms_local
        do i1=1,num_atoms_local
          read(symunit,*)idummy
        enddo
        read(symunit,*)zdummy     
        read(trainstructunit,*)idummy,lperiodic_local
        if(lperiodic_local)then
          do i1=1,3
            read(trainstructunit,*)zdummy
          enddo
        endif
        do i1=1,num_atoms_local
          read(trainstructunit,*)idummy
        enddo
!! if we have completed reading relevant structures
      elseif(nstruct1.gt.nend1)then  
        nstruct1=nstruct1-1
        goto 21
!! read relevant set of reference structures
      else 
!!        write(ounit,*)' reading nstruct1 ',nstruct1
        nref=nref+1
        index1(nref)=nstruct1
!!        write(ounit,*)'structure found ',nref,index1(nref)
        read(symunit,*)num_atoms_list1(nref)
        do i1=1,num_atoms_list1(nref)
          read(symunit,*)idummy
          backspace(symunit)
          read(symunit,*)zelem_list1(nref,i1),(symfunction1(i2,i1,nref),&
            i2=1,num_funcvalues_local(elementindex(idummy)))
        enddo
        read(symunit,*)zdummy             ! energy line in function.data
!! scale the symmetry functions if requested
        zelem_local(:)=zelem_list1(nref,:)
        do i1=1,num_atoms_list1(nref)
          do i2=1,num_funcvalues_local(elementindex(zelem_list1(nref,i1)))
            symfunction_local(i2,i1)=symfunction1(i2,i1,nref)
          enddo
        enddo
        call scalesymone(nelem,&
          maxnum_funcvalues_local,num_funcvalues_local,num_atoms_list1(nref),&
          zelem_local,symfunction_local,&
          minvalue_local,maxvalue_local,avvalue_local,&
          scmin_local,scmax_local)
        symfunction1(:,:,nref)=symfunction_local(:,:)

!! read the forces
        read(trainstructunit,*)idummy,lperiodic_list1(nref)
        if(lperiodic_list1(nref))then
          do i1=1,3
            read(trainstructunit,*)(lattice_list1(i1,i2,nref),i2=1,3)
          enddo
        endif
        do i1=1,num_atoms_list1(nref)
          read(trainstructunit,*)zdummy,(xyzstruct_list1(i2,i1,nref),i2=1,3),&
            zdummy1,zdummy2,(totalforce_list1(i2,i1,nref),i2=1,3)
        enddo
      endif
      goto 20
 21   continue
      close(trainstructunit)
      close(symunit)
      nstart1=nstart1+npoints1
!!================================================================
!! reading reference block done
!!================================================================

!!================================================================
!! read second block for comparison
!!================================================================
      ncount2=ntrain
 31   continue
      nstruct2=0
      if(ncount2.gt.nblock)then
        npoints2=nblock
        ncount2=ncount2-nblock
        nend2 =nstart2+npoints2-1
      else
        npoints2=ncount2
        ncount2=ncount2-npoints2
        nend2 =nstart2+npoints2-1
      endif
!!      write(ounit,*)'reading comparison structures ',nstart2,nend2,npoints2

      nref=0
      open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
      rewind(trainstructunit)
      open(symunit,file='function.data',form='formatted',status='old')
      rewind(symunit)

 40   continue
      nstruct2=nstruct2+1
!! read previous dummy structures which are not used
      if(nstruct2.gt.ntrain)goto 41
      if(nstruct2.lt.nstart2)then
        read(symunit,*)num_atoms_local
        do i1=1,num_atoms_local
          read(symunit,*)idummy
        enddo
        read(symunit,*)zdummy
        read(trainstructunit,*)idummy,lperiodic_local
        if(lperiodic_local)then
          do i1=1,3
            read(trainstructunit,*)zdummy
          enddo
        endif
        do i1=1,num_atoms_local
          read(trainstructunit,*)idummy
        enddo
!! if we have completed reading relevant structures
      elseif(nstruct2.gt.nend2)then
        nstruct2=nstruct2-1
        goto 41
!! read relevant set of reference structures
      else
!!        write(ounit,*)'reading nstruct2 ',nstruct2
        nref=nref+1
        index2(nref)=nstruct2
!!        write(ounit,*)'structure found ',nref,index2(nref)
        read(symunit,*)num_atoms_list2(nref)
        do i1=1,num_atoms_list2(nref)
          read(symunit,*)idummy
          backspace(symunit)
          read(symunit,*)zelem_list2(nref,i1),(symfunction2(i2,i1,nref),&
            i2=1,num_funcvalues_local(elementindex(idummy)))
        enddo
        read(symunit,*)zdummy             ! energy line in function.data
!! scale the symmetry functions if requested
        zelem_local(:)=zelem_list2(nref,:)
        symfunction_local(:,:)=symfunction2(:,:,nref)
        call scalesymone(nelem,&
          maxnum_funcvalues_local,num_funcvalues_local,&
          num_atoms_list2(nref),&
          zelem_local,symfunction_local,&
          minvalue_local,maxvalue_local,avvalue_local,&
          scmin_local,scmax_local)
        symfunction2(:,:,nref)=symfunction_local(:,:)

!! read the forces
        read(trainstructunit,*)idummy,lperiodic_list2(nref)
        if(lperiodic_list2(nref))then
          do i1=1,3
            read(trainstructunit,*)(lattice_list2(i1,i2,nref),i2=1,3)
          enddo
        endif
        do i1=1,num_atoms_list2(nref)
          read(trainstructunit,*)zdummy,(xyzstruct_list2(i2,i1,nref),i2=1,3),&
            zdummy1,zdummy2,(totalforce_list2(i2,i1,nref),i2=1,3)
        enddo
      endif
      goto 40
 41   continue
      close(trainstructunit)
      close(symunit)
      nstart2=nstart2+npoints2
!!================================================================
!! reading comparison block done
!!================================================================

!!================================================================
!! compare now all atoms of all structures in both blocks
!!================================================================
!!      write(ounit,*)'now comparing ',npoints1,npoints2
      do i1=1,npoints1
        do i2=1,npoints2
!!          write(ounit,*)' structure pair ',index1(i1),index2(i2)
!! avoid double counting of structure pairs
          if(index1(i1).gt.index2(i2))then
!! compare all atoms of all structures
            do i3=1,num_atoms_list1(i1)
              do i4=1,num_atoms_list2(i2)
!! compare environments only if atoms have the same element
                if(zelem_list1(i1,i3).eq.zelem_list2(i2,i4))then
                  deltaf1=(totalforce_list1(1,i3,i1))**2&
                    +(totalforce_list1(2,i3,i1))**2+(totalforce_list1(3,i3,i1))**2
                  deltaf1=dsqrt(deltaf1)
                  deltaf2=(totalforce_list2(1,i4,i2))**2&
                    +(totalforce_list2(2,i4,i2))**2+(totalforce_list2(3,i4,i2))**2
                  deltaf2=dsqrt(deltaf2)
                  deltaf=abs(deltaf1-deltaf2)
                  deltag = 0.0d0
                  do i5=1,num_funcvalues_local(elementindex(zelem_list1(i1,i3)))
                    deltag = deltag + (symfunction1(i5,i3,i1)-symfunction2(i5,i4,i2))**2 
                  enddo ! i5
                  deltag = dsqrt(deltag)
                  deltagmax(elementindex(zelem_list1(i1,i3)))&
                    =max(deltag,deltagmax(elementindex(zelem_list1(i1,i3))))
!! write output
                  if((deltag.lt.deltagthres).and.(deltaf.gt.deltafthres))then
                    ncontra(elementindex(zelem_list1(i1,i3)))=ncontra(elementindex(zelem_list1(i1,i3)))+1
!                    write(ounit,'(a,a2,2f14.6)')' CONTRADICTION ',&
!                      element(elementindex(zelem_list1(i1,i3))),&
!                      deltag,deltaf
                    write(ounit,'(a,a2,4x,2i8,3x,2i8,2f14.6)')' CONTRADICTION ',&
                      element(elementindex(zelem_list1(i1,i3))),& !'
                      index1(i1),i3,index2(i2),i4,&
                      deltag,deltaf
                  endif
                endif
              enddo ! i4
            enddo ! i3
          elseif(index1(i1).eq.index2(i2))then
!! compare all atoms of all structures
            do i3=1,num_atoms_list1(i1)-1
              do i4=i3+1,num_atoms_list2(i2)
!! compare environments only if atoms have the same element
                if(zelem_list1(i1,i3).eq.zelem_list2(i2,i4))then

!! this is wrong because only the absolute values of the forces can be compared
                  deltaf=(totalforce_list1(1,i3,i1)-totalforce_list2(1,i4,i2))**2 &
                       + (totalforce_list1(2,i3,i1)-totalforce_list2(2,i4,i2))**2 &
                       + (totalforce_list1(3,i3,i1)-totalforce_list2(3,i4,i2))**2
                  deltaf1=(totalforce_list1(1,i3,i1))**2+(totalforce_list1(2,i3,i1))**2+(totalforce_list1(3,i3,i1))**2
                  deltaf1=dsqrt(deltaf1)
                  deltaf2=(totalforce_list2(1,i4,i2))**2+(totalforce_list2(2,i4,i2))**2+(totalforce_list2(3,i4,i2))**2
                  deltaf2=dsqrt(deltaf2)
                  deltaf=abs(deltaf1-deltaf2)
!                  write(ounit,'(3f14.8)')totalforce_list1(1,i3,i1),totalforce_list2(1,i4,i2),totalforce_list1(1,i3,i1)-totalforce_list2(1,i4,i2) 
!                  write(ounit,'(3f14.8)')totalforce_list1(2,i3,i1),totalforce_list2(2,i4,i2),totalforce_list1(2,i3,i1)-totalforce_list2(2,i4,i2) 
!                  write(ounit,'(3f14.8)')totalforce_list1(3,i3,i1),totalforce_list2(3,i4,i2),totalforce_list1(3,i3,i1)-totalforce_list2(3,i4,i2)
                  deltag = 0.0d0
                  do i5=1,num_funcvalues_local(elementindex(zelem_list1(i1,i3)))
                    deltag = deltag + (symfunction1(i5,i3,i1)-symfunction2(i5,i4,i2))**2 
                  enddo ! i5
                  deltag = dsqrt(deltag)
                  deltagmax(elementindex(zelem_list1(i1,i3)))&
                    =max(deltag,deltagmax(elementindex(zelem_list1(i1,i3))))
!! write output
                  if((deltag.lt.deltagthres).and.(deltaf.gt.deltafthres))then
                    ncontra(elementindex(zelem_list1(i1,i3)))=ncontra(elementindex(zelem_list1(i1,i3)))+1
                    write(ounit,'(a,a2,4x,2i8,3x,2i8,2f14.6)')' CONTRADICTION ',&
                      element(elementindex(zelem_list1(i1,i3))),&
                      index1(i1),i3,index2(i2),i4,&
                      deltag,deltaf
                  endif
                endif
              enddo ! i4
            enddo ! i3
          endif
        enddo ! i2
      enddo ! i1
!!================================================================
!! comparison done 
!!================================================================

!! read next block of comparison points
      if(ncount2.gt.0) goto 31 ! do next block of comparison points
!!
!! read next block of reference points
      if(ncount1.gt.0) goto 11 ! do next block of reference points
!!
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)' Total contraditions found:'
      do i1=1,nelem
        write(ounit,'(a2,x,i5)')element(i1),ncontra(i1)
      enddo
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)' find_contradictions: Maximum distance of atoms in symmetry function space:'
      do i1=1,nelem
        write(ounit,'(a2,x,f14.10)')element(i1),deltagmax(i1)
      enddo
      write(ounit,*)'-------------------------------------------------------------'
!!
      return
      end
