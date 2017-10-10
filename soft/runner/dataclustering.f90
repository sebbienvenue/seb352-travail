!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: Analyze the atomic environments in symmetry function space and group them into clusters, i.e. groups of atomic environments
!! This might help to remove redundant information from the training set by using only a few similar points per cluster

!! called by:
!!
      subroutine dataclustering(iswitch,nstruct_local,natoms_local,atomselem_local,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        minvalue_local,maxvalue_local,avvalue_local,&
        scmin_local,scmax_local)
!!
      use fileunits
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer iswitch                                           ! switch for training (0) or test (1) set 
      integer nstruct_local                                     ! total number of training or test points 
      integer num_clusters(nelem)                               ! internal
      integer nassigned(nelem)                                  ! internal
      integer npass                                             ! internal
      integer cluster_atoms(nelem,1000)                         ! internal (number of atoms in each cluster, dimension hard-coded!)
      integer natoms_local                                      ! in (total number of atoms in training or test set)
      integer atomselem_local(nelem)                            ! in (number of atoms of each element in training or test set)
      integer clusterindex(nstruct_local,max_num_atoms)         ! internal (label for cluster identity number)
      integer nstruct                                           ! internal
      integer iatom                                             ! internal
      integer idummy                                            ! internal
      integer istruct                                           ! internal
      integer refelem                                           ! internal (element index of working atom)
      integer num_atoms_local                                   ! internal
      integer num_atoms_ref                                     ! internal
      integer i1,i2,i3,i4,i5
      integer zelem_local(max_num_atoms)                        ! internal
      integer maxnum_funcvalues_local                           ! in
      integer num_funcvalues_local(nelem)                       ! in
      integer ncount,ncount1,ncount2                                   ! internal
      integer nstruct1,nstruct2                                 ! internal
      integer npoints,npoints1,npoints2                                 ! internal
      integer nstart1,nstart2                                   ! internal
      integer nend1,nend2                                       ! internal
      integer nref                                              ! internal
      integer localunit                                         ! unit for symmetry function file of training or test symfunctions
      integer temp_atom(nelem) 
      integer temp_struct(nelem) 
!!
      real*8 zdummy
      real*8 minvalue_local(nelem,maxnum_funcvalues_local)                ! in
      real*8 maxvalue_local(nelem,maxnum_funcvalues_local)                ! in
      real*8 avvalue_local(nelem,maxnum_funcvalues_local)                 ! in
      real*8 scmin_local                                                  ! in
      real*8 scmax_local                                                  ! in
      real*8 symfunction_local(maxnum_funcvalues_local,max_num_atoms)     ! internal
      real*8 symfunction1(maxnum_funcvalues_local,max_num_atoms,nblock)   ! internal
      real*8 symfunction2(maxnum_funcvalues_local,max_num_atoms,nblock)   ! internal
      real*8 deltag
      real*8 distance 
      real*8 clustercenter(nelem,1000,maxnum_funcvalues_local) 
      integer clustercentercounter(nelem,1000) 
!!
!! arrays for first block of points
      integer     num_atoms_list1(nblock)
      integer     zelem_list1(nblock,max_num_atoms)
      integer     index1(nblock)
!!
!! arrays for second block of points
      integer     num_atoms_list2(nblock)
      integer     zelem_list2(nblock,max_num_atoms)
      integer     index2(nblock)

      logical     ldone(nelem)             ! by element
      logical     lfound(nelem)            ! an atom of an element has been assigned to a cluster in this pass 
      logical     lnewcluster(nelem)       ! start new cluster for a given element

      character*20 filename


      write(ounit,*)'============================================================='
      if(iswitch.eq.0)then
        write(ounit,*)'Performing data clustering of the training set'
        localunit=symunit
        filename='function.data'
      elseif(iswitch.eq.1)then
        write(ounit,*)'Performing data clustering of the test set'
        localunit=tymunit
        filename='testing.data'
      else
        write(ounit,*)'ERROR: unknown iswitch in subroutine dataclustering'
        stop
      endif
      write(ounit,*)'Threshold for distance in symmetry function space: ',&
        dataclusteringthreshold1
      write(ounit,*)'============================================================='
!'
!!===========================================================
!! initializations
!!===========================================================
      num_clusters(:)=0 ! number of identified different clusters for each element
      clusterindex(:,:)=0 ! index number of the cluster the atom is part of (caution: cluster 1 for element A is a different cluster than cluster 1 for element B)
      cluster_atoms(:,:)=0 ! number of atoms included in each cluster for each element
      nstart1        =1
      nstruct2       =0
      ldone(:)       =.false. ! true if all atoms have been assigned to clusters
      lnewcluster(:) =.true.
      clustercenter(:,:,:)=0.0d0 
      clustercentercounter(:,:)=0 

!!===========================================================
!! find the first atom of each element not assigned to a cluster to create the new nucleus 
!!===========================================================
 50   continue
!      write(ounit,*)'Searching for new cluster nucleus'
!      do i1=1,nelem
!        write(ounit,*)'ldone ',element(i1),ldone(i1)
!      enddo
      lfound(:)=.false.
      ncount=nstruct_local ! remaining structures to be read
      istruct=0
      open(localunit,file=filename,form='formatted',status='old')
      rewind(localunit)
 51   continue             ! next block of structures
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
      do i2=1,npoints
        read(localunit,*)num_atoms_local
        istruct=istruct+1
!        write(ounit,*)'Read structure ',istruct
        do i1=1,num_atoms_local
          read(localunit,*)zelem_local(i1)  ! nuclear charge
          backspace(localunit)
          read(localunit,*)idummy,(symfunction_local(i3,i1),&
            i3=1,num_funcvalues_local(elementindex(zelem_local(i1))))
        enddo
        read(localunit,*)zdummy             ! energy line 
!! check atoms
        do i1=1,num_atoms_local
          if(.not.lfound(elementindex(zelem_local(i1))))then ! no nucleus for this element yet
            if(clusterindex(istruct,i1).eq.0)then ! this atom does not yet belong to a cluster    
              num_clusters(elementindex(zelem_local(i1)))=num_clusters(elementindex(zelem_local(i1)))+1
              clusterindex(istruct,i1)=num_clusters(elementindex(zelem_local(i1)))
              lfound(elementindex(zelem_local(i1)))=.true.
              temp_atom(elementindex(zelem_local(i1)))=i1
              temp_struct(elementindex(zelem_local(i1)))=istruct
              if(num_clusters(zelem_local(i1)).gt.1000)then
                write(ounit,*)'ERROR in dataclustering: redimension cluster_atoms array'
                stop !'
              endif
              cluster_atoms(elementindex(zelem_local(i1)),&
                num_clusters(elementindex(zelem_local(i1))))=&
                cluster_atoms(elementindex(zelem_local(i1)),&
                num_clusters(elementindex(zelem_local(i1))))+1
              nassigned(elementindex(zelem_local(i1)))&
                =nassigned(elementindex(zelem_local(i1)))+1
! include this nucleus in the calculation of the cluster's center of mass
              clustercentercounter(elementindex(zelem_local(i1)),&
                num_clusters(elementindex(zelem_local(i1))))&
                =clustercentercounter(elementindex(zelem_local(i1)),&
                num_clusters(elementindex(zelem_local(i1))))+1
              do i3=1,num_funcvalues_local(elementindex(zelem_local(i1)))
                clustercenter(elementindex(zelem_local(i1)),&
                  num_clusters(elementindex(zelem_local(i1))),i3)&
                  =clustercenter(elementindex(zelem_local(i1)),&
                  num_clusters(elementindex(zelem_local(i1))),i3)&
                  +symfunction_local(i3,i1)
              enddo ! i3
            endif
          endif
        enddo ! i1 loop over points
      enddo ! i2 loop over structures
      if(ncount.gt.0) goto 51 ! go to next block of structures 
 52   continue                ! full data set done
      close(localunit)

!! print identified new nuclei
      do i1=1,nelem
        if(lfound(i1))then
!          write(ounit,'(a8,x,a2,x,a40,i8,a7,i5)')'Element ',element(i1),&
!            ': new cluster nucleus found: structure ',temp_struct(i1),&
!            ' atom: ',temp_atom(i1)
        else
!          write(ounit,'(a8,x,a2,x,a)')'Element ',&
!            element(i1),': all atoms have been assigned to clusters'
          ldone(i1)=.true.
        endif
      enddo

!!===========================================================
!! Now we have a new nucleus for each element or all atoms of an element have been assigned already
!!===========================================================
!! now find all atoms belonging to the new nucleus. This requires several passes until no new
!! atom is found in a pass because clusters (e.g. with a chain shape) need to be identified atom by atom 
!!===========================================================
      npass=0
!      write(ounit,*)&
!        'Searching for all atoms belonging to the following clusters:'
      do i1=1,nelem
        if(.not.ldone(i1))then
          write(ounit,'(a,x,a2,x,a9,i8)')&
            'Now searching: Element ',element(i1),' cluster ',num_clusters(i1)
        endif
      enddo

!! jump here for the next pass through all structures for the same cluster
 60   continue
      nstart1        =1
      nstruct2       =0
      npass=npass+1
!      write(ounit,*)'Starting pass ',npass
      lfound(:)=.false. 
      nassigned(:)=0

!!===========================================================
!! Now we need to compare the symmetry function vector of a given atom (reference structure)
!! with the symmetry function vectors of all other atoms (trial structures) (avoid double counting)
!! for the current cluster (element-specific)
!!===========================================================
!!
!!===========================================================
!! read a block of reference structures
!!===========================================================
      ncount1=nstruct_local ! remaining structures to be read
 11   continue
      nstruct1=0 ! counter for structures in reference set
      nstart2=1  ! if we read new reference structures, start with all trial structures from the beginning
      if(ncount1.gt.nblock)then
        npoints1=nblock
        ncount1=ncount1-nblock
        nend1 =nstart1+npoints1-1
      else
        npoints1=ncount1
        ncount1=ncount1-npoints1
        nend1 =nstart1+npoints1-1
      endif
!      write(ounit,'(a8,i8,a,i8,a4,i8)')&
!        'reading ',npoints1,' reference structures from ',nstart1,' to ',nend1

!! Caution: the file needs to be reopened for each block of structures, a simple continuation
!! of reading is not possible, because the file will also be accessed for reading the trial 
!! structures below
      nref=0  ! absolute number of reference structure in full data set
      open(localunit,file=filename,form='formatted',status='old')
      rewind(localunit)
 20   continue
      nstruct1=nstruct1+1
      if(nstruct1.gt.nstruct_local)goto 21 ! all structures have been read
!! read previous reference structures as dummies to reach the correct starting position 
      if(nstruct1.lt.nstart1)then ! do not use these structures 
        read(localunit,*)num_atoms_local
        do i1=1,num_atoms_local
          read(localunit,*)idummy
        enddo
        read(localunit,*)zdummy     
!! if we have completed reading relevant structures of this block of structures
      elseif(nstruct1.gt.nend1)then  
        nstruct1=nstruct1-1
        goto 21
!! read relevant set of reference structures
      else 
        nref=nref+1
        index1(nref)=nstruct1 ! remember which structure of the data set this is
        read(localunit,*)num_atoms_list1(nref)
        do i1=1,num_atoms_list1(nref)
          read(localunit,*)idummy
          backspace(localunit)
          read(localunit,*)zelem_list1(nref,i1),(symfunction1(i2,i1,nref),&
            i2=1,num_funcvalues_local(elementindex(idummy)))
        enddo
        read(localunit,*)zdummy             ! energy line in function.data
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
      endif
      goto 20
 21   continue
      close(localunit)
      nstart1=nstart1+npoints1 ! prepare starting point for next block
!!================================================================
!! reading reference block done
!!================================================================

!!================================================================
!! read block of trial structures 
!!================================================================
      ncount2=nstruct_local
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
!      write(ounit,'(a8,i8,a,i8,a4,i8)')&
!        'reading ',npoints2,' trial structures from ',nstart2,' to ',nend2

      nref=0
      open(localunit,file=filename,form='formatted',status='old')
      rewind(localunit)

 40   continue
      nstruct2=nstruct2+1
!! read previous dummy structures which are not used
      if(nstruct2.gt.nstruct_local)goto 41
      if(nstruct2.lt.nstart2)then
        read(localunit,*)num_atoms_local
        do i1=1,num_atoms_local
          read(localunit,*)idummy
        enddo
        read(localunit,*)zdummy
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
        read(localunit,*)num_atoms_list2(nref)
        do i1=1,num_atoms_list2(nref)
          read(localunit,*)idummy
          backspace(localunit)
          read(localunit,*)zelem_list2(nref,i1),(symfunction2(i2,i1,nref),&
            i2=1,num_funcvalues_local(elementindex(idummy)))
        enddo
        read(localunit,*)zdummy             ! energy line 
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

      endif
      goto 40
 41   continue
      close(localunit)
      nstart2=nstart2+npoints2
!!================================================================
!! reading trial block done
!!================================================================

!!================================================================
!! compare now
!!================================================================
!! loop over all structures and atoms of reference block 
!!================================================================
!      write(ounit,*)'comparison is in progress '
      do i1=1,npoints1 ! reference structures
        do i2=1,npoints2 ! trial structures
          if(index1(i1).ne.index2(i2))then
!            write(ounit,'(a,2i6)')'Comparing structures gt ',index1(i1),index2(i2)
!! compare all atoms of all structures
            do i3=1,num_atoms_list1(i1) ! atoms of reference structure
!! compare only if reference structure atom belongs to the current cluster
              if(clusterindex(index1(i1),i3).eq.&
                num_clusters(elementindex(zelem_list1(i1,i3))))then
!                write(ounit,'(a,2i5,a,i5,x,a2)')'Reference atom ',index1(i1),&
!                  i3,' belongs to cluster ',&
!                  num_clusters(elementindex(zelem_list1(i1,i3))),&
!                  element(elementindex(zelem_list1(i1,i3)))
                do i4=1,num_atoms_list2(i2) ! atoms of trial structure
!! compare environments only if atoms have the same element and if atom has not been assigned yet
                  if(zelem_list1(i1,i3).eq.zelem_list2(i2,i4))then
!                    write(ounit,'(a,i6,x,i6,x,a2,x,i6,x,i6,x,a2)')&
!                      'Comparing atoms ',&
!                      index1(i1),&
!                      i3,element(elementindex(zelem_list1(i1,i3))),&
!                      index2(i2),&
!                      i4,element(elementindex(zelem_list2(i2,i4)))
                    if(clusterindex(index2(i2),i4).eq.0)then
                      deltag = 0.0d0
                      do i5=1,num_funcvalues_local(elementindex(zelem_list1(i1,i3)))
                        deltag = deltag + (symfunction1(i5,i3,i1)-symfunction2(i5,i4,i2))**2 
                      enddo ! i5
                      deltag = dsqrt(deltag)
!                      write(ounit,'(a,f14.6)')'distance is ',deltag
!! check if this trial belongs to the current cluster 
                      if(deltag.le.dataclusteringthreshold1)then
!                        write(ounit,*)npass,' assigning atom to cluster ',&
!                          num_clusters(elementindex(zelem_list1(i1,i3))),& 
!                          element(elementindex(zelem_list1(i1,i3))) !'
                        clusterindex(index2(i2),i4)&
                          =num_clusters(elementindex(zelem_list2(i2,i4)))
                        lfound(elementindex(zelem_list1(i1,i3)))=.true.
                        nassigned(elementindex(zelem_list1(i1,i3)))&
                          =nassigned(elementindex(zelem_list1(i1,i3)))+1
                        cluster_atoms(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list2(i2,i4))))&
                          =cluster_atoms(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list2(i2,i4))))+1
                        clustercentercounter(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list1(i1,i3))))&
                          =clustercentercounter(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list1(i1,i3))))+1
                        do i5=1,num_funcvalues_local(elementindex(zelem_list2(i2,i4)))
                          clustercenter(elementindex(zelem_list2(i2,i4)),&
                            num_clusters(elementindex(zelem_list2(i2,i4))),i5)&
                            =clustercenter(elementindex(zelem_list2(i2,i4)),&
                            num_clusters(elementindex(zelem_list2(i2,i4))),i5)&
                            +symfunction2(i5,i4,i2)
                        enddo ! i5
                      else
!                        write(ounit,*)'Atom already assigned ',clusterindex(index2(i2),i4)
                      endif
                    else
!                      write(ounit,'(a9,f14.6,a34,i6,a,a2)')&
!                        'distance ',deltag,&
!                        ' atom does not belong to cluster ',&
!                        num_clusters(elementindex(zelem_list1(i1,i3))),&
!                        ' for element ',element(elementindex(zelem_list1(i1,i3)))
                    endif
                  else
!                    write(ounit,*)'wrong element'
                  endif ! zelem_list
                enddo ! i4
              endif ! clusterindex
            enddo ! i3
          elseif(index1(i1).eq.index2(i2))then
!            write(ounit,'(a,2i6)')'Comparing structures eq ',index1(i1),index2(i2)
!! compare all atoms within the reference structure (here identical to the trial structure) 
            do i3=1,num_atoms_list1(i1)-1
!! compare only if reference structure atom belongs to the current cluster
              if(clusterindex(index1(i1),i3).eq.&
                num_clusters(elementindex(zelem_list1(i1,i3))))then
!                write(ounit,'(a,2i5,a,i5,x,a2)')'Reference atom ',index1(i1),&
!                  i3,' belongs to cluster ',&
!                  num_clusters(elementindex(zelem_list1(i1,i3))),&
!                  element(elementindex(zelem_list1(i1,i3)))
                do i4=i3+1,num_atoms_list2(i2)
!! compare environments only if atoms have the same element and if atom has not been assigned yet
                  if(zelem_list1(i1,i3).eq.zelem_list2(i2,i4))then
!                    write(ounit,'(a,i6,x,i6,x,a2,x,i6,x,i6,x,a2)')&
!                      'Comparing atoms ',&
!                      index1(i1),&
!                      i3,element(elementindex(zelem_list1(i1,i3))),&
!                      index2(i2),&
!                      i4,element(elementindex(zelem_list2(i2,i4)))
                    if(clusterindex(index2(i2),i4).eq.0)then
                      deltag = 0.0d0
                      do i5=1,num_funcvalues_local(elementindex(zelem_list1(i1,i3)))
                        deltag = deltag + (symfunction1(i5,i3,i1)-symfunction2(i5,i4,i2))**2 
                      enddo ! i5
                      deltag = dsqrt(deltag)
!                      write(ounit,'(a,f14.6)')'distance is ',deltag
!! check if this trial belongs to the current cluster 
                      if(deltag.le.dataclusteringthreshold1)then
!                        write(ounit,*)npass,&
!                          ' assigning atom to cluster ',&
!                          num_clusters(elementindex(zelem_list1(i1,i3))),&
!                          element(elementindex(zelem_list1(i1,i3)))
                        clusterindex(index2(i2),i4)&
                          =num_clusters(elementindex(zelem_list2(i2,i4)))
                        lfound(elementindex(zelem_list1(i1,i3)))=.true.
                        nassigned(elementindex(zelem_list1(i1,i3)))&
                          =nassigned(elementindex(zelem_list1(i1,i3)))+1
                        cluster_atoms(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list2(i2,i4))))&
                          =cluster_atoms(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list2(i2,i4))))+1
                        clustercentercounter(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list1(i1,i3))))&
                          =clustercentercounter(elementindex(zelem_list1(i1,i3)),&
                          num_clusters(elementindex(zelem_list1(i1,i3))))+1
                        do i5=1,num_funcvalues_local(elementindex(zelem_list2(i2,i4)))
                          clustercenter(elementindex(zelem_list2(i2,i4)),&
                            num_clusters(elementindex(zelem_list2(i2,i4))),i5)&
                            =clustercenter(elementindex(zelem_list2(i2,i4)),&
                            num_clusters(elementindex(zelem_list2(i2,i4))),i5)&
                            +symfunction2(i5,i4,i2)
                        enddo ! i5
                      else
!                        write(ounit,'(a9,f14.6,a34,i6,a,a2)')&
!                          'distance ',deltag,&
!                          ' atom does not belong to cluster ',&
!                          num_clusters(elementindex(zelem_list1(i1,i3))),&
!                          ' for element ',element(elementindex(zelem_list1(i1,i3)))
                      endif !'
                    else
!                      write(ounit,*)'Atom already assigned ',clusterindex(index2(i2),i4)
                    endif
                  else
!                    write(ounit,*)'wrong element'
                  endif ! zelem_list1
                enddo ! i4
              endif ! clusterindex
            enddo ! i3
          endif
        enddo ! i2
      enddo ! i1
!!================================================================
!! end loop over all structures and atoms of reference block 
!!================================================================

!!================================================================
!! read next block of trial structures 
!!================================================================
      if(ncount2.gt.0) goto 31 ! do next block of comparison points
!!
!!================================================================
!! read next block of reference points
!!================================================================
      if(ncount1.gt.0) goto 11 ! do next block of reference points

!! Print how many atoms have been assigned in this pass
!      do i1=1,nelem
!        write(ounit,'(a,x,a2,x,l,x,i6)')&
!          ' Assignment ',element(i1),lfound(i1),&
!          nassigned(i1)
!      enddo

!! If a successful assignment has been made in this pass, do it again 
      do i1=1,nelem
        if(lfound(i1))then
          goto 60
        endif
      enddo

!! If no new atom has been assigned in this pass but there are
!! still atoms to be assigned generate a new nucleus
      do i1=1,nelem
        if(.not.ldone(i1))then
          goto 50
        endif
      enddo

!!
!!================================================================
!! print final results 
!!================================================================
      do i1=1,nelem
        write(ounit,*)'-------------------------------------------------------------'
        if(iswitch.eq.0)then
          write(ounit,'(a,a2,a)')'Element ',element(i1),': (Training set)'
        elseif(iswitch.eq.1)then
          write(ounit,'(a,a2,a)')'Element ',element(i1),': (Test set)'
        endif
        write(ounit,*)'            Cluster  atoms'
        do i2=1,num_clusters(i1)
          write(ounit,'(a8,x,a2,x,i8,x,i6)')&
            'Cluster ',element(i1),i2,cluster_atoms(i1,i2)
        enddo ! i2
      enddo ! i1
      write(ounit,*)'-------------------------------------------------------------'

!!================================================================
!! calculate final cluster centers 
!!================================================================
      if(iswitch.eq.0)then
        write(ounit,*)'Calculating the centers of the clusters (Training set):'
      elseif(iswitch.eq.1)then
        write(ounit,*)'Calculating the centers of the clusters (Test set):'
      endif
      do i1=1,nelem
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,a2,a1)')'Element ',element(i1),':'
        do i2=1,num_clusters(i1)
          clustercenter(i1,i2,:)=clustercenter(i1,i2,:)&
            /dble(clustercentercounter(i1,i2))
          if(iswitch.eq.0)then
            write(ounit,'(a24,x,a2,x,i8,x,100f14.6)')&
              'Training Cluster center ',element(i1),i2,&
              (clustercenter(i1,i2,i3),i3=1,num_funcvalues_local(i1))
          elseif(iswitch.eq.1)then
            write(ounit,'(a24,x,a2,x,i8,x,100f14.6)')&
              'Test Cluster center ',element(i1),i2,&
              (clustercenter(i1,i2,i3),i3=1,num_funcvalues_local(i1))
          endif
        enddo ! i2
      enddo ! i1
      write(ounit,*)'-------------------------------------------------------------'
     
!!================================================================
!! calculate distances between cluster centers 
!!================================================================
      if(iswitch.eq.0)then
        write(ounit,*)'Calculating the distances between the centers of the clusters (Training set):'
      elseif(iswitch.eq.1)then
        write(ounit,*)'Calculating the distances between the centers of the clusters (Test set):'
      endif
      do i1=1,nelem
        if(num_clusters(i1).gt.0)then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,*)'Element ',element(i1)
          do i2=1,num_clusters(i1)-1
            distance=0.0d0
            do i3=i2+1,num_clusters(i1)
              do i4=1,num_funcvalues_local(i1)
                distance=distance+(clustercenter(i1,i2,i4)-clustercenter(i1,i3,i4))**2.0d0
              enddo ! i4
              distance=dsqrt(distance)
              if(iswitch.eq.0)then
                write(ounit,'(a17,x,2i6,x,a11,f14.6)')&
                  'Training Centers ',i2,i3,' distance: ',distance
              elseif(iswitch.eq.1)then
                write(ounit,'(a17,x,2i6,x,a11,f14.6)')&
                  'Test Centers ',i2,i3,' distance: ',distance
              endif
            enddo ! i3
          enddo ! i2
        endif
      enddo ! i1
      write(ounit,*)'-------------------------------------------------------------'

!!================================================================
!! write full list of all atoms and their assignments 
!!================================================================
      if(dataclusteringthreshold2.gt.1.0d0)then
        if(iswitch.eq.0)then
          write(ounit,*)'Cluster assignment of all atoms in the training set:'
        elseif(iswitch.eq.1)then
          write(ounit,*)'Cluster assignment of all atoms in the test set:'
        endif
        do i1=1,nstruct_local
          do i2=1,max_num_atoms
            if(clusterindex(i1,i2).gt.0)then ! only existing atoms
              if(iswitch.eq.0)then
                write(ounit,*)'Training Structure, Atom, Cluster: ',i1,i2,clusterindex(i1,i2)
              elseif(iswitch.eq.1)then
                write(ounit,*)'Test Structure, Atom, Cluster: ',i1,i2,clusterindex(i1,i2)
              endif
            endif
          enddo
        enddo
      endif


!!================================================================
!! Now we could characterize the clusters: centers, distance matrix between centers (allows 3D plot), standard deviation and diameters of clusters etc. 
!!================================================================

      return
      end
