!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!! - main.f90
!!
      subroutine getpairsymfunctions(iseed)
!!
      use fileunits
      use globaloptions
      use mode1options
      use symfunctions
      use nnpair
      use nnewald
      use structures
!!
      implicit none
!!
      integer npoints                                                      ! internal
      integer ncount                                                       ! internal
      integer iseed                                                        ! in/out
      integer numtrain                                                     ! internal 
      integer numtest                                                      ! internal
      integer numrej                                                       ! internal 
      integer pointnumber                                                  ! internal
      integer num_atoms_element_list(nblock,nelem)                         ! internal 
      integer i1,i2                                                        ! internal
      integer ndone                                                        ! internal

      real*8 chargemin(nelem)                                              ! internal 
      real*8 chargemax(nelem)                                              ! internal 
      real*8 minvaluee(nelem,maxnum_funcvaluese)                           ! internal 
      real*8 maxvaluee(nelem,maxnum_funcvaluese)                           ! internal 
      real*8 avvaluee(nelem,maxnum_funcvaluese)                            ! internal 
      real*8 dummy                                                         ! internal

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      write(ounit,*)'Calculating Symmetry Functions'
      write(ounit,*)'for ',totnum_structures,' structures'
      write(ounit,*)'-------------------------------------------------------------'
!!'
!!    allocate arrays of structures module
      allocate (num_atoms_list(nblock))
      allocate (num_pairs_list(nblock))
      allocate (zelem_list(nblock,max_num_atoms))
      allocate (zelemp_list(2,nblock,max_num_pairs))
      allocate (lattice_list(3,3,nblock))
      allocate (xyzstruct_list(3,max_num_atoms,nblock))
      allocate (totalcharge_list(nblock))
      allocate (totalenergy_list(nblock))
      allocate (shortenergy_list(nblock))
      allocate (ewaldenergy_list(nblock))
      allocate (xyzforce_list(3,max_num_atoms,nblock))
      allocate (ewaldforce_list(3,max_num_atoms,nblock))
      allocate (shortforce_list(3,max_num_atoms,nblock))
      allocate (totforce_list(3,max_num_atoms,nblock))
      allocate (atomcharge_list(nblock,max_num_atoms))
      allocate (atomenergy_list(nblock,max_num_atoms))
      allocate (lperiodic_list(nblock))
      allocate (elementsymbol_list(nblock,max_num_atoms))
!!
!! initializations
      numtrain           = 0
      numtest            = 0
      numrej             = 0
      pointnumber        = 0
      maxcutoffp         = 0.0d0
      maxcutoffe         = 0.0d0
      weights_ewald(:,:) = 0.0d0
      minvaluee(:,:)     = 0.0d0
      maxvaluee(:,:)     = 0.0d0
      avvaluee(:,:)      = 0.0d0
!!
!! check if the number of function values  
!! the same as input nodes in input.nn
      do i1=1,npairs
        if(num_funcvaluesp(i1).ne.nodes_pair(0,i1))then
          write(ounit,*)'Error: num_funcvaluesp .ne. nodes_short(0)'
          write(ounit,*)i1,num_funcvaluesp(i1),nodes_pair(0,i1)
          stop
        endif
      enddo ! i1
      do i1=1,nelem
        if(num_funcvaluese(i1).ne.nodes_ewald(0,i1))then
          write(ounit,*)'Error: num_funcvaluese .ne. nodes_ewald(0)'
          write(ounit,*)i1,num_funcvaluese(i1),nodes_ewald(0,i1)
          stop
        endif
      enddo ! i1
!!
!!    determine maxcutoffp
      do i2=1,npairs
        do i1=1,num_funcvaluesp(i2)
          maxcutoffp=max(maxcutoffp,funccutoffp(i1,i2))
        enddo ! i1
      enddo

!!    determine maxcutoffe
      do i2=1,nelem
        do i1=1,num_funcvaluese(i2)
          maxcutoffe=max(maxcutoffe,funccutoffe(i1,i2))
        enddo ! i1
      enddo ! i2
!!
!! get NN data for charges if electrostatic forces need to be calculated
      if(lelec.and.(etype.eq.1).and.luseforces.and.(.not.lfixedcharges))then
        write(ounit,*)'Reading charge fit data for determination of electrostatic forces'
        write(ounit,*)'-------------------------------------------------------------'
!!'
!! read scalinge data for electrostatics
        call readscale(nelem,3,&
          maxnum_funcvaluese,num_funcvaluese,&
          minvaluee,maxvaluee,avvaluee,&
          dummy,dummy,chargemin,chargemax)
!!
!! read electrostatic weights
        call readweights(1,nelem,&
          maxnum_weightsewald,num_weightsewald,&
          weights_ewald)
!!
      endif ! lelec
!!
!! open files
      open(dataunit,file='input.data',form='formatted')
      rewind(dataunit)
      open(symunit,file='function.data',form='formatted',status='replace')
      rewind(symunit)
      open(tymunit,file='testing.data',form='formatted',status='replace')
      rewind(tymunit)
      open(trainstructunit,file='trainstruct.data',form='formatted',status='replace')
      rewind(trainstructunit)
      open(teststructunit,file='teststruct.data',form='formatted',status='replace')
      rewind(teststructunit)
      open(symeunit,file='functione.data',form='formatted',status='replace')
      rewind(symeunit)
      open(tymeunit,file='testinge.data',form='formatted',status='replace')
      rewind(tymeunit)
      open(trainfunit,file='trainforces.data',form='formatted',status='replace')
      rewind(trainfunit)
      open(trainfeunit,file='trainforcese.data',form='formatted',status='replace')
      rewind(trainfeunit)
      open(testfunit,file='testforces.data',form='formatted',status='replace')
      rewind(testfunit)
      open(testfeunit,file='testforcese.data',form='formatted',status='replace')
      rewind(testfeunit)
!!
!! calculate the symmetry functions for all structures
      ncount=totnum_structures
      ndone=0
!! read a set of points
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! read data file(s)
      call readstructures(npoints,num_atoms_element_list)
!!
!! remove atomic energies from total energies if requested
      if(lremoveatomenergies)then
        call removeatoms(nblock,npoints, &
          num_atoms_list,zelem_list,&
          num_atoms_element_list,&
          totalenergy_list,atomenergy_list)
      endif
!!
!! if charges are fixed, overwrite reference charges
!! this is needed to subtract the correct electrostatic energies from the total energy
      if(lfixedcharges)then
        do i1=1,npoints
          do i2=1,num_atoms_list(i1)
            atomcharge_list(i1,i2)=fixedcharge(elementindex(zelem_list(i1,i2)))
          enddo
        enddo
      endif
!!
!! calculate and write symmetry functions for a block of structures
      call calcpairfunctions(&
        npoints,ndone,&
        iseed,numtrain,numtest,numrej,pointnumber,&
        minvaluee,maxvaluee,avvaluee)
!!
      ndone=ndone+npoints
      if(ncount.gt.0) goto 10
!!
!! close files
      close(dataunit)
      close(symunit)
      close(tymunit)
      close(trainstructunit)
      close(teststructunit)
      close(symeunit)
      close(tymeunit)
      close(trainfunit)
      close(trainfeunit)
      close(testfunit)
      close(testfeunit)
!!
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Number of fitting points: ',numtrain
      write(ounit,*)'Number of testing points: ',numtest
      write(ounit,*)'Number of rejected points:',numrej
!! '
!!    deallocate arrays of structures module
      deallocate (num_atoms_list)
      deallocate (num_pairs_list)
      deallocate (zelem_list)
      deallocate (zelemp_list)
      deallocate (lattice_list)
      deallocate (xyzstruct_list)
      deallocate (totalcharge_list)
      deallocate (totalenergy_list)
      deallocate (shortenergy_list)
      deallocate (ewaldenergy_list)
      deallocate (xyzforce_list)
      deallocate (shortforce_list)
      deallocate (ewaldforce_list)
      deallocate (totforce_list)
      deallocate (atomcharge_list)
      deallocate (lperiodic_list)
      deallocate (elementsymbol_list)
!!
      return
      end
