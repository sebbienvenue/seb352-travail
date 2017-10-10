!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine erroranalysis(ntrain,ntest)
!!
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer ntrain                                  ! in
      integer ntest                                   ! in
      integer i1,i2,i3,i4                             ! internal
      integer idummy                                  ! internal
      integer icount                                  ! internal
      integer istruct                                 ! internal
      integer iatom                                   ! internal
      integer zelemtrain(ntrain,max_num_atoms)        ! internal
      integer zelemtest(ntest,max_num_atoms)          ! internal
      integer nstepse
      integer nstepsf
      integer nstepsq
      integer, dimension(:), allocatable :: counte 
      integer, dimension(:,:), allocatable :: countf 
      integer, dimension(:,:), allocatable :: countq 
      integer abovee
      integer abovef(nelem)
      integer aboveq(nelem)
      integer num_atoms
      integer num_pairs
!!   
      real*8 rdummy1                                  ! internal
      real*8 rdummy2                                  ! internal
      real*8 etrainref(ntrain)                        ! internal
      real*8 etrainnn(ntrain)                         ! internal
      real*8 etestref(ntest)                          ! internal
      real*8 etestnn(ntest)                           ! internal
      real*8 ftrainref(3,max_num_atoms,ntrain)        ! internal
      real*8 ftrainnn(3,max_num_atoms,ntrain)         ! internal
      real*8 ftestref(3,max_num_atoms,ntest)          ! internal
      real*8 ftestnn(3,max_num_atoms,ntest)           ! internal
      real*8 qtrainref(max_num_atoms,ntrain)          ! internal
      real*8 qtrainnn(max_num_atoms,ntrain)           ! internal
      real*8 qtestref(max_num_atoms,ntest)            ! internal
      real*8 qtestnn(max_num_atoms,ntest)             ! internal
      real*8 de                                       
      real*8 df                                       
      real*8 dq                                       
      real*8 erroremin                                       
      real*8 errorfmin                                       
      real*8 errorqmin                                       
      real*8, dimension(:), allocatable :: egrid 
      real*8, dimension(:), allocatable :: fgrid 
      real*8, dimension(:), allocatable :: qgrid 
      real*8 error                                    ! internal
      real*8 tounit                                   ! internal
      real*8 errormaxe
      real*8 errormaxf(nelem)
      real*8 errormaxq(nelem)
!! 
      logical lperiodic_local
!!
      character*25 filename                           ! internal
      character*50 dummy                              ! internal
!! 
!!
!! Definition of grids
      de       =analyze_error_energy_step
      df       =analyze_error_force_step
      dq       =analyze_error_charge_step
!!      de       =0.01d0
!!      df       =0.01d0
!!      dq       =0.001d0
      nstepse  =101
      nstepsf  =101
      nstepsq  =101
      erroremin=0.00d0                                       
      errorfmin=0.00d0                                       
      errorqmin=0.00d0                                       
!! allocate grids
      allocate(counte(nstepse)) 
      allocate(countf(nstepsf,nelem)) 
      allocate(countq(nstepsq,nelem)) 
      counte(:)=0
      countf(:,:)=0
      countq(:,:)=0
      allocate(egrid(nstepse))
      allocate(fgrid(nstepsf))
      allocate(qgrid(nstepsq))
      egrid(1)=erroremin
      do i1=2,nstepse
        egrid(i1)=egrid(i1-1)+de
      enddo
      fgrid(1)=errorfmin
      do i1=2,nstepsf
        fgrid(i1)=fgrid(i1-1)+df
      enddo
      qgrid(1)=errorqmin
      do i1=2,nstepsq
        qgrid(i1)=qgrid(i1-1)+dq
      enddo
!!
      write(ounit,*)'============================================================='
      write(ounit,*)'Analysis of fitting errors:'
      write(ounit,*)'Note: maxenergy and maxforce are not yet taken into account!'
      write(ounit,*)'-------------------------------------------------------------'
!! initializations
      etrainref(:)    =0.0d0 
      etrainnn(:)     =0.0d0
      etestref(:)     =0.0d0 
      etestnn(:)      =0.0d0
      ftrainref(:,:,:)=0.0d0
      ftrainnn(:,:,:) =0.0d0
      ftestref(:,:,:) =0.0d0
      ftestnn(:,:,:)  =0.0d0
      qtrainref(:,:)  =0.0d0
      qtrainnn(:,:)   =0.0d0
      qtestref(:,:)   =0.0d0
      qtestnn(:,:)    =0.0d0
      zelemtrain(:,:) =0 
      zelemtest(:,:)  =0 
      abovee          =0
      abovef(:)       =0
      aboveq(:)       =0
      errormaxe       =0.0d0
      errormaxf(:)    =0.0d0
      errormaxq(:)    =0.0d0
!! setting the energy and force unit converter
      if(fitting_unit.eq.1)then
        tounit               = 27.211d0   ! energy conversion Ha to eV
      elseif(fitting_unit.eq.2)then
        tounit               = 1.d0   !  stay with Ha
      endif
!!
      if(lshort.and.lwritetrainpoints)then
        write(ounit,*)'Short range energy (training set):'
        if(fitting_unit.eq.1)then
          write(ounit,'(2a)')'                                eV/atom:',&
           '               Points '
        elseif(fitting_unit.eq.2)then
          write(ounit,'(2a)')'                                Ha/atom:',&
           '               Points '
        endif
!! get the file name for the training data '
        filename='trainpoints.000000.out' 
        if(nepochs.gt.999)then
          write(filename(15:18),'(i4)')nepochs
        elseif(nepochs.gt.99)then
          write(filename(16:18),'(i3)')nepochs
        elseif(nepochs.gt.9)then
          write(filename(17:18),'(i2)')nepochs
        else
          write(filename(18:18),'(i1)')nepochs
        endif
!!
        open(trainpxunit,file=filename,form='formatted',status='old')
!!          write(ounit,*)'Opening ',filename
          rewind(trainpxunit)
          read(trainpxunit,*)dummy
          do i1=1,ntrain
            read(trainpxunit,*)idummy,etrainref(i1),etrainnn(i1) 
          enddo
        close(trainpxunit)
!! analyze data
        do i1=1,ntrain
          error=abs(etrainref(i1)-etrainnn(i1))*tounit
          errormaxe=max(errormaxe,error)
          do i2=2,nstepse
            if(error.lt.egrid(i2))then
              counte(i2)=counte(i2)+1
              goto 33
            endif
          enddo ! i2
          abovee=abovee+1
 33       continue 
        enddo ! i1
!! write results
        do i1=2,nstepse
          write(ounit,'(a,i5,f10.4,a11,f10.4,i10)')' NNerrorEtrain ',&
            i1-1,egrid(i1-1),' < error < ',egrid(i1),counte(i1)
        enddo
        write(ounit,'(a,f10.4,a,i8)')' Points with error above ',&
          egrid(nstepse),' in training set: ',abovee
        if(fitting_unit.eq.1)then
          write(ounit,'(a,f10.4,a)')' Maximum error Eshort in training set:   ',errormaxe,' eV/atom'
        elseif(fitting_unit.eq.2)then
          write(ounit,'(a,f10.4,a)')' Maximum error Eshort in training set:   ',errormaxe,' Ha/atom'
        endif
        write(ounit,*)'-------------------------------------------------------------'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! reinitialize counters:
        errormaxe=0.0d0
        abovee   =0
        counte(:)=0 
!! get the file name for the testing data
        write(ounit,*)'Short range energy (test set):'
        if(fitting_unit.eq.1)then
          write(ounit,'(2a)')'                               eV/atom:',&
           '               Points '
        elseif(fitting_unit.eq.2)then
          write(ounit,'(2a)')'                               Ha/atom:',&
           '               Points '
        endif
!! get the file name for the test data '
        filename='testpoints.000000.out'
        if(nepochs.gt.999)then
          write(filename(14:17),'(i4)')nepochs
        elseif(nepochs.gt.99)then
          write(filename(15:17),'(i3)')nepochs
        elseif(nepochs.gt.9)then
          write(filename(16:17),'(i2)')nepochs
        else
          write(filename(17:17),'(i1)')nepochs
        endif
!!
        open(testpxunit,file=filename,form='formatted',status='old')
          rewind(testpxunit)
          read(testpxunit,*)dummy
          do i1=1,ntest
            read(testpxunit,*)idummy,etestref(i1),etestnn(i1)
          enddo
        close(testpxunit)
!! analyze data
        do i1=1,ntest
          error=abs(etestref(i1)-etestnn(i1))*tounit
          errormaxe=max(errormaxe,error)
          do i2=2,nstepse
            if(error.lt.egrid(i2))then
              counte(i2)=counte(i2)+1
              goto 34
            endif
          enddo ! i2
          abovee=abovee+1
 34       continue
        enddo ! i1
!! write results
        do i1=2,nstepse
          write(ounit,'(a,i5,f10.4,a11,f10.4,i10)')' NNerrorEtest ',&
            i1-1,egrid(i1-1),' < error < ',egrid(i1),counte(i1)
        enddo
        write(ounit,'(a,f10.4,a,i10)')' Points with error above ',&
          egrid(nstepse),' in test set:  ',abovee
        if(fitting_unit.eq.1)then
          write(ounit,'(a,f10.4,a)')' Maximum error Eshort in test set:      ',errormaxe,' eV/atom'
        elseif(fitting_unit.eq.2)then
          write(ounit,'(a,f10.4,a)')' Maximum error Eshort in test set:      ',errormaxe,' Ha/atom'
        endif
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FORCES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(lshort.and.lwritetrainforces)then
        write(ounit,*)'Short range forces:'
        write(ounit,*)'-------------------------------------------------------------'
!! get the file name for the training data '
        filename='trainforces.000000.out'
        if(nepochs.gt.999)then
          write(filename(15:18),'(i4)')nepochs
        elseif(nepochs.gt.99)then
          write(filename(16:18),'(i3)')nepochs
        elseif(nepochs.gt.9)then
          write(filename(17:18),'(i2)')nepochs
        else
          write(filename(18:18),'(i1)')nepochs
        endif
!!
        open(trainfxunit,file=filename,form='formatted',status='old')
          rewind(trainfxunit)
          read(trainfxunit,*)dummy
          icount=1
 10       continue
          read(trainfxunit,*,END=11)istruct,iatom
          backspace(trainfxunit)
          read(trainfxunit,*)idummy,idummy,dummy,&
            ftrainref(icount,iatom,istruct),ftrainnn(icount,iatom,istruct)
          icount=icount+1
          if(icount.eq.4)icount=1
          goto 10
 11       continue
        close(trainfxunit)
!! we also need the nuclear charges to distinguish between elements
        open(symunit,file='function.data',form='formatted',status='old')
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        icount=1 ! counter for structures
        rewind(symunit)
        rewind(trainstructunit)
        zelemtrain(:,:)=0
 40     continue
!! get num_atoms
        if(nn_type_short.eq.1)then
          read(symunit,*,END=41)num_atoms 
          do i1=1,num_atoms
            read(symunit,*)idummy
          enddo
        elseif(nn_type_short.eq.2)then
          read(symunit,*,END=41)num_atoms,num_pairs 
          do i1=1,num_pairs
            read(symunit,*)idummy
          enddo
        endif
        read(symunit,*)rdummy1
!! now get nuclear charges of atoms
        read(trainstructunit,*)idummy,lperiodic_local
        if(lperiodic_local)then
          do i1=1,3
            read(trainstructunit,*)rdummy1
          enddo 
        endif
        do i1=1,num_atoms
          read(trainstructunit,*)zelemtrain(icount,i1)
        enddo !i1
        icount=icount+1
        goto 40
 41     continue
        close(symunit)
        close(trainstructunit)
!! analyze data
        do i1=1,ntrain
          do i2=1,max_num_atoms
            do i3=1,3
              if(zelemtrain(i1,i2).gt.0)then ! check if atom exists
                error=abs(ftrainref(i3,i2,i1)-ftrainnn(i3,i2,i1))
                errormaxf(elementindex(zelemtrain(i1,i2)))=&
                  max(errormaxf(elementindex(zelemtrain(i1,i2))),error)
                do i4=2,nstepsf
                  if(error.lt.fgrid(i4))then
                    countf(i4,elementindex(zelemtrain(i1,i2)))=&
                      countf(i4,elementindex(zelemtrain(i1,i2)))+1
                    goto 53
                  endif
                enddo ! i4
                abovef(elementindex(zelemtrain(i1,i2)))=&
                  abovef(elementindex(zelemtrain(i1,i2)))+1
 53             continue
              endif
            enddo ! i3
          enddo ! i2
        enddo ! i1
!! write results
        do i1=1,nelem
          write(ounit,*)'Training set force error analysis for element ',element(i1),':'
          if(fitting_unit.eq.1)then
            write(ounit,*)'                                 eV/Bohr                Points'
          elseif(fitting_unit.eq.2)then
            write(ounit,*)'                                 Ha/Bohr                Points'
          endif
          do i2=2,nstepsf
            write(ounit,'(a,a2,i5,f10.4,a11,f10.4,i10)')' NNerrorFtrain ',element(i1),&
              i2-1,fgrid(i2-1),' < error < ',fgrid(i2),countf(i2,i1)
          enddo
          write(ounit,'(a,f9.4,a,i6)')' Points with force error above ',&
            fgrid(nstepsf),' in training set:',abovef(i1)
          if(fitting_unit.eq.1)then
            write(ounit,'(a,f10.4,a)')' Maximum force error in training set:      ',errormaxf(i1),' eV/Bohr'
          elseif(fitting_unit.eq.2)then
            write(ounit,'(a,f10.4,a)')' Maximum force error in training set:      ',errormaxf(i1),' Ha/Bohr'
          endif
          write(ounit,*)'-------------------------------------------------------------'
        enddo ! i1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! test set force error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! reinitialize counters:
        errormaxf(:)=0.0d0
        abovef(:)   =0
        countf(:,:) =0 
!! get the file name for the testing data
        filename='testforces.000000.out'
        if(nepochs.gt.999)then
          write(filename(14:17),'(i4)')nepochs
        elseif(nepochs.gt.99)then
          write(filename(15:17),'(i3)')nepochs
        elseif(nepochs.gt.9)then
          write(filename(16:17),'(i2)')nepochs
        else
          write(filename(17:17),'(i1)')nepochs
        endif
!!
        open(testfunit,file=filename,form='formatted',status='old')
          rewind(testfunit)
          read(testfunit,*)dummy
          icount=1
 60       continue
          read(testfunit,*,END=61)istruct,iatom
          backspace(testfunit)
          read(testfunit,*)idummy,idummy,dummy,&
            ftestref(icount,iatom,istruct),ftestnn(icount,iatom,istruct)
          icount=icount+1
          if(icount.eq.4)icount=1
          goto 60
 61       continue
        close(testfunit)
!! we also need the nuclear charges to distinguish between elements
        open(tymunit,file='testing.data',form='formatted',status='old')
        open(teststructunit,file='teststruct.data',form='formatted',status='old')
        icount=1 ! counter for structures
        rewind(tymunit)
        rewind(teststructunit)
        zelemtest(:,:)=0
 70     continue
!! read num_atoms
        if(nn_type_short.eq.1)then
          read(tymunit,*,END=71)num_atoms
          do i1=1,num_atoms
            read(tymunit,*)idummy
          enddo
        elseif(nn_type_short.eq.2)then
          read(tymunit,*,END=71)num_atoms,num_pairs
          do i1=1,num_pairs
            read(tymunit,*)idummy
          enddo
        endif
        read(tymunit,*)rdummy1
!! now read nuclear charges
        read(teststructunit,*)idummy,lperiodic_local
        if(lperiodic_local)then
          do i1=1,3
            read(teststructunit,*)rdummy1
          enddo 
        endif
        do i1=1,num_atoms
         read(teststructunit,*)zelemtest(icount,i1)
        enddo
        icount=icount+1
        goto 70
 71     continue
        close(tymunit)
        close(teststructunit)
!! analyze data
        do i1=1,ntest
          do i2=1,max_num_atoms
            do i3=1,3
              if(zelemtest(i1,i2).gt.0)then ! check if atom exists
                error=abs(ftestref(i3,i2,i1)-ftestnn(i3,i2,i1))
                errormaxf(elementindex(zelemtest(i1,i2)))=&
                  max(errormaxf(elementindex(zelemtest(i1,i2))),error)
                do i4=2,nstepsf
                  if(error.lt.fgrid(i4))then
                    countf(i4,elementindex(zelemtest(i1,i2)))=&
                      countf(i4,elementindex(zelemtest(i1,i2)))+1
                    goto 73
                  endif
                enddo ! i4
                abovef(elementindex(zelemtest(i1,i2)))=&
                  abovef(elementindex(zelemtest(i1,i2)))+1
 73             continue
              endif
            enddo ! i3
          enddo ! i2
        enddo ! i1
!! write results
        do i1=1,nelem
          write(ounit,*)'Test set force error analysis for element ',element(i1),':'
          if(fitting_unit.eq.1)then
            write(ounit,*)'                                eV/Bohr                Points'
          elseif(fitting_unit.eq.2)then
            write(ounit,*)'                                Ha/Bohr                Points'
          endif
          do i2=2,nstepsf
            write(ounit,'(a,a2,i5,f10.4,a11,f10.4,i10)')' NNerrorFtest ',element(i1),&
              i2-1,fgrid(i2-1),' < error < ',fgrid(i2),countf(i2,i1)
          enddo
          write(ounit,'(a,f9.4,a,i6)')' Points with force error above ',&
            fgrid(nstepsf),' in test set:   ',abovef(i1)
          if(fitting_unit.eq.1)then
            write(ounit,'(a,f10.4,a)')' Maximum force error in test set:         ',errormaxf(i1),' eV/Bohr'
          elseif(fitting_unit.eq.1)then
            write(ounit,'(a,f10.4,a)')' Maximum force error in test set:         ',errormaxf(i1),' Ha/Bohr'
          endif
          write(ounit,*)'-------------------------------------------------------------'
        enddo ! i1
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CHARGES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(lelec.and.lwritetraincharges)then
        write(ounit,*)'Atomic charges:'
        write(ounit,*)'-------------------------------------------------------------'
!! get the file name for the training data '
        filename='traincharges.000000.out'
        if(nepochs.gt.999)then
          write(filename(16:19),'(i4)')nepochs
        elseif(nepochs.gt.99)then
          write(filename(17:19),'(i3)')nepochs
        elseif(nepochs.gt.9)then
          write(filename(18:19),'(i2)')nepochs
        else
          write(filename(19:19),'(i1)')nepochs
        endif
!!
        open(trainqxunit,file=filename,form='formatted',status='old')
          rewind(trainqxunit)
          zelemtrain(:,:)=0
          read(trainqxunit,*)dummy
 20       continue
          read(trainqxunit,*,END=21)istruct,iatom
          backspace(trainqxunit)
          read(trainqxunit,*)idummy,idummy,&
            zelemtrain(istruct,iatom),&
            qtrainref(iatom,istruct),qtrainnn(iatom,istruct)
          goto 20
 21       continue
        close(trainqxunit)
!! analyze data
        do i1=1,ntrain
          do i2=1,max_num_atoms
            if(zelemtrain(i1,i2).gt.0)then ! check if atom exists
              error=abs(qtrainref(i2,i1)-qtrainnn(i2,i1))
              errormaxq(elementindex(zelemtrain(i1,i2)))=&
                max(errormaxq(elementindex(zelemtrain(i1,i2))),error)
              do i3=2,nstepsq
                if(error.lt.qgrid(i3))then
                  countq(i3,elementindex(zelemtrain(i1,i2)))=&
                    countq(i3,elementindex(zelemtrain(i1,i2)))+1
                  goto 23
                endif
              enddo ! i3
              aboveq(elementindex(zelemtrain(i1,i2)))=&
                aboveq(elementindex(zelemtrain(i1,i2)))+1
 23           continue
            endif
          enddo ! i2
        enddo ! i1
!! write results
        do i1=1,nelem
          write(ounit,*)'Training set charge error analysis for element ',element(i1),':'
          write(ounit,*)'                                    e                   Points'
          do i2=2,nstepsq
            write(ounit,'(a,a2,i5,f10.4,a11,f10.4,i10)')' NNerrorQtrain ',element(i1),&
              i2-1,qgrid(i2-1),' < error < ',qgrid(i2),countq(i2,i1)
          enddo
          write(ounit,'(a,f8.4,a,i6)')' Points with charge error above ',&
            qgrid(nstepsq),' in training set:',aboveq(i1)
          write(ounit,'(a,f10.4,a)')' Maximum charge error in training set:     ',errormaxq(i1),' e'
          write(ounit,*)'-------------------------------------------------------------'
        enddo ! i1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        aboveq(:)   =0
        errormaxq(:)=0.0d0
        countq(:,:) =0
!! get the file name for the testing data
        filename='testcharges.000000.out'
        if(nepochs.gt.999)then
          write(filename(15:18),'(i4)')nepochs
        elseif(nepochs.gt.99)then
          write(filename(16:18),'(i3)')nepochs
        elseif(nepochs.gt.9)then
          write(filename(17:18),'(i2)')nepochs
        else
          write(filename(18:18),'(i1)')nepochs
        endif
!!
        open(testqxunit,file=filename,form='formatted',status='old')
          rewind(testqxunit)
          zelemtest(:,:)=0
          read(testqxunit,*)dummy
 24       continue
          read(testqxunit,*,END=25)istruct,iatom
          backspace(testqxunit)
          read(testqxunit,*)idummy,idummy,&
            zelemtest(istruct,iatom),&
            qtestref(iatom,istruct),qtestnn(iatom,istruct)
          goto 24
 25       continue
        close(testqxunit)
!! analyze data
        do i1=1,ntest
          do i2=1,max_num_atoms
            if(zelemtest(i1,i2).gt.0)then ! check if atom exists
              error=abs(qtestref(i2,i1)-qtestnn(i2,i1))
              errormaxq(elementindex(zelemtest(i1,i2)))=&
                max(errormaxq(elementindex(zelemtest(i1,i2))),error)
              do i3=2,nstepsq
                if(error.lt.qgrid(i3))then
                  countq(i3,elementindex(zelemtest(i1,i2)))=&
                    countq(i3,elementindex(zelemtest(i1,i2)))+1
                  goto 26
                endif
              enddo ! i3
              aboveq(elementindex(zelemtest(i1,i2)))=&
                aboveq(elementindex(zelemtest(i1,i2)))+1
 26           continue
            endif
          enddo ! i2
        enddo ! i1
!! write results
        do i1=1,nelem
          write(ounit,*)'Test set charge error analysis for element ',element(i1),':'
          write(ounit,*)'                                   e                   Points'
          do i2=2,nstepsq
            write(ounit,'(a,a2,i5,f10.4,a11,f10.4,i10)')' NNerrorQtest ',element(i1),&
              i2-1,qgrid(i2-1),' < error < ',qgrid(i2),countq(i2,i1)
          enddo
          write(ounit,'(a,f10.4,a,i6)')' Points with charge error above ',&
            qgrid(nstepsq),' in test set: ',aboveq(i1)
          write(ounit,'(a,f10.4,a)')' Maximum charge error in test set:        ',errormaxq(i1),' e'
          write(ounit,*)'-------------------------------------------------------------'
        enddo ! i1
      endif
!!
      deallocate(counte) 
      deallocate(countf) 
      deallocate(countq) 
      deallocate(egrid)
      deallocate(fgrid)
      deallocate(qgrid)
!!
      return
      end
