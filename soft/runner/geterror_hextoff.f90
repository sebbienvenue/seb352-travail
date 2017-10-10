!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine geterror_hextoff(iswitch,countepoch,ntrain,&
        imaxerror_hextoff,minvalue_hextoff,&
        maxvalue_hextoff,avvalue_hextoff,&
        rmse_hextoff,mad_hextoff,maxerror_hextoff)
!!
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use symfunctions
      use nnham
      use structures
      use basismod
!!
      implicit none
!!
      integer i1,i2,i3                                     ! internal
      integer iswitch                                      ! in 0=train, 1=test
      integer nforces_hextoff                              ! internal
      integer ntrain                                       ! in
      integer ndone                                        ! internal
      integer ndonepara                                    ! internal
      integer ncount                                       ! internal
      integer npoints                                      ! internal
      integer nhextoff                                     ! internal
      integer nstruct                                      ! internal
      integer n_start                                      ! internal
      integer n_end                                        ! internal
      integer imaxerror_hextoff                            ! out
      integer, dimension(:), allocatable :: imaxerror_temp ! internal
      integer icount                                       ! internal
      integer countepoch                                   ! in
      integer nenergies                                    ! internal
      integer tempunit                                     ! internal
      integer pxunit                                       ! internal
      integer fxunit                                       ! internal
      integer qxunit                                       ! internal
      integer matrixsize                                   ! internal
!!
      real*8, dimension(:), allocatable :: maxerror_temp           ! internal
      real*8 minvalue_hextoff(1,maxnum_funcvalues_hextoff)       ! in
      real*8 maxvalue_hextoff(1,maxnum_funcvalues_hextoff)       ! in
      real*8 avvalue_hextoff(1,maxnum_funcvalues_hextoff)        ! in
      real*8 rmse_hextoff                             ! out
      real*8 mad_hextoff                              ! out
      real*8 maxerror_hextoff                         ! out
      real*8 nnhextoff_list(nblock,&
               num_basis(elementindex(hextoff_training_triplet(1))),&
               num_basis(elementindex(hextoff_training_triplet(2)))) ! internal
!      real*8 nnhextoff_listtemp(npoints,&
!               num_basis(elementindex(hextoff_training_triplet(1))),&
!               num_basis(elementindex(hextoff_training_triplet(2)))) ! internal
      real*8 edummy                                   ! internal
      real*8 forcesum(3)                              ! internal
      real*8 maxforce_dummy                           ! internal
!!
      character*40 filename                           ! internal
!!
!!
!!==============================================
!! initializations
!!==============================================
!! counters
      nhextoff                 = 0
      ncount                   = ntrain
      ndone                    = 0
      n_start                  = 1
      ndonepara                = 0     
!! data
      nnhextoff_list(:,:,:)   = 0.0d0
!!      nnhextoff_listtemp(:,:,:) = 0.0d0
!! others
      edummy                   = 1.d12
      matrixsize               = num_basis(elementindex(hextoff_training_triplet(1)))*&
                         num_basis(elementindex(hextoff_training_triplet(2)))

!!
!!==============================================
!! open files and write headers
!!==============================================
        if(iswitch.eq.0)then !! We need to FIXME the format statements of the output of the testpoint/trainpoint info
          write(filename,'(A,I3.3,A,I3.3,A,I3.3,A)')&
            'function_hextoff.',hextoff_training_triplet(1),'.',&
            hextoff_training_triplet(2),'.',&
            hextoff_training_triplet(3),'.data'
          open(symhextoffunit,file=filename,form='formatted',status='old')
          rewind(symhextoffunit) !'
          if(lwritetrainhextoff)then
            write(filename,'(A,I6.6,A)')&
              'trainhextoff.',countepoch,'.out'
            open(trainpxunit,file=filename,form='formatted',status='replace')
            rewind(trainpxunit) !'
            write(trainpxunit,'(a)')' structure element     hextoff(DFT)      hextoff(NN)  '
          endif
          open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
          rewind(trainstructunit) !'
        elseif(iswitch.eq.1)then
          write(filename,'(A,I3.3,A,I3.3,A,I3.3,A)')&
            'testing_hextoff.',hextoff_training_triplet(1),'.',&
            hextoff_training_triplet(2),'.',&
            hextoff_training_triplet(3),'.data'
          open(tymhextoffunit,file=filename,form='formatted',status='old')
          rewind(tymhextoffunit) !'
          if(lwritetrainhextoff)then
            write(filename,'(A,I6.6,A)')&
              'testhextoff.',countepoch,'.out'
            open(testpxunit,file=filename,form='formatted',status='replace')
            rewind(testpxunit) !'
            write(testpxunit,'(a)')' structure element     hextoff(DFT)      hextoff(NN)  '
          endif
          open(teststructunit,file='teststruct.data',form='formatted',status='old')
          rewind(teststructunit) !'
        endif ! iswitch
!!
!!==============================================
!! loop block-wise over all structures 
!!==============================================
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!!==============================================
!!==============================================
!! do all file reading here for nstruct structures at one place to allow for parallelization
!!==============================================
!!==============================================
!!==============================================
!! read npoint short range symmetry function sets (train or test)
!!==============================================
        if(iswitch.eq.0)then ! train
          tempunit=symhextoffunit          
        elseif(iswitch.eq.1)then ! test
          tempunit=tymhextoffunit          
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
        call readfunctions_hextoff(1,tempunit,npoints,1,&
          max_num_atoms,maxnum_funcvalues_hextoff,num_funcvalues_hextoff,&
          symfunction_hextoff_list)
!!
!!==============================================
!! read the structures needed for the calculation of the electrostatic energy
!! and get reference forces from DFT
!! must be called after readfunctions because it needs num_atoms_list
!!==============================================
        if(iswitch.eq.0)then ! train
          tempunit=trainstructunit          
        elseif(iswitch.eq.1)then ! test
          tempunit=teststructunit          
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
        call getstructures(tempunit,npoints)
!!
!!==============================================
!!==============================================
!! end of file reading 
!!==============================================
!!==============================================
!!
!!==============================================
!! fill local arrays for structures n_start to n_end
!!==============================================
        icount=0
!!        do i1=n_start,n_end
!        num_atoms_mpi(icount)       = num_atoms_list(i1)
!!        hextoff_listtemp(icount,:,:)     = hextoff_list(i1,:,:)
!        xyzstruct_mpi(:,:,icount)   = xyzstruct_list(:,:,i1)
!!        enddo
!!
!!==============================================
!! get the energies, forces and charges 
!!==============================================
         
        ndonepara=ndone+n_start-1 !XXX
!! FIXME: getatomicoutput_para is not completed yet
        call gethextoffoutput(npoints,ndonepara,&
          nhextoff,imaxerror_hextoff,&
          minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
          rmse_hextoff,mad_hextoff,maxerror_hextoff,&
          nnhextoff_list)
!!
        ndone = ndone + npoints
!!============================================================
!! if there are structures left go to next block of structures
!!============================================================
        if(ncount.gt.0) goto 10
!============================================================
!! end of block wise loop
!!============================================================
!!
!!============================================================
!! close files
!!============================================================

!!
        if(iswitch.eq.0)then
          close(symunit)
          close(trainstructunit)
          if(lwritetrainhextoff)                close(trainpxunit)
        elseif(iswitch.eq.1)then
          close(tymunit)
          close(teststructunit)
          if(lwritetrainhextoff)                close(testpxunit)
        else
          write(ounit,*)'ERROR: unknown iswitch in geterror ',iswitch
          stop
        endif
        
!!
!!==============================================
!! calculate the final RMSEs
!!==============================================

      matrixsize = num_basis(elementindex(hextoff_training_triplet(1)))*&
                   num_basis(elementindex(hextoff_training_triplet(2)))

      call getrmse_hextoff(matrixsize,ndone,&
         rmse_hextoff,mad_hextoff)
!!
!      return
      end
