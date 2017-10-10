!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate nneshort_mpi and nnshortforces_mpi for nstruct structures

!! called by:
!!
      subroutine gethextoffoutput(npoints,ndone,&
        nhextoff,imaxerror_hextoff,&
        minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
        rmse_hextoff,mad_hextoff,maxerror_hextoff,&
        nnhextoff_list)
!!
      use fileunits
      use fittingoptions
      use nnflags 
      use globaloptions
      use structures
      use symfunctions
      use nnham
      use basismod
!!
      implicit none
!!
      integer nhextoff                                   ! in/out
      integer npoints                                    ! in
      integer ndone                                      ! in
      integer imaxerror_hextoff                         ! in/out
      integer ndummy                                     ! internal
      integer matrixsize                                 ! internal
!! in
      real*8 minvalue_hextoff(1,maxnum_funcvalues_hextoff)           ! in
      real*8 maxvalue_hextoff(1,maxnum_funcvalues_hextoff)           ! in
      real*8 avvalue_hextoff(1,maxnum_funcvalues_hextoff)            ! in
      real*8 nnhextoff_list(nblock,&
               num_basis(elementindex(hextoff_training_triplet(1))),&
               num_basis(elementindex(hextoff_training_triplet(2)))) ! out
!! errors
      real*8 rmse_hextoff                                ! in/out
      real*8 rmse_hextoffcombo                           ! in/out
      real*8 mad_hextoff                                 ! in/out
!! symmetry function parameters
      real*8 maxerror_hextoff                           ! in/out
      real*8 edummy                                      ! internal
!!
!!
!!
!!========================================================
!! scale symmetry functions 
!!========================================================
      call scalesym_hextoff(npoints,&
        maxnum_funcvalues_hextoff,num_funcvalues_hextoff,&
        symfunction_hextoff_list,&
        minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
        scmin_hextoff,scmax_hextoff)
!!
!!========================================================
!! predict the NN output for npoint data sets
!!========================================================
      matrixsize = num_basis(elementindex(hextoff_training_triplet(1)))*&
                   num_basis(elementindex(hextoff_training_triplet(2)))
      call gethextoff(matrixsize,npoints,& 
        symfunction_hextoff_list,nnhextoff_list)
!!
!!
!!========================================================
!! calculate rmse_short: in/out rmse_short
!!========================================================
       call calcrmse_hextoff(matrixsize,npoints,ndone,&
         imaxerror_hextoff,rmse_hextoff,mad_hextoff,&
         maxerror_hextoff,hextoff_list,nnhextoff_list)
!!
!!
      return
      end
