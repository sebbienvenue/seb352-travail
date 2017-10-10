!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - ewaldenergies_para.f90
!! - geterror.f90
!! - geterrorpair.f90
!! - getshortenergies_para.f90
!! - getshortenergies_parapair.f90
!!
      subroutine calcrmse_hextoff(matrixsize,npoints,ndone,&
          imaxerror_local,rmse_local,mad_local,&
          maxerror_local,hextoff_local,nnhextoff_local)
!!
      use fileunits
      use globaloptions
      use fittingoptions
      use basismod
!!
      implicit none
!!
      integer matrixsize                  ! in
      integer npoints                     ! in
      integer ndone                       ! in
      integer i1,i2, i3                   ! internal
      integer imaxerror_local             ! in/out
!!
      real*8 rmse_local                   ! in/out
      real*8 mad_local                    ! in/out
      real*8 hextoff_local(nblock,maxnum_basis,maxnum_basis)   ! in
      real*8 nnhextoff_local(nblock,&
                num_basis(elementindex(hextoff_training_triplet(1))),&
                num_basis(elementindex(hextoff_training_triplet(2))))   ! in
!!      real*8 maxenergy_local              ! in
      real*8 maxerror_local               ! in/out
!!
!!
      do i1=1,npoints
        do i2=1,num_basis(elementindex(hextoff_training_triplet(1)))
          do i3=1,num_basis(elementindex(hextoff_training_triplet(2)))
!!           write(*,*) i1, hextoff_local(i1,i2,i3), nnhextoff_local(i1,i2,i3)
            rmse_local= rmse_local +(hextoff_local(i1,i2,i3)-nnhextoff_local(i1,i2,i3))**2
            mad_local = mad_local  +abs(hextoff_local(i1,i2,i3)-nnhextoff_local(i1,i2,i3))
            if(abs(hextoff_local(i1,i2,i3)-nnhextoff_local(i1,i2,i3)).gt.maxerror_local)then
              maxerror_local = abs(hextoff_local(i1,i2,i3)-nnhextoff_local(i1,i2,i3)) 
              imaxerror_local= ndone + i1
            endif
          enddo
        enddo
      enddo
!!
      return
      end
