!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! MG: Added: Subroutine to analyze average trace of correlation matrix
!! MG: Added: and print information to output file

!! called by:
!! - fitting_short_atomic.f90
!! 
      subroutine analyzecorrmat(maxcorrdim,corrdim,corrmatrix_list,num_weights_short_atomic)
!!
      use globaloptions
      use fileunits
!!
      implicit none
!!
      integer i1,i2
      integer num_weights_short_atomic(nelem)  ! in
      integer maxcorrdim                       ! in
      integer corrdim(nelem)                   ! in
      integer currpos                          ! internal
      real*8 trace_sum                         ! internal
      real*8 corrmatrix_list(maxcorrdim,nelem) ! in
!!
      do i1=1,nelem
        trace_sum=0d0
        currpos=1
        !! MG: Sum over the elements in the packed triangular matrix
        !! MG: corresponding to the diagonal positions.
        do i2=1,num_weights_short_atomic(i1)
          trace_sum=trace_sum+corrmatrix_list(currpos,i1)
          currpos=currpos+num_weights_short_atomic(i1)-i2+1
        enddo ! i2
        trace_sum=trace_sum/dble(num_weights_short_atomic(i1)) 
        !!write(ounit,'(a10,x,a2,x,1f20.2)') ' CORRMAT  ',element(i1),trace_sum 
!!        write(ounit,'(a10,x,a2,x,1e20.5e3)') ' CORRMAT  ',element(i1),trace_sum 
      enddo ! i1
      return
      end
