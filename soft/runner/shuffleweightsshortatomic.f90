!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:

      subroutine shuffleweightsshortatomic(sseed)

      use fileunits
      use globaloptions
      use fittingoptions
      use nnshort_atomic

      implicit none

      integer num_shuffled
      integer i1,i2
      integer sseed

      real*8 ran0
      real*8 z 


      num_shuffled=0

!! loop over all weights
      do i1=1,nelem
        do i2=1,num_weights_short_atomic(i1)
!! get random number
          z=ran0(sseed)
          if(z.lt.shuffle_weights_short_atomic)then
            weights_short_atomic(i2,i1)=weights_min+(weights_max-weights_min)*ran0(sseed)
            num_shuffled=num_shuffled+1
          endif
        enddo
      enddo

      write(ounit,'(a,i5,a)')' Shuffling ',num_shuffled,' weights'
      write(ounit,*)'-------------------------------------------------------------------------------'

      return
      end
