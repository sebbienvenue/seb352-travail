!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - prediction.f90
!!
      subroutine enforcecharge(nntotalcharge,nnatomcharge)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1                          ! internal
!!
      real*8 nntotalcharge                ! in
      real*8 nnatomcharge(max_num_atoms)  ! in/out
      real*8 sumpos                       ! internal
      real*8 sumneg                       ! internal
      real*8 scalefactorpos               ! internal
      real*8 scalefactorneg               ! internal
      real*8 chargeerror                  ! internal
      real*8 checksum                     ! internal
!!
!! initializations
      sumpos=0.0d0
      sumneg=0.0d0
      scalefactorpos=0.0d0
      scalefactorneg=0.0d0
      chargeerror=0.0d0
      checksum=0.0d0
!!
      write(ounit,'(a)')' Charge rescaling enforced for neutrality'
      write(ounit,'(a,f14.8)')' Total system charge before rescaling ',nntotalcharge
      write(ounit,'(a,f14.8)')' Average charge modification per atom ',-1.0d0*nntotalcharge/dble(max_num_atoms)
!!
!! calculate sum of positive and negative charges
      do i1=1,max_num_atoms
        if(nnatomcharge(i1).gt.0.0d0)then
          sumpos=sumpos+nnatomcharge(i1)
        else
          sumneg=sumneg+nnatomcharge(i1)
        endif
      enddo
!!      write(ounit,'(a,f14.8)')'Sum of positive charges ',sumpos
!!      write(ounit,'(a,f14.8)')'Sum of negative charges ',sumneg
!!      write(ounit,'(a,f14.8)')'Sum ',sumneg+sumpos
!!
!! check
      if(abs(nntotalcharge-(sumpos+sumneg)).gt.0.000001d0)then
        write(ounit,*)'Error: total charge not correct in enforcetotcharge'
        stop !'
      endif
!!
      chargeerror=nntotalcharge/2.0d0
!!
!! rescale charges
      scalefactorneg=(sumneg-chargeerror)/sumneg
      scalefactorpos=(sumpos-chargeerror)/sumpos
!!      write(ounit,'(a,2f14.8)')'scalefactors ',scalefactorneg,scalefactorpos
      do i1=1,max_num_atoms
        if(nnatomcharge(i1).gt.0.0d0)then
          nnatomcharge(i1)=nnatomcharge(i1)*scalefactorpos
        else
          nnatomcharge(i1)=nnatomcharge(i1)*scalefactorneg
        endif
      enddo      
!!
!! final check 
      do i1=1,max_num_atoms
        checksum=checksum+nnatomcharge(i1)
      enddo
      if(abs(checksum).gt.0.0000001d0)then
        write(ounit,*)'Error: charge rescaling did not work!'
        stop
      endif
!!
!! set new total charge
      nntotalcharge=0.0d0
      do i1=1,max_num_atoms
        nntotalcharge=nntotalcharge+nnatomcharge(i1)
      enddo
      write(ounit,'(a,f14.8)')' Total system charge after rescaling  ',nntotalcharge
!!
      write(ounit,*)'-------------------------------------------------------------'
      return
      end
