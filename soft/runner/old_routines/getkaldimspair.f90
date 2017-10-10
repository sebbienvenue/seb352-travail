!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine getkaldimspair(nelem,npairs,&
        kaldim,kaledim,kalcdim,corrdim,corrfdim,corredim,corrcdim,&
        maxkaldim,maxkaledim,maxcorrdim,maxcorrfdim,maxcorredim,&
        num_weightsewaldfree,num_weightspairfree,&
        lshort,lewald,lchargeconstraint)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
!!
      implicit none
!!
      integer nelem             ! in
      integer npairs            ! in
      integer n_start           ! internal, just dummy here
      integer n_end             ! internal, just dummy here
!! Kalman matrix dimensions:
      integer corrdim(npairs)   ! out
      integer corrfdim(npairs)  ! out
      integer corredim(nelem)   ! out
      integer corrcdim          ! out
      integer kaldim(npairs)    ! out
      integer kaledim(nelem)    ! out
      integer kalcdim           ! out
      integer i1
      integer isum              ! internal
      integer maxkaldim         ! out
      integer maxkaledim        ! out
      integer maxcorrdim        ! out
      integer maxcorrfdim       ! out
      integer maxcorredim       ! out
      integer num_weightsewaldfree(nelem) ! in
      integer num_weightspairfree(npairs)
!!
      logical lshort            ! in
      logical lewald            ! in
      logical lchargeconstraint ! in
!!
!! initializations
      isum=0
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)

!!
!! non-parallel case
      if(mpisize.eq.1)then
         do i1=1,npairs !nelem
           if(lshort)then
             if(optmodee.eq.1)then
               corrdim(i1)=num_weightspairfree(i1)*(num_weightspairfree(i1)+1)/2  
             else
               corrdim(i1)=1
             endif
             if(optmodef.eq.1)then
               corrfdim(i1)=num_weightspairfree(i1)*(num_weightspairfree(i1)+1)/2 
             else
               corrfdim(i1)=1
             endif
             if((optmodee.ne.1).and.(optmodef.ne.1))then
               kaldim(i1)=1
             else
               kaldim(i1)=num_weightspairfree(i1) 
             endif
           endif ! lshort
         enddo ! i1

!!
         do i1=1,nelem
           if(lewald.and.(optmodeq.eq.1))then
             corredim(i1)=num_weightsewaldfree(i1)*(num_weightsewaldfree(i1)+1)/2
             kaledim(i1)=num_weightsewaldfree(i1)
             if(lchargeconstraint)then
!!              corrcdim=nelem*num_weightsewaldfree*(nelem*num_weightsewaldfree+1)/2
                isum=isum+num_weightsewaldfree(i1)
!!              kalcdim=nelem*num_weightsewaldfree
             endif
           else
             corredim(i1)=1
             kaledim(i1)=1
           endif ! lewald
         enddo ! i1
!!
!! get dimensions for charge constraint
         if(lchargeconstraint.and.lewald)then
           corrcdim=isum*(isum+1)/2
           kalcdim =isum
           isum=0
         else
           corrcdim=1
           kalcdim=1
           isum=0
         endif
!!
      else ! parallel case
!!
         do i1=1, npairs !nelem
           if(lshort)then
             if(optmodee.eq.1)then
               call mpifitdistribution(num_weightspairfree(i1),kaldim(i1),n_start,n_end)
               corrdim(i1)=kaldim(i1)*num_weightspairfree(i1)
             else
               corrdim(i1)=1
             endif
             if(optmodef.eq.1)then
               call mpifitdistribution(num_weightspairfree(i1),kaldim(i1),n_start,n_end)
               corrfdim(i1)=kaldim(i1)*num_weightspairfree(i1)
             else
               corrfdim(i1)=1
             endif
             if((optmodee.ne.1).and.(optmodef.ne.1))then
               kaldim(i1) =1
             endif
           endif ! lshort
         enddo ! i1
!!
         do i1=1, nelem
           if(lewald.and.(optmodeq.eq.1))then
             call mpifitdistribution(num_weightsewaldfree(i1),kaledim(i1),n_start,n_end)
             corredim(i1)=kaledim(i1)*num_weightsewaldfree(i1)
!! this still needs to be split for the different processes 
             if(lchargeconstraint)then
!!              corrcdim=nelem*num_weightsewaldfree*(nelem*num_weightsewaldfree+1)/2
                isum=isum+num_weightsewaldfree(i1)
!!              kalcdim=nelem*num_weightsewaldfree
             endif
           else
             kaledim(i1)=1
             corredim(i1)=1
           endif ! lewald
         enddo ! i1

         if(lchargeconstraint.and.lewald)then
           corrcdim=isum*(isum+1)/2
           kalcdim =isum
           isum=0
         else
           corrcdim=1
           kalcdim=1
           isum=0
         endif
       endif ! mpisize
!!
       maxkaldim=0
       maxkaledim=0
       maxcorrdim=0
       maxcorrfdim=0
       maxcorredim=0
       do i1=1,npairs
         maxkaldim=max(maxkaldim,kaldim(i1))
         maxcorrdim=max(maxcorrdim,corrdim(i1))
         maxcorrfdim=max(maxcorrfdim,corrfdim(i1))
       enddo

       do i1=1,nelem
         maxkaledim=max(maxkaledim,kaledim(i1))
         maxcorredim=max(maxcorredim,corredim(i1))
       enddo
!!
!! debug
!!      write(ounit,*)'kaldim   ',kaldim
!!      write(ounit,*)'kaledim  ',kaledim
!!      write(ounit,*)'kalcdim  ',kalcdim
!!      write(ounit,*)'corrdim  ',corrdim
!!      write(ounit,*)'corredim ',corredim
!!      write(ounit,*)'corrcdim ',corrcdim
!!
!!     
      return
      end      
