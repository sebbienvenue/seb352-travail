!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################


!! called by:
!!
      subroutine getkaldims_elec(&
        kaledim,kalcdim,corredim,corrcdim,&
        maxkaledim,maxcorredim,&
        num_weightsewaldfree)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer num_weightsewaldfree(nelem) ! in
      integer n_start           ! internal, just dummy here
      integer n_end             ! internal, just dummy here
!! Kalman matrix dimensions:
      integer corredim(nelem)   ! out
      integer corrcdim          ! out
      integer kaledim(nelem)    ! out
      integer kalcdim           ! out
      integer i1
      integer isum              ! internal
      integer maxkaledim        ! out
      integer maxcorredim       ! out
!!
!!
!! initializations
      isum=0
!!
!! non-parallel case
      if((mpisize.eq.1).or.lompmkl)then
!!
        do i1=1,nelem
          if(lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))then
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
          endif ! lelec
        enddo ! i1
!!
!! get dimensions for charge constraint
        if(lchargeconstraint.and.lelec)then
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
        do i1=1,nelem
          if(lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))then
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
          endif ! lelec
        enddo ! i1

        if(lchargeconstraint.and.lelec)then
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
      maxkaledim=0
      maxcorredim=0
      do i1=1,nelem
        maxkaledim=max(maxkaledim,kaledim(i1))
        maxcorredim=max(maxcorredim,corredim(i1))
      enddo
!!
      return
      end      
