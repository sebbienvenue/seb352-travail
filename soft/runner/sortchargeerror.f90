!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - optimize_ewald.f90
!!
      subroutine sortchargeerror(npoints,&
        chargeerror_list,luseq)
!!
      use fileunits
      use fittingoptions
      use globaloptions
      use structures
!!
      implicit none
!!
      integer npoints
      integer i1,i2
      integer icount
      integer iatom(nblock*max_num_atoms)                                    ! internal 
      integer istruct(nblock*max_num_atoms)                                  ! internal 
!!
      real*8 chargeerror_list(nblock,max_num_atoms)                          ! in
      real*8 chargeerror_temp(nblock*max_num_atoms)                          ! internal
      real*8 chargeerror_copy(nblock*max_num_atoms)                          ! internal
      real*8 temp                                                            ! internal
      real*8 qthres                                                          ! out
      real*8 select1                                                         ! internal
!!
      logical luseq(nblock,max_num_atoms)                                    ! out 
!!
!! initializations
      luseq(:,:)=.false.
!!
!! transform array
      chargeerror_temp(:)=0
      icount=0
      do i1=1,npoints
        do i2=1,num_atoms_list(i1)
          if(lupdatebyelement.and.(elemupdate.eq.zelem_list(i1,i2)))then
            icount=icount+1
            chargeerror_temp(icount)=chargeerror_list(i1,i2)
            iatom(icount)=i2
            istruct(icount)=i1
          elseif(.not.lupdatebyelement)then
            icount=icount+1
            chargeerror_temp(icount)=chargeerror_list(i1,i2)
            iatom(icount)=i2
            istruct(icount)=i1
          endif
        enddo
      enddo
!!
!! debug
!!      write(ounit,*)'chargeerror_temp before'
!!      do i1=1,icount
!!        write(ounit,'(i6,f14.8,2i6)')i1,chargeerror_temp(i1),istruct(i1),iatom(i1)
!!      enddo
!!
!! determine the threshold for update
      temp=(1.d0-worstq)*dble(icount)
      i2=int(temp)
      chargeerror_copy(:)=chargeerror_temp(:)
!! caution: select changes order of chargeerror_copy
      qthres=select1(i2,icount,chargeerror_copy)
!!      write(ounit,*)' select found ',i1,chargeerror_copy(i1)
!!
!! set luseq array
      do i1=1,icount
        if(chargeerror_temp(i1).gt.qthres)then
!!          write(ounit,*)i1,chargeerror_temp(i1),istruct(i1),iatom(i1)
          luseq(istruct(i1),iatom(i1))=.true.
        endif
      enddo
!! debug
!!      write(ounit,*)'luseq array'
!!      do i1=1,npoints
!!        do i2=1,num_atoms_list(i1)
!!          write(ounit,'(3i6,f14.8,l)')i1,i2,zelem_list(i1,i2),&
!!            chargeerror_list(i1,i2),luseq(i1,i2)
!!        enddo
!!      enddo
!!
!!      stop
!!
      return
      end
