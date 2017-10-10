!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - optimize_ewald.f90
!!
      subroutine sortforceerror(npoints,forceerror_list,lusef)
!!
      use fileunits
      use fittingoptions
      use globaloptions
      use structures
!!
      implicit none
!!
      integer npoints
      integer i1,i2,i3
      integer icount
      integer ixyz(3*nblock*max_num_atoms)                                   ! internal 
      integer iatom(3*nblock*max_num_atoms)                                  ! internal 
      integer istruct(3*nblock*max_num_atoms)                                ! internal 
!!
      real*8 forceerror_list(3,max_num_atoms,nblock)                         ! in
      real*8 forceerror_temp(3*max_num_atoms*nblock)                         ! internal
      real*8 forceerror_copy(3*max_num_atoms*nblock)                         ! internal
      real*8 temp                                                            ! internal
      real*8 fthres                                                          ! out
      real*8 select1                                                         ! internal
!!
      logical lusef(3,max_num_atoms,nblock)                                  ! out 
!!
!! initializations
      lusef(:,:,:)=.false.
!!
!! transform array
      forceerror_temp(:)=0
      icount=0
      do i1=1,npoints
        do i2=1,num_atoms_list(i1)
          do i3=1,3
            if(lupdatebyelement.and.(elemupdate.eq.zelem_list(i1,i2)))then
              icount=icount+1
              forceerror_temp(icount)=forceerror_list(i3,i2,i1)
              iatom(icount)=i2
              istruct(icount)=i1
              ixyz(icount)=i3
            elseif(.not.lupdatebyelement)then
              icount=icount+1
              forceerror_temp(icount)=forceerror_list(i3,i2,i1)
              iatom(icount)=i2
              istruct(icount)=i1
              ixyz(icount)=i3
            endif
          enddo
        enddo
      enddo
!!
!! determine the threshold for update
      temp=(1.d0-worstf)*dble(icount)
      i2=int(temp)
      forceerror_copy(:)=forceerror_temp(:)
!! caution: select changes order of forceerror_copy
      fthres=select1(i2,icount,forceerror_copy)
!!      write(ounit,*)' select found ',i1,forceerror_copy(i1)
!!
!! set luseq array
      do i1=1,icount
        if(forceerror_temp(i1).gt.fthres)then
!!          write(ounit,*)i1,forceerror_temp(i1),istruct(i1),iatom(i1),ixyz(i1)
          lusef(ixyz(i1),iatom(i1),istruct(i1))=.true.
        endif
      enddo
!!
!! debug
!!      write(ounit,*)'lusef array'
!!      do i1=1,npoints
!!        do i2=1,num_atoms_list(i1)
!!          do i3=1,3
!!            write(ounit,'(4i6,f14.8,l)')i1,i2,i3,zelem_list(i1,i2),&
!!              forceerror_list(i3,i2,i1),lusef(i3,i2,i1)
!!          enddo
!!        enddo
!!      enddo
!!
!!
!!      stop
!!
      return
      end
