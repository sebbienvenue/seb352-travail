!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getsymmetryfunctions.f90
!! - prediction.f90
!!
      subroutine getrefatoms(ounit,iunit,nelem,atomrefenergies,&
         lremoveatomenergies,ldebug)
!!
      implicit none
!!
      integer iunit                      ! in
      integer ounit                      ! in
      integer nelem                      ! in
      integer zelem(nelem)               ! internal
      integer ztemp                      ! internal
      integer i1,i2
!!
      real*8 atomrefenergies(nelem)      ! out
      real*8 etemp                       ! internal
!!
      character*2 elementsymbol(nelem)
      character*2 elementtemp            ! internal
!!
      logical lremoveatomenergies        ! in
      logical ldebug                     ! in
      logical lexist                     ! internal
!!
!! initialization
      atomrefenergies(:)=0.0d0
!!
      if(lremoveatomenergies)then
!! first check if input file is present
        inquire(file='atomenergies.in',exist=lexist)
        if(.not.lexist)then
          write(ounit,*)'Error: file atomenergies.in not found'
          stop
        endif
        open(iunit,file='atomenergies.in',form='formatted',status='old')
        rewind(iunit)
        do i1=1,nelem 
          read(iunit,*)elementsymbol(i1),atomrefenergies(i1)
        enddo ! i1
        close(iunit)
!! get the nuclear charges 
        do i1=1,nelem
          call nuccharge(elementsymbol(i1),zelem(i1))
        enddo
!! sort according to nuclear charge
        if(nelem.gt.1)then
 10       continue
          do i1=1,nelem-1
            if(zelem(i1).gt.zelem(i1+1))then
              ztemp=zelem(i1)
              elementtemp=elementsymbol(i1)
              etemp=atomrefenergies(i1)
              zelem(i1)=zelem(i1+1)
              elementsymbol(i1)=elementsymbol(i1+1)
              atomrefenergies(i1)=atomrefenergies(i1+1)
              zelem(i1+1)=ztemp
              elementsymbol(i1+1)=elementtemp
              atomrefenergies(i1+1)=etemp
              goto 10
            endif
          enddo
        endif
!!
      endif
!!
!! debugging
      if(ldebug.and.lremoveatomenergies)then
        write(ounit,*)'atomic reference energies read from file:'
        do i1=1,nelem
          write(ounit,'(a1,a2,x,f18.8)')' ',elementsymbol(i1),atomrefenergies(i1)
        enddo
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      return
      end
