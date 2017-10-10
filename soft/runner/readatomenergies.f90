!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90
!!
      subroutine readatomenergies()
!!
      use fileunits
      use globaloptions 
!! 
      implicit none
!!
      integer i1
      integer icount  
      integer ztemp
      integer zelem(nelem)                              ! internal
!!
      real*8 etemp
!!
      character*2 elementtemp                           ! internal
      character*40 dummy                                ! internal
      character*2 elementsymbol(nelem)                  ! internal
      character*40 keyword                              ! internal
!!
      logical lfound(nelem)                             ! internal
!!
        open(nnunit,file='input.nn',form='formatted',status='old')
        rewind(nnunit)
        icount=0
 31     continue       
          read(nnunit,*,END=32)keyword
          if(keyword.eq.'atom_energy')then
            backspace(nnunit)
            read(nnunit,*,err=99)dummy,elementtemp
            call nuccharge(elementtemp,ztemp)
            do i1=1,nelem
              if(ztemp.eq.nucelem(i1))then
                backspace(nnunit)
                icount=icount+1
                read(nnunit,*,err=99)dummy,elementsymbol(icount),atomrefenergies(icount)
                goto 31
              endif
            enddo
            write(ounit,*)'### WARNING ### atom_energy for element ',elementtemp,' is ignored'
          endif
          goto 31
 32     continue      
        close(nnunit)
!! get the nuclear charges
        do i1=1,icount
          call nuccharge(elementsymbol(i1),zelem(i1))
        enddo
!! sort according to nuclear charge
        if(nelem.gt.1)then
 30       continue
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
              goto 30
            endif
          enddo
        endif
!! check if all elements have been found
        lfound(:)=.false.
        do i1=1,nelem
          lfound(elementindex(zelem(i1)))=.true.
        enddo
        do i1=1,nelem
          if(lfound(i1).eqv..false.)then
            write(ounit,*)'Error: atom_energy not found for element ',nucelem(i1)
            stop
          endif
        enddo
!! write reference energies to runner.out 
        write(ounit,*)'atomic reference energies read from input.nn:'
        do i1=1,nelem
          write(ounit,'(a1,a2,x,f18.8)')' ',elementsymbol(i1),atomrefenergies(i1)
        enddo
        write(ounit,*)'-------------------------------------------------------------'
!!
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop

      end

