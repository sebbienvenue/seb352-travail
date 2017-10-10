!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose:
!! determine num_basis, basis

!! called by:
!! - readinput.f90
!!
      subroutine readbasis()
!!
      use fileunits
      use globaloptions
      use basismod
!!
      implicit none
!!
      integer ztemp                                     ! internal
      integer itemp                                     ! internal
      integer icount                                    ! internal
      integer bcount                                    ! internal
      integer i1,i2                                     ! internal
!!
      character*40 keyword                              ! internal
      character*40 dummy                                ! internal
      character*80 line                                 ! internal
      character*80 elementline(nelem)                   ! internal
      character*2 elementtemp                           ! internal
      character*2 cbasis(50)                            ! internal
      character*2 ctemp                                 ! internal
!!
      logical lfound(nelem)                             ! internal
!!
!! initializations
      lfound(:)=.false.
      maxnum_basis=0
!!
      allocate(num_basis(nelem))
      num_basis(:)=0
!!
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
!!
 10   continue
!! read first 4 comment lines
      read(nnunit,*,END=20) keyword 
!!
      if(keyword.eq.'basis')then
        backspace(nnunit)
        read(nnunit,'(a80)',ERR=99)line
!!        write(ounit,*)'line ',line
        read(line,*)dummy,elementtemp
        call nuccharge(elementtemp,ztemp)
!!        write(ounit,*)'found element ',elementtemp,ztemp
!! check if basis for this element has already been specified
        if(lfound(elementindex(ztemp)))then
          write(ounit,'(a,a,i3)')'ERROR: basis specified twice for element ',elementtemp,ztemp
          stop !'
        endif
        lfound(elementindex(ztemp))=.true.
!! store basis information for this line for later processing
        elementline(elementindex(ztemp))=line
!! determine the number of basis function blocks
        icount=-1
 30     continue
        icount=icount+1
        read(line,*,END=31)dummy,elementtemp,(cbasis(i1),i1=1,icount)
        goto 30
 31     continue
        icount=icount-1
!        write(ounit,*)'Number of basis blocks found ',icount,(cbasis(i1),i1=1,icount)
!! count basis for this element:
        do i1=1,icount  
          ctemp=cbasis(i1)
!          write(ounit,'(i4,a2,x,a2,x,a2)')i1,cbasis(i1),ctemp(1:1),ctemp(2:2)     
          if(ctemp(2:2).eq.'s')then
            num_basis(elementindex(ztemp))=num_basis(elementindex(ztemp))+1
          elseif(ctemp(2:2).eq.'p')then
            num_basis(elementindex(ztemp))=num_basis(elementindex(ztemp))+3
          elseif(ctemp(2:2).eq.'d')then
            num_basis(elementindex(ztemp))=num_basis(elementindex(ztemp))+5
          elseif(ctemp(2:2).eq.'f')then
            num_basis(elementindex(ztemp))=num_basis(elementindex(ztemp))+7
          elseif(ctemp(2:2).eq.'g')then
            num_basis(elementindex(ztemp))=num_basis(elementindex(ztemp))+9
          else
            write(ounit,*)'ERROR: unknown basis function ',cbasis(i1)
            stop
          endif          
        enddo
!!        write(ounit,*)'Number of basis functions for ',elementtemp,num_basis(elementindex(ztemp))
!!
      endif
      goto 10
 20   continue
      close(nnunit)
!!
!! check if basis for all elements has been found
      do i1=1,nelem
        if(lfound(i1).eqv..false.)then
          write(ounit,*)'ERROR: basis is missing for element ',element(i1)
          stop
        endif
      enddo
!!
!! determine maxnum_basis
      do i1=1,nelem
        if(num_basis(i1).eq.0)then
          write(ounit,*)'ERROR: zero basis functions for element ',element(i1)
          stop
        endif
        maxnum_basis=max(maxnum_basis,num_basis(i1))
      enddo
      maxsize_subblock=maxnum_basis*maxnum_basis
!!
!! analyze basis in detail 
      allocate(basis(nelem,maxnum_basis,3))
      basis(:,:,:)=0
      do i1=1,nelem
        icount=-1
        bcount=0
 40     continue
        icount=icount+1
        line=elementline(i1)
        read(line,*,END=41)dummy,elementtemp,(cbasis(i2),i2=1,icount)
        goto 40
 41     continue
        icount=icount-1
!! count basis for this element:
        do i2=1,icount
          ctemp=cbasis(i2)
          read(ctemp(1:1),*)itemp
          if(ctemp(2:2).eq.'s')then
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=0
            basis(i1,bcount,3)=0
          elseif(ctemp(2:2).eq.'p')then
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=1
            basis(i1,bcount,3)=-1
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=1
            basis(i1,bcount,3)=0
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=1
            basis(i1,bcount,3)=1
          elseif(ctemp(2:2).eq.'d')then
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=2
            basis(i1,bcount,3)=-2
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=2
            basis(i1,bcount,3)=-1
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=2
            basis(i1,bcount,3)=0
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=2
            basis(i1,bcount,3)=1
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=2
            basis(i1,bcount,3)=2
          elseif(ctemp(2:2).eq.'f')then
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=3
            basis(i1,bcount,3)=-3
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=3
            basis(i1,bcount,3)=-2
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=3
            basis(i1,bcount,3)=-1
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=3
            basis(i1,bcount,3)=0
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=3
            basis(i1,bcount,3)=1
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=3
            basis(i1,bcount,3)=2
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=3
            basis(i1,bcount,3)=3
          elseif(ctemp(2:2).eq.'g')then
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=-4
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=-3
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=-2
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=-1
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=0
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=1
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=2
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=3
            bcount=bcount+1
            basis(i1,bcount,1)=itemp
            basis(i1,bcount,2)=4
            basis(i1,bcount,3)=4
          else
            write(ounit,*)'ERROR: unknown basis function ',cbasis(i1)
            stop
          endif
        enddo
      enddo
!!
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop

      end
