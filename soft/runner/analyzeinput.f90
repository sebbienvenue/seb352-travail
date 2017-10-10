!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine analyzeinput(iswitch)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer idummy                                              ! internal
      integer num_atoms                                           ! internal
      integer i1,i2,i3,i4,i5                                      ! internal
      integer j1,j2,j3,j4,j5                                      ! internal
      integer nstruct                                             ! internal
      integer, dimension(:), allocatable :: zelem                 ! internal 
      integer num_atomselement(nelem)                             ! internal
      integer dim1,dim2,dim3,dim4,dim5                            ! internal
      integer memsize                                             ! internal
      integer memthres                                            ! internal
      integer, dimension(:,:,:,:,:,:), allocatable :: temparray   ! internal 
      integer, dimension(:,:,:,:,:,:), allocatable :: temparray2  ! internal 
      integer iswitch                                             ! in
      integer itemp
      integer num_pairs                                           ! internal
      integer unit1                                               ! internal
      integer unit2                                               ! internal
      integer nsystemsize(max_num_atoms,2)                        ! internal
!!
      real*8 zdummy                                               ! internal
!!
      character*14 ctemp
!!
      logical lperiodic                                           ! internal
!!
!! initializations
      nstruct=0
      nsystemsize(:,:)=0
!! FIXME: We need a more clever way to dimension the array temparray
!! It needs too much memory, and not all combinations of element numbers will be present
!! (maybe introduce pointer array?)
      dim1=max_num_atoms
      dim2=max_num_atoms
      dim3=max_num_atoms
      dim4=max_num_atoms
      dim5=max_num_atoms
      if(nelem.lt.5)dim5=0
      if(nelem.lt.4)dim4=0
      if(nelem.lt.3)dim3=0
      if(nelem.lt.2)dim2=0

      memsize=(dim1+1)*(dim2+1)*(dim3+1)*(dim4+1)*(dim5+1)*2
      memthres=1024*1024*1024  ! roughly 1 GB

      if(memsize.gt.memthres)then
        write(ounit,*)'### WARNING ### : analyzeinput.f90: Analysis of input.data cannot be done (insufficient memory) ',memsize
        return !'
      endif

!! this is extremely memory consuming and needs to be changed
      allocate(temparray(0:dim1,0:dim2,0:dim3,0:dim4,0:dim5,2))
      temparray(:,:,:,:,:,:)=0
!!
!! allocate array for summary
      allocate(temparray2(2,2,2,2,2,2))
      temparray2(:,:,:,:,:,:)=0
!!
      if(iswitch.eq.0)then
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
        open(symunit,file='function.data',form='formatted',status='old')
        rewind(symunit)
        unit1=trainstructunit
        unit2=symunit
      elseif(iswitch.eq.1)then
        open(teststructunit,file='teststruct.data',form='formatted',status='old')
        rewind(teststructunit)
        open(tymunit,file='testing.data',form='formatted',status='old')
        rewind(tymunit)
        unit1=teststructunit
        unit2=tymunit
      elseif(iswitch.eq.2)then
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
        open(symeunit,file='functione.data',form='formatted',status='old')
        rewind(symeunit)
        unit1=trainstructunit
        unit2=symeunit
      elseif(iswitch.eq.3)then
        open(teststructunit,file='teststruct.data',form='formatted',status='old')
        rewind(teststructunit)
        open(tymeunit,file='testinge.data',form='formatted',status='old')
        rewind(tymeunit)
        unit1=teststructunit
        unit2=tymeunit
      endif
 30   continue
      if((lshort.and.(nn_type_short.eq.1)).or.(lelec.and.(nn_type_elec.eq.1)))then
        read(unit2,*,END=40)num_atoms
        nstruct=nstruct+1
        num_atomselement(:)=0
        do i1=1,num_atoms
          read(unit2,*)idummy
        enddo
      elseif(nn_type_short.eq.2)then
        read(unit2,*,END=40)num_atoms,num_pairs
        nstruct=nstruct+1
        num_atomselement(:)=0
        do i1=1,num_pairs
          read(unit2,*)idummy
        enddo
      endif
      read(unit2,*)zdummy
      read(unit1,*)idummy,lperiodic
      allocate(zelem(num_atoms))
      if(lperiodic)then
        do i1=1,3
          read(unit1,*)zdummy
        enddo
      endif
      do i1=1,num_atoms
        read(unit1,*)zelem(i1)
      enddo
!!
!! analyze structure
      do i1=1,num_atoms
        num_atomselement(elementindex(zelem(i1)))&
          =num_atomselement(elementindex(zelem(i1)))+1
      enddo
      j1=1
      j2=1
      j3=1
      j4=1
      j5=1
      if(num_atomselement(1).gt.0)then
        j1=2
      endif
      if((nelem.gt.1).and.(num_atomselement(2).gt.0))then
        j2=2
      endif
      if((nelem.gt.2).and.(num_atomselement(3).gt.0))then
        j3=2
      endif
      if((nelem.gt.3).and.(num_atomselement(4).gt.0))then
        j4=2
      endif
      if((nelem.gt.4).and.(num_atomselement(5).gt.0))then
        j5=2
      endif
      if(nelem.eq.1)then
        if(.not.lperiodic)then
          temparray(num_atomselement(1),0,0,0,0,1)&
            =temparray(num_atomselement(1),0,0,0,0,1)+1
          temparray2(j1,j2,j3,j4,j5,1)=temparray2(j1,j2,j3,j4,j5,1)+1
        else
          temparray(num_atomselement(1),0,0,0,0,2)&
            =temparray(num_atomselement(1),0,0,0,0,2)+1
          temparray2(j1,j2,j3,j4,j5,2)=temparray2(j1,j2,j3,j4,j5,2)+1
        endif
      elseif(nelem.eq.2)then
        if(.not.lperiodic)then
          temparray(num_atomselement(1),num_atomselement(2),0,0,0,1)&
            =temparray(num_atomselement(1),num_atomselement(2),0,0,0,1)+1
          temparray2(j1,j2,j3,j4,j5,1)=temparray2(j1,j2,j3,j4,j5,1)+1
        else
          temparray(num_atomselement(1),num_atomselement(2),0,0,0,2)&
            =temparray(num_atomselement(1),num_atomselement(2),0,0,0,2)+1
          temparray2(j1,j2,j3,j4,j5,2)=temparray2(j1,j2,j3,j4,j5,2)+1
        endif
      elseif(nelem.eq.3)then
        if(.not.lperiodic)then
          temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),0,0,1)&
            =temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),0,0,1)+1
          temparray2(j1,j2,j3,j4,j5,1)=temparray2(j1,j2,j3,j4,j5,1)+1
        else
          temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),0,0,2)&
            =temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),0,0,2)+1
          temparray2(j1,j2,j3,j4,j5,2)=temparray2(j1,j2,j3,j4,j5,2)+1
        endif
      elseif(nelem.eq.4)then
        if(.not.lperiodic)then
          temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),num_atomselement(4),0,1)&
            =temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),num_atomselement(4),0,1)+1
          temparray2(j1,j2,j3,j4,j5,1)=temparray2(j1,j2,j3,j4,j5,1)+1
        else
          temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),num_atomselement(4),0,2)&
            =temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),num_atomselement(4),0,2)+1
          temparray2(j1,j2,j3,j4,j5,2)=temparray2(j1,j2,j3,j4,j5,2)+1
        endif
      elseif(nelem.eq.5)then
        if(.not.lperiodic)then
          temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),&
            num_atomselement(4),num_atomselement(5),1)&
            =temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),&
            num_atomselement(4),num_atomselement(5),1)+1
          temparray2(j1,j2,j3,j4,j5,1)=temparray2(j1,j2,j3,j4,j5,1)+1
        else
          temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),&
            num_atomselement(4),num_atomselement(5),2)&
            =temparray(num_atomselement(1),num_atomselement(2),num_atomselement(3),&
            num_atomselement(4),num_atomselement(5),2)+1
          temparray2(j1,j2,j3,j4,j5,2)=temparray2(j1,j2,j3,j4,j5,2)+1
        endif
      else
        write(ounit,*)'Error: analyzeinput works only for up to 5 elements'
        stop !'
      endif
!!
      deallocate(zelem)
      goto 30
 40   continue
      close(unit2)
      close(unit1)
!!
      write(ounit,*)'-------------------------------------------------------------'
      if(iswitch.eq.0)then
        write(ounit,*)'Training data set statistics:'
      elseif(iswitch.eq.1)then
        write(ounit,*)'Testing data set statistics:'
      elseif(iswitch.eq.2)then
        write(ounit,*)'Training data set statistics:'
      elseif(iswitch.eq.3)then
        write(ounit,*)'Testing data set statistics:'
      endif
      write(ounit,*)'Composition:              non-periodic:      periodic:      total:'
      write(ounit,'(2x,a2,3x,a2,3x,a2,3x,a2,3x,a2)')(trim(element(i1)),i1=1,nelem)
!!      if(nelem.eq.1)then
!!        do i1=1,max_num_atoms
!!          if((temparray(i1,1,1,1,1).ne.0).or.(temparray(i1,1,1,1,2).ne.0))then
!!            write(ounit,'(i4,15x,2i20)')i1,temparray(i1,1,1,1,1),temparray(i1,1,1,1,2)
!!          endif
!!        enddo
!!      endif
!!      if(nelem.eq.2)then
!!        do i1=1,max_num_atoms
!!          do i2=1,max_num_atoms
!!            if((temparray(i1,i2,1,1,1).ne.0).or.(temparray(i1,i2,1,1,2).ne.0))then
!!              write(ounit,'(i4,x,i4,10x,2i20)')i1,i2,temparray(i1,i2,1,1,1),temparray(i1,i2,1,1,2)
!!            endif
!!          enddo
!!        enddo
!!      endif
!!      if(nelem.eq.3)then
!!        do i1=1,max_num_atoms
!!          do i2=1,max_num_atoms
!!            do i3=1,max_num_atoms
!!              if((temparray(i1,i2,i3,1,1).ne.0).or.(temparray(i1,i2,i3,1,2).ne.0))then
!!                write(ounit,'(i4,x,i4,x,i4,6x,2i20)')i1,i2,i3,temparray(i1,i2,i3,1,1),temparray(i1,i2,i3,1,2)
!!              endif
!!            enddo
!!          enddo
!!        enddo
!!      endif
!!      if(nelem.eq.4)then
!!        do i1=1,max_num_atoms
!!          do i2=1,max_num_atoms
!!            do i3=1,max_num_atoms
!!              do i4=1,max_num_atoms
!!                if((temparray(i1,i2,i3,i4,1).ne.0).or.(temparray(i1,i2,i3,i4,2).ne.0))then
!!                  write(ounit,'(i4,x,i4,x,i4,x,i4,x,2i20)')i1,i2,i3,i4,temparray(i1,i2,i3,i4,1),temparray(i1,i2,i3,i4,2)
!!                endif
!!              enddo
!!            enddo
!!          enddo
!!        enddo
!!      endif
!!
!! general output      
        do i1=0,dim1
          do i2=0,dim2
            do i3=0,dim3
              do i4=0,dim4
                do i5=0,dim5
                if((i1+i2+i3+i4+i5).gt.0)then ! avoid output for structures with no elements
                  if((temparray(i1,i2,i3,i4,i5,1).ne.0).or.(temparray(i1,i2,i3,i4,i5,2).ne.0))then     ! write
                    if(nelem.eq.1)then
                      write(ounit,'(i4,16x,3i14)')i1,&
                        temparray(i1,i2,i3,i4,i5,1),temparray(i1,i2,i3,i4,i5,2),&
                        temparray(i1,i2,i3,i4,i5,1)+temparray(i1,i2,i3,i4,i5,2)
                    elseif(nelem.eq.2)then
                      write(ounit,'(i4,x,i4,11x,3i14)')i1,i2,&
                        temparray(i1,i2,i3,i4,i5,1),temparray(i1,i2,i3,i4,i5,2),&
                        temparray(i1,i2,i3,i4,i5,1)+temparray(i1,i2,i3,i4,i5,2)
                    elseif(nelem.eq.3)then
                      write(ounit,'(i4,x,i4,x,i4,6x,3i14)')i1,i2,i3,&
                        temparray(i1,i2,i3,i4,i5,1),temparray(i1,i2,i3,i4,i5,2),&
                        temparray(i1,i2,i3,i4,i5,1)+temparray(i1,i2,i3,i4,i5,2)
                    elseif(nelem.eq.4)then
                      write(ounit,'(i4,x,i4,x,i4,x,i4,x,3i14)')i1,i2,i3,i4,&
                        temparray(i1,i2,i3,i4,i5,1),temparray(i1,i2,i3,i4,i5,2),&
                        temparray(i1,i2,i3,i4,i5,1)+temparray(i1,i2,i3,i4,i5,2)
                    elseif(nelem.eq.5)then
                      write(ounit,'(i4,x,i4,x,i4,x,i4,x,i4,x,3i14)')i1,i2,i3,i4,i5,&
                        temparray(i1,i2,i3,i4,i5,1),temparray(i1,i2,i3,i4,i5,2),&
                        temparray(i1,i2,i3,i4,i5,1)+temparray(i1,i2,i3,i4,i5,2)
                    else
                      write(ounit,*)'ERROR in analyzeinput.f90, valid only up to 5 elements'
                      stop
                    endif
                  endif
! analyze also by system size
                  itemp=i1+i2+i3+i4+i5
                  nsystemsize(itemp,1)=nsystemsize(itemp,1)+temparray(i1,i2,i3,i4,i5,1)
                  nsystemsize(itemp,2)=nsystemsize(itemp,2)+temparray(i1,i2,i3,i4,i5,2)
                endif
                enddo
              enddo
            enddo
          enddo
        enddo
      write(ounit,*)'-------------------------------------------------------------'
!!
!! write summary
      if(iswitch.eq.0)then
        write(ounit,*)'Training set summary by element combinations: '
      elseif(iswitch.eq.1)then
        write(ounit,*)'Testing set summary by element combinations:'
      endif
      write(ounit,*)'   Subsystem:     Non-periodic:       Periodic:          Total:'
      do i1=1,2
       do i2=1,2
        do i3=1,2
         do i4=1,2
          do i5=1,2
            ctemp='              '
            if(i1.eq.2)then
              write(ctemp(1:2),'(a2)')element(1)
            endif
            if((i2.eq.2).and.(nelem.gt.1))then
              write(ctemp(4:5),'(a2)')element(2)
            endif
            if((i3.eq.2).and.(nelem.gt.2))then
              write(ctemp(7:8),'(a2)')element(3)
            endif
            if((i4.eq.2).and.(nelem.gt.3))then
              write(ctemp(10:11),'(a2)')element(4)
            endif
            if((i5.eq.2).and.(nelem.gt.4))then
              write(ctemp(13:14),'(a2)')element(5)
            endif
            if((temparray2(i1,i2,i3,i4,i5,1).gt.0).or.(temparray2(i1,i2,i3,i4,i5,2).gt.0))then
              write(ounit,'(x,a14,x,10x,i6,10x,i6,10x,i6)')ctemp,&
              temparray2(i1,i2,i3,i4,i5,1),temparray2(i1,i2,i3,i4,i5,2),&
              temparray2(i1,i2,i3,i4,i5,1)+temparray2(i1,i2,i3,i4,i5,2)
            endif
          enddo
         enddo
        enddo
       enddo
      enddo
      deallocate(temparray2)

!      do i1=0,nelem
!       do i2=0,nelem
!        do i3=0,nelem
!         do i4=0,nelem
!          do i5=0,nelem
!            if((temparray2(i1,i2,i3,i4,i5,1).gt.0).or.(temparray2(i1,i2,i3,i4,i5,2).gt.0))then
!              if(nelem.eq.1)then
!                write(ounit,'(a2,x,i6,x,i6)')element(i1),&
!                  temparray2(i1,i2,i3,i4,i5,1),temparray2(i1,i2,i3,i4,i5,2)
!              endif
!              if(nelem.eq.2)then
!                write(ounit,'(a2,x,a2,x,i6,x,i6)')element(i1),element(i2),&
!                  temparray2(i1,i2,i3,i4,i5,1),temparray2(i1,i2,i3,i4,i5,2)
!              endif
!              if(nelem.eq.3)then
!                write(ounit,'(a2,x,a2,x,a2,x,i6,x,i6)')element(i1),element(i2),element(i3),&
!                  temparray2(i1,i2,i3,i4,i5,1),temparray2(i1,i2,i3,i4,i5,2)
!              endif
!              if(nelem.eq.4)then
!                write(ounit,'(a2,x,a2,x,a2,x,a2,x,i6,x,i6)')&
!                  element(i1),element(i2),element(i3),element(i4),&
!                  temparray2(i1,i2,i3,i4,i5,1),temparray2(i1,i2,i3,i4,i5,2)
!              endif
!              if(nelem.eq.5)then
!                write(ounit,'(a2,x,a2,x,a2,x,a2,x,a2,x,i6,x,i6)')&
!                  element(i1),element(i2),element(i3),element(i4),element(i5),&
!                  temparray2(i1,i2,i3,i4,i5,1),temparray2(i1,i2,i3,i4,i5,2)
!              endif
!            endif
!          enddo
!         enddo
!        enddo
!       enddo
!      enddo
!      deallocate(temparray2)
!!
      deallocate(temparray)

! write system size summary
      write(ounit,*)'-------------------------------------------------------------'
      if(iswitch.eq.0)then
        write(ounit,*)'Training set summary by system size (number of atoms): '
      elseif(iswitch.eq.1)then
        write(ounit,*)'Testing set summary by system size (number of atoms):'
      endif
      write(ounit,*)' System Size:     Non-periodic:       Periodic:          Total:'
      do i1=1,max_num_atoms
        if((nsystemsize(i1,1).gt.0).or.(nsystemsize(i1,2).gt.0))then
          write(ounit,'(i14,i18,2i16)')i1,nsystemsize(i1,1),nsystemsize(i1,2),&
            nsystemsize(i1,1)+nsystemsize(i1,2)
        endif
      enddo


      return
      end

