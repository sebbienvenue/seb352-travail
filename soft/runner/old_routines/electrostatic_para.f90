!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!!            - prediction.f90
!!            - predictionpair.f90
!!
      subroutine electrostatic_para(num_atoms,atomcharge,xyzstruct,elecenergy)
!!
      use mpi_mod
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer num_atoms
      integer i1,i2
      integer i1_start,i2_start
      integer i1_end,i2_end
      integer icount,jcount
      integer kcount
      integer num,totnum                              ! internal
      integer n_start,n_end                           ! internal
!!
      real*8 xyzstruct(3,max_num_atoms)               ! in
      real*8 atomcharge(max_num_atoms)                ! in
      real*8 elecenergy                               ! out
      real*8 distance                                 ! internal
      real*8 fscreen                                  ! internal
      real*8 fscreenderiv                             ! internal
!!
      logical lstart
      logical lend
!!
!! For this subroutine the use of Hartree energy and Bohr length units is mandatory!
!! 1 Hartree = \frac{\hbar^2}{m_e a_0^2} = \frac{e^2}{4\pi \epsilon_0 a_0}
!!
      elecenergy     =0.0d0
!!
      if(num_atoms.gt.1) then
!!
!! determine the total number of pair interactions to be calculated
        totnum=num_atoms*(num_atoms+1)/2-num_atoms
! THE SAME:        totnum=num_atoms*(num_atoms-1)/2
!      write(ounit,*)'number of interactions ',totnum
!! determine the number of loops for this process
        call mpifitdistribution(totnum,num,n_start,n_end)
!!
!! i1 specifies column
!! i2 specifies row
        icount=num_atoms-1 ! initially number of elements in longest column
        jcount=0
        kcount=0 ! column counter
        lstart=.false.
        lend=.false.
!!        write(ounit,*)'n_start,n_end ',n_start,n_end
 9      continue
        kcount=kcount+1

!!
!! CHANGE ANDI: GFORTRAN: gfortran forbids .eq. for logicals, use .eqv. instead.
!!                        ifort allows both.
!!
       !if(lstart.eq..false.)then
        if(lstart.eqv..false.)then
!! END CHANGE

          if(n_start.le.(jcount+icount))then
            i1_start=kcount !! determine starting column
            i2_start=n_start-jcount+kcount !! determine row in column          
!! => this process starts at point (i2_start,i1_start)=(row,column)
            lstart=.true.
          endif
        endif ! lstart

!!
!! CHANGE ANDI: GFORTRAN: gfortran forbids .eq. for logicals, use .eqv. instead.
!!                        ifort allows both.
!!
       !if(lend.eq..false.)then
        if(lend.eqv..false.)then
!! END CHANGE

          if(n_end.le.(jcount+icount))then
            i1_end=kcount !! determine starting column
            i2_end=n_end-jcount+kcount !! determine row in column          
!! => this process ends at point (i2_end,i1_end)=(row,column)
            lend=.true.
          endif
        endif ! lend
        if(icount.gt.1) then
          jcount=jcount+icount ! sum of elements in previous columns
          icount=icount-1 ! number of elements in next column
          goto 9
        endif
!! checks
        if(.not.lstart)then
          write(ounit,*)'ERROR: electrostatic_para: start not found'
          stop
        endif
        if(.not.lend)then
          write(ounit,*)'ERROR: electrostatic_para: end not found'
          stop
        endif
!!
!!        write(ounit,'(a,i4,a,i6,a1,i6,a4,i6,a1,i6,a1)')&
!!          ' Process ',mpirank,' calculates (',i2_start,',',i1_start,') to (',&
!!          i2_end,',',i1_end,')'
!!'
!! loop over all interaction pairs of this process
        do i1=i1_start,i1_end  ! column
          do i2=i2_start,i2_end ! row
!! decide if we are inside the lower triangle => only calculate interaction in this case
            if(i2.gt.i1)then
!!              write(ounit,*)'calculating pair ',(i2,i1)
              distance=(xyzstruct(1,i1)-xyzstruct(1,i2))**2 + &
                       (xyzstruct(2,i1)-xyzstruct(2,i2))**2 + &
                       (xyzstruct(3,i1)-xyzstruct(3,i2))**2
              distance=dsqrt(distance)
!!              write(ounit,*)'TEST distance ',distance
!!
              if(lscreen) then
                call getscreenfunctionforelectrostatics(&
                       distance,fscreen,fscreenderiv,0.0d0)
!!
                elecenergy=elecenergy + atomcharge(i1)*atomcharge(i2)/ &
                              distance*fscreen
              else
                elecenergy=elecenergy + atomcharge(i1)*atomcharge(i2)/distance
              endif
            endif
          enddo ! i2
        enddo ! i1
!!
!!
!!
      endif ! num_atoms.gt.1
!!
 20   continue
      return
!!
      end
