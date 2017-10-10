!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: get array num_atoms_all, contains number of atoms in each training point

!! called by:
!! - fitting.f90
!! - fittingpair.f90 
!!
      subroutine getnumatomsall(ntrain,num_atoms_all)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!     
      integer ntrain                                                ! in
      integer num_atoms_all(ntrain)                                 ! out
      integer i1,i2                                                 ! internal
      integer ndummy                                                ! internal
      integer idummy                                                ! internal
!!
      real*8 rdummy                                                 ! internal 
      character*200 :: filenametemp                                 ! internal
!!
!!
      num_atoms_all(:)=0
!!
      if(lshort)then
        open(symunit,file='function.data',form='formatted',status='old')
        rewind(symunit)
        do i1=1,ntrain
          if(nn_type_short.eq.1)then
            read(symunit,*)num_atoms_all(i1)
            do i2=1,num_atoms_all(i1)
              read(symunit,*)idummy
            enddo ! i2
          elseif(nn_type_short.eq.2)then
            read(symunit,*)num_atoms_all(i1),ndummy
            do i2=1,ndummy
              read(symunit,*)idummy
            enddo ! i2
          endif
          read(symunit,*)rdummy 
        enddo ! i1
        close(symunit)
      elseif(lelec.and.(nn_type_elec.eq.1))then  ! we also need num_atoms_all if lshort is switched off
        open(symeunit,file='functione.data',form='formatted',status='old')
        rewind(symeunit) !'
        do i1=1,ntrain
          read(symeunit,*)num_atoms_all(i1)
          do i2=1,num_atoms_all(i1)
            read(symeunit,*)idummy
          enddo ! i2
          read(symeunit,*)rdummy 
        enddo ! i1
      elseif(lnntb)then
        if(nntb_flag(3))then
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
          'function_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          open(symhextoffunit,file=filenametemp,form='formatted',status='old')
          rewind(symhextoffunit)
          do i1=1,ntrain
            num_atoms_all(i1) = 3
          enddo
        endif
      else
      
        write(ounit,*)'ERROR in getnumatomsall'
        stop
      endif
!!
      return
      end
