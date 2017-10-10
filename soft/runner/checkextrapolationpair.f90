!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - predictionpair.f90
!!
!! This subroutine checks the symmetry functions for extrapolation
!!
      subroutine checkextrapolationpair(num_pairs,pindex,npairs,pairindex,&
             max_num_pairs,maxnum_funcvaluesp,num_funcvaluesp,&
             pairs_charge,symfunctionp,minvalue_short_pair,maxvalue_short_pair,slabel)
!!
      use fileunits
!!   
      implicit none
!!
      integer pairindex(102,102)
      integer max_num_pairs                    ! in
      integer maxnum_funcvaluesp            ! in
      integer num_funcvaluesp(npairs)        ! in
      integer npairs                        ! in
      integer i1,i3                     ! internal
      integer num_pairs                       ! in
      integer pindex(max_num_pairs)               ! in
!!      integer pindex(num_pairs_para)               ! in
      integer pairs_charge(2,max_num_pairs)
!!
      real*8 symfunctionp(maxnum_funcvaluesp,max_num_pairs) ! in
!!      real*8 symfunctionp(maxnum_funcvaluesp,num_pairs_para) ! in
      real*8 minvalue_short_pair(npairs,maxnum_funcvaluesp)            ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvaluesp)            ! in
      real*8 threshold                                 ! internal
!!
      character*5 slabel                    ! in
!!
!!      
      threshold=1.0d-8
!!
      do i1=1,num_pairs
        do i3=1,num_funcvaluesp(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))))
          if((symfunctionp(i3,i1)-maxvalue_short_pair(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))),i3))&
            .gt.threshold)then
            write(ounit,'(a,a5,x,2i5,x,a,2f18.8)')'### EXTRAPOLATION WARNING ### ',slabel,&
            pindex(i1),i3,'too large ',symfunctionp(i3,i1),&
            maxvalue_short_pair(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))),i3)
          elseif((-symfunctionp(i3,i1)+minvalue_short_pair(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))),i3))&
            .gt.threshold)then
            write(ounit,'(a,a5,x,2i5,x,a,2f18.8)')'### EXTRAPOLATION WARNING ### ',slabel,&
            pindex(i1),i3,'too small ',symfunctionp(i3,i1),&
            minvalue_short_pair(pairindex(pairs_charge(1,pindex(i1)),pairs_charge(2,pindex(i1))),i3)
          endif
        enddo ! i3
      enddo ! i1
!!
      return
      end
