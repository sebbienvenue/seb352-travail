!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90 
!! - geterrorpair.f90
!! - preconditionpair.f90
!!
      subroutine readfunctionspair(unit,npoints,npairs,&
         max_num,num_atoms_list,size_pairs_list,&
         pairindex,maxnum_funcvaluesp,num_funcvaluesp,&
         zelemp_list,symfunction_list,totalcharge_list,&
         totalenergy_list,shortenergy_list,&
         ewaldenergy_list,ldebug)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!     
      integer unit
      integer i1,i2,i3                                                    ! internal
      integer max_num                                                     ! in
      integer maxnum_funcvaluesp                                          ! in
      integer npairs                                                      ! in
      integer npoints                                                     ! in
      integer num_atoms_list(nblock)                                      ! out
      integer num_funcvaluesp(npairs)                                     ! out
      integer pairindex(102,102)                                          ! in
      integer size_pairs_list(nblock)                                     ! out
      integer zelemp_list(2,nblock,max_num)                               ! out
!!
      real*8  ewaldenergy_list(nblock)                                    ! out
      real*8  shortenergy_list(nblock)                                    ! out
      real*8  totalcharge_list(nblock)                                    ! out
      real*8  totalenergy_list(nblock)                                    ! out 
      real*8  symfunction_list(maxnum_funcvaluesp,max_num,nblock)        ! out
!!
      logical ldebug                                                      ! in
!!
     do i1=1,npoints
        read(unit,*)num_atoms_list(i1),size_pairs_list(i1)
        do i2=1,size_pairs_list(i1) 
          read(unit,*)zelemp_list(1,i1,i2),zelemp_list(2,i1,i2),&
                       (symfunction_list(i3,i2,i1),i3=1,&
                       num_funcvaluesp(pairindex(zelemp_list(1,i1,i2),&
                                                 zelemp_list(2,i1,i2))))
        enddo ! i2
        read(unit,*)totalcharge_list(i1),totalenergy_list(i1),&
                     shortenergy_list(i1),ewaldenergy_list(i1)
     enddo ! i1

!!
      return
      end
