!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! This is a multipurpose subroutine to scale general symfunctions

!! called by:
!! gethextoffoutput.f90
!! CMH FIXME mode 3 calls required here
!! In mode 3 we will be passing one triplet at a time to this process
!! for all points, therefore rather than passing full arrays of info
!! e.g. num_funcvalues_hextoff(ntriplets), we pass for a particular value of ntriplets
!! so above in mode 3 we will require temp containers to pass from
!! and pass back to for error information

      subroutine scalesym_hextoff(npoints,&
         maxnum_funcvalues_local,num_funcvalues_local,&
         symfunction_local,&
         minvalue_local,maxvalue_local,avvalue_local,&
         scmin_local,scmax_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer npoints                                                        ! in
      integer maxnum_funcvalues_local                                        ! in
      integer num_funcvalues_local                                           ! in
      integer i1,i3                                                          ! internal
!!
      real*8 symfunction_local(maxnum_funcvalues_local,nblock)               ! in/out
      real*8 minvalue_local(maxnum_funcvalues_local)                         ! in
      real*8 maxvalue_local(maxnum_funcvalues_local)                         ! in
      real*8 avvalue_local(maxnum_funcvalues_local)                          ! in
      real*8 scmin_local                                                     ! in
      real*8 scmax_local                                                     ! in
!!
      do i1=1,npoints
        do i3=1,num_funcvalues_local
          if(lcentersym.and..not.lscalesym)then
!! For each symmetry function remove the CMS of the respective element 
            symfunction_local(i3,i1)=symfunction_local(i3,i1) &
            -avvalue_local(i3)
          elseif(lscalesym.and..not.lcentersym)then
!! Scale each symmetry function value for the respective element
            symfunction_local(i3,i1)=&
           (symfunction_local(i3,i1)-minvalue_local(i3))/ &
           (maxvalue_local(i3)-minvalue_local(i3))&
            *(scmax_local-scmin_local) + scmin_local
          elseif(lscalesym.and.lcentersym)then
            symfunction_local(i3,i1)=&
           (symfunction_local(i3,i1)-avvalue_local(i3))&
           / (maxvalue_local(i3)-minvalue_local(i3))
          else
          endif
        enddo ! i3
      enddo ! i1
!!
      return
      end
