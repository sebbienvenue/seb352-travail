!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! Purpose: get a part of the squared force errors for the RMSE calculation of the forces
!!          Each force component is considered independently.

!! called by:
!! - geterror.f90  
!! - geterrorpair.f90  
!!
      subroutine calcrmse_forces(npoints,&
        nforces,rmse_forces,mad_forces,&
        force_list_local,nnforce_list_local,maxforce_local)
!!
      use fileunits
      use fittingoptions
      use globaloptions
      use structures
!!
      implicit none
!!
      integer npoints                                 ! in
      integer nforces                                 ! in/out
!!
      integer i1,i2,i3
!!
      real*8 force_list_local(3,max_num_atoms,nblock)
      real*8 nnforce_list_local(3,max_num_atoms,nblock)
      real*8 rmse_forces                              ! in/out
      real*8 mad_forces                               ! in/out
      real*8 maxforce_local                           ! in
!!
!!

      do i1=1,npoints                   ! all structures
        do i2=1,num_atoms_list(i1)      ! all atoms of this structure
          do i3=1,3                     ! fx,fy,fz
!! don't include force in the RMSE calculation if it is not used for updating:
!! lupdatebyelement does make sense only for nn_type_short = 1
            if((nn_type_short.eq.1).and.lupdatebyelement.and.(zelem_list(i1,i2).ne.elemupdate)) goto 99
!! don't include force component in RMSE if reference force is larger than maxforce
!!            write(ounit,*)force_list_local(i3,i2,i1),maxforce_local
            if(abs(force_list_local(i3,i2,i1)).gt.maxforce_local) goto 99
!!              write(ounit,*)'using force'
              nforces=nforces+1
              rmse_forces=rmse_forces &
                +(force_list_local(i3,i2,i1)-nnforce_list_local(i3,i2,i1))**2

              mad_forces =mad_forces &
                +abs(force_list_local(i3,i2,i1)-nnforce_list_local(i3,i2,i1))
 99         continue
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
      return
      end
