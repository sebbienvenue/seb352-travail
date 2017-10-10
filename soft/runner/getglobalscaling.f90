!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!> MG: Subroutine for the computation of the global scaling factor A
!! MG: required in the element-decoupled Kalman filter

!! called by:
!! - optimize_short_combined.f90
!!
      subroutine getglobalscaling(num_weights_short_atomic_free,&
                deshortdw,corrmatrix_list,maxcorrdim,corrdim,&
                num_atoms_element,num_atoms,global_scaling)
!!
      use fileunits
      use fittingoptions !! lupdatebyelement, nucelem, elemupdate
      use globaloptions  !! maxnum_weights_short_atomic
!!
      implicit none

      integer i2                                                  ! internal
      integer i3                                                  ! internal
      integer maxcorrdim                                          ! in
      integer num_atoms                                           ! in
      integer num_atoms_element(nelem)                            ! in
      integer num_weights_short_atomic_free(nelem)                ! in
      integer corrdim(nelem)                                      ! in
!!
      real*8 global_scaling                                       ! in/out
      real*8 corrmatrix_list(maxcorrdim,nelem)                    ! in
      real*8 deshortdw(maxnum_weights_short_atomic,1,nelem)       ! in
      real*8 ddot                                                 ! internal
      real*8, dimension(:), allocatable :: deshortdw_temp         ! internal
      real*8, dimension(:), allocatable :: coh                    ! internal

!! MG: Initialize
      global_scaling = 0d0

!! MG: 1) Compute contribution to global scaling factor for every element
!! MG:       all_contributions = SUM^nelem[ J^T(elem)*P(elem)*J(elem) ]
      do i2=1,nelem
!!
!! MG: Check if update for element is requested [*]
        if(lupdatebyelement.and.(nucelem(i2).ne.elemupdate))then
          continue
!!
        else
!! MG: Do update only if elements are present [*]
          if(num_atoms_element(i2).gt.0)then
            allocate(deshortdw_temp(num_weights_short_atomic_free(i2)))
!! MG: reduce deshortdw to array containing only the free weights [*]
            do i3=1,num_weights_short_atomic_free(i2)
              deshortdw_temp(i3) = deshortdw(i3,1,i2)
            enddo ! i3
!! MG: Fix2: Rescale derivatives by number of atoms
            deshortdw_temp = deshortdw_temp*dble(num_atoms)
!!
!! MG: Use DSPMV subroutine to compute coh=P_elem(n-1)*J_elem(n)
            allocate(coh(num_weights_short_atomic_free(i2)))
            coh(:) = 0d0
            call dspmv('l',num_weights_short_atomic_free(i2),1.d0,&
              corrmatrix_list(1,i2),deshortdw_temp,1,0.d0,coh,1)
!!
!! MG: Use DDOT function to compute J^T*P*J for element
            global_scaling = global_scaling + ddot(num_weights_short_atomic_free(i2),&
              deshortdw_temp,1,coh,1)
!!
            deallocate(deshortdw_temp,coh)
!! MG: Fix2: Extended reach of if-statement.
          endif ! num_atoms_element(i2).gt.0
        endif ! lupdatebyelement.and.(nucelem(i2).ne.elemupdate)
      enddo ! i2
!! MG: 2) Compute scaling factor as A=[lambda + all_contributions]^-1
!! MG: Fixed: Adressed the kalmanlambda array incorrectly, resulting in a
!! MG: Fixed: lambda of 0. 
      global_scaling = 1.d0/(kalmanlambda(1)+global_scaling)
!!
      return
      end
