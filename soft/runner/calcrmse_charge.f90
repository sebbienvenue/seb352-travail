!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - ewaldenergies_para.f90  
!!
      subroutine calcrmse_charge(ndim,npoints,ncharges,&
           zelem_mpi,num_atoms_mpi,rmse_charge,mad_charge,&
           totalcharge_mpi,rmse_totalcharge,mad_totalcharge,&
           atomcharge_mpi,nnatomcharge_mpi,nnchargesum_mpi)
!!
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer npoints                                 ! in
      integer ndim                                    ! in
      integer num_atoms_mpi(ndim)                     ! in
      integer ncharges                                ! in/out
      integer i1,i2                                   ! internal
      integer zelem_mpi(ndim,max_num_atoms)           ! in
!!
      real*8 atomcharge_mpi(ndim,max_num_atoms)       ! in
      real*8 nnatomcharge_mpi(ndim,max_num_atoms)     ! in
      real*8 rmse_charge                              ! in/out
      real*8 rmse_totalcharge                         ! in/out
      real*8 mad_charge                               ! in/out
      real*8 mad_totalcharge                          ! in/out
      real*8 totalcharge_mpi(ndim)                    ! in
      real*8 chargesum                                ! internal
      real*8 nnchargesum_mpi(ndim)                    ! out 
!!
!!
      nnchargesum_mpi(:)=0.0d0
!!
      do i1=1,npoints
        chargesum=0.0d0
        do i2=1,num_atoms_mpi(i1)
          if(.not.lupdatebyelement)then
            ncharges=ncharges+1
            rmse_charge=rmse_charge +(atomcharge_mpi(i1,i2)-nnatomcharge_mpi(i1,i2))**2
            mad_charge =mad_charge  +abs(atomcharge_mpi(i1,i2)-nnatomcharge_mpi(i1,i2))
          else
            if(elemupdate.eq.zelem_mpi(i1,i2))then
              ncharges=ncharges+1
              rmse_charge=rmse_charge +(atomcharge_mpi(i1,i2)-nnatomcharge_mpi(i1,i2))**2
              mad_charge =mad_charge +abs(atomcharge_mpi(i1,i2)-nnatomcharge_mpi(i1,i2))
            endif
          endif ! lupdatebyelement
          chargesum=chargesum+nnatomcharge_mpi(i1,i2)
        enddo ! i2
        nnchargesum_mpi(i1)=chargesum
        rmse_totalcharge=rmse_totalcharge +(totalcharge_mpi(i1)-chargesum)**2
        mad_totalcharge =mad_totalcharge +abs(totalcharge_mpi(i1)-chargesum)
      enddo ! i1
!!
      return
      end
