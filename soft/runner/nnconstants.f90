      module nnconstants 
      implicit none

      real*8 pi 
      real*8 rad2deg 

      contains
      subroutine get_nnconstants()
      use mpi_mod
      implicit none
!!
      pi        = 4.d0*datan(1.d0)
      rad2deg   = 180.d0/pi

      end subroutine get_nnconstants 
      end module nnconstants 
