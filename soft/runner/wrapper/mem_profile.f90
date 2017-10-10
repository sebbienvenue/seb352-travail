module mem_profile
  use fileunits
  use mpi_mod

  implicit none

  ! explicitely make stuff public; so everything is private by default
  private
  public :: array_size

  ! generic format for outputs
  character(len=*), parameter :: OUT_FMT = '("!! MEM_PROFILE !! Size of ",A," = ",F8.2," MiB")'
  ! size of integer/float/whatever
  real(kind=8), parameter :: INT_SIZE = 4.0d0
  real(kind=8), parameter :: DBL_SIZE = 8.0d0

  ! interface for generic routines
  interface array_size
    module procedure array_size_dbl_1
    module procedure array_size_dbl_2
    module procedure array_size_dbl_3
    module procedure array_size_dbl_4
    module procedure array_size_int_1
    module procedure array_size_int_2
    module procedure array_size_int_3
  end interface array_size

  contains
  subroutine print_size(array_name, array_size)
    implicit none

    character(len=*), intent(in) :: array_name
    real(kind=8), intent(in) :: array_size

    if (mpirank == 0) then
      write(ounit, OUT_FMT) array_name, array_size
    end if
  end subroutine

  subroutine array_size_dbl_1(array_name, array)
    implicit none

    character(len=*), intent(in) :: array_name
    real(kind=8), intent(in) :: array(:)

    call print_size(array_name, real(size(array), kind=8) * DBL_SIZE / 1048576.0d0)
  end subroutine

  subroutine array_size_dbl_2(array_name, array)
    implicit none

    character(len=*), intent(in) :: array_name
    real(kind=8), intent(in) :: array(:,:)

    call print_size(array_name, real(size(array), kind=8) * DBL_SIZE / 1048576.0D0)
  end subroutine

  subroutine array_size_dbl_3(array_name, array)
    implicit none

    character(len=*), intent(in) :: array_name
    real(kind=8), intent(in) :: array(:,:,:)

    call print_size(array_name, real(size(array), kind=8) * DBL_SIZE / 1048576.0D0)
  end subroutine

  subroutine array_size_dbl_4(array_name, array)
    implicit none

    character(len=*), intent(in) :: array_name
    real(kind=8), intent(in) :: array(:,:,:,:)

    call print_size(array_name, real(size(array), kind=8) * DBL_SIZE / 1048576.0D0)
  end subroutine

  subroutine array_size_int_1(array_name, array)
    implicit none

    character(len=*), intent(in) :: array_name
    integer, intent(in) :: array(:)

    call print_size(array_name, real(size(array), kind=8) * INT_SIZE / 1048576.0D0)
  end subroutine

  subroutine array_size_int_2(array_name, array)
    implicit none

    character(len=*), intent(in) :: array_name
    integer, intent(in) :: array(:,:)

    call print_size(array_name, real(size(array), kind=8) * INT_SIZE / 1048576.0D0)
  end subroutine

  subroutine array_size_int_3(array_name, array)
    implicit none

    character(len=*), intent(in) :: array_name
    integer, intent(in) :: array(:,:,:)

    call print_size(array_name, real(size(array), kind=8) * INT_SIZE / 1048576.0D0)
  end subroutine
end module mem_profile
