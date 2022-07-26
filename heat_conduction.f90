!========================Program Description========================
! 1-D heat conduction solver using FDM on Fourier Law
!
! Function to solve: q = k*dT/dx
!
!===================================================================
module heat_transfer
  implicit none
  integer, public, parameter :: p = selected_real_kind(15,100)
  integer :: i
  
  contains
!========================SUBROUTINE fourier_oned====================
! calculates temperature field on 1-D rod due to conduction
!===================================================================
  subroutine fourier_oned(T, q, k, nx, dx, T_west)
    real(kind=p), dimension(nx), intent(inout) :: T
    real(kind=p), intent(in) :: q, T_west, k, dx
    integer, intent(in) :: nx
    real(kind=p), dimension(:), allocatable :: T_last !temporary array to hold previos iteration values
    real(kind=p) :: res = 0.5 !tolerance
    integer :: count = 0

    allocate(T_last(nx))
    T_last(:) = T_west

    residual_loop: do
      steady_state_loop: do i = 1, size(T)
        if (i .eq. size(T)) then !last node
          T_last(i) = 0.5*(q*k/dx + T(i-1))
        else if (i .eq. 1) then !first node
          T_last(i) = 0.5*(T(i+1) + T_west)
        else !interior nodes
          T_last(i) = 0.5*(T(i+1) + T(i-1))
        end if
      end do steady_state_loop
      count = count + 1
      print *, "iteration = ", count, ",  residual = ", maxval(abs(T_last-T))

      if (maxval(abs(T_last - T)) .le. res) then
        T = T_last
        exit
      end if

      T = T_last
    end do residual_loop
  end subroutine fourier_oned
end module heat_transfer

program heat_conduction
  use heat_transfer
  implicit none
  integer :: nx
  real(kind=p) :: q, k, T_west, dx
  real(kind=p), dimension(:), allocatable :: T
  
  nx = 5
  allocate(T(nx))
  T_west = 300.0 !room temp in Kelvin
  T(:) = T_west 
  k = 1.2_p !W/m^2-K
  dx = 0.5_p !distance delta
  q = 500.0_p

  call fourier_oned(T, q, k, nx, dx, T_west)
  
  print *, T
  open(1, file="temperature_profile.dat")
  write(1,*) T
  close(1)

end program heat_conduction