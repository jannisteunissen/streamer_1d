! Requires Fortran 2008 features: using internal procedure as argument
program test_m_transport_schemes
  use m_transport_schemes
  use m_time_steppers

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  real(dp), parameter   :: test_velocity   = -1.0_dp
  real(dp), parameter   :: test_diff_coeff = 0.0_dp
  real(dp), parameter   :: test_length     = 1.0_dp
  real(dp), parameter   :: end_time        = 0.1_dp * test_length / abs(test_velocity)
  integer, parameter    :: max_n_cells     = 5120

  real(dp)              :: delta_x, delta_t, time
  integer               :: n, n_cells, i_scheme
  integer               :: time_cntr, max_time_cntr
  real(dp)              :: diff, max_diff, sum_diff
  real(dp), allocatable :: dens_c(:)

  procedure(STEP_1d_type), pointer  :: time_stepper
  procedure(TS_dd_1d_type), pointer :: transport_scheme

  time_stepper => STEP_rk4_1d

  print *, "Testing m_transport_schemes.f90 implementation"

  do i_scheme = 1, 4

     print *, ""

     select case (i_scheme)
     case (1)
        transport_scheme => TS_dd_1d_up
        print *, " Scheme: first order upwind"
     case (2)
        transport_scheme => TS_dd_1d
        print *, " Scheme: second order centered differences"
     case (3)
        transport_scheme => TS_dd_1d_up_fl
        print *, " Scheme: ~third order upwind with Koren limiter"
     case (4)
        transport_scheme => TS_dd_1d_up_fl
        print *, " Scheme: Willem ~third order upwind with Koren limiter"
     end select

     print *, "    N_cells        max abs error            2-norm error"

     n_cells = 20

     do
        n_cells = n_cells * 2
        if (n_cells > max_n_cells) exit
        allocate(dens_c(n_cells))

        delta_x = test_length / (n_cells - 1)
        delta_t        = 0.9_dp * delta_x / abs(test_velocity) ! Cfl number 0.9
        max_time_cntr  = nint(end_time / delta_t)

        do n = 1, n_cells
           dens_c(n) = analytical_sol((n-1) * delta_x, 0.0_dp)
        end do

        do time_cntr = 1, max_time_cntr
           call time_stepper(dens_c, time, delta_t, set_derivs)
        end do

        ! Check the error
        sum_diff = 0.0_dp
        max_diff = 0.0_dp
        do n = 1, n_cells
           diff = dens_c(n) - analytical_sol((n-1) * delta_x, max_time_cntr * delta_t)
           sum_diff = sum_diff + diff**2
           max_diff = max(max_diff, abs(diff))
        end do
        print *, n_cells, max_diff, sqrt(sum_diff/n_cells)
        deallocate(dens_c)
     end do
  end do

contains

  subroutine set_derivs(dens_c, time, derivs)
    real(dp), intent(in) :: dens_c(:), time
    real(dp), intent(out) :: derivs(:)
    real(dp) :: vel_f(size(dens_c)-1)
    real(dp) :: dif_f(size(dens_c)-1)
    real(dp) :: flux_f(size(dens_c)-1)

    vel_f = test_velocity
    dif_f = test_diff_coeff / delta_x
    call transport_scheme(dens_c, vel_f, dif_f, flux_f)

    flux_f = flux_f / delta_x
    derivs = 0.0_dp
    derivs(1:size(derivs)-1) = -flux_f
    derivs(2:size(derivs)) = derivs(2:size(derivs)) + flux_f
  end subroutine set_derivs

  real(dp) function analytical_sol(x, t)
    real(dp), intent(in) :: x, t
    real(dp)             :: adv_x
    real(dp), parameter  :: wave_width = 0.15_dp, grad_width = 0.05_dp

    adv_x = (x - test_velocity * t) / test_length

    if (adv_x < 0.5_dp - wave_width) then
       analytical_sol = 0.0_dp
    else if (adv_x < 0.5_dp - wave_width + grad_width) then
       analytical_sol = smootherstep(0.5_dp - wave_width, 0.5_dp - wave_width + grad_width, adv_x)
    else if (adv_x < 0.5_dp + wave_width) then
       analytical_sol = 1.0_dp
    else if (adv_x < 0.5_dp + wave_width + grad_width) then
       analytical_sol = smootherstep(0.5_dp + wave_width + grad_width, 0.5_dp + wave_width, adv_x)
    else
       analytical_sol = 0.0_dp
    end if

  end function analytical_sol

  real(dp) function smootherstep(edge0, edge1, xx)
    real(dp) :: edge0, edge1, xx

    xx           = (xx - edge0)/(edge1 - edge0)
    smootherstep =  xx**3 * (xx * (xx * 6 - 15) + 10)
  end function smootherstep


end program test_m_transport_schemes
