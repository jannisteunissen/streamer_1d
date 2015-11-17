program test_m_time_steppers
  use m_time_steppers

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n_vars = 10, max_num_steps = 1000, max_adap_steps = 30
  real(dp), parameter :: end_time = 1.0_dp, max_dt = huge(1.0_dp)
  integer :: n, n_steps, i_stepper
  real(dp) :: var_list(n_vars), sol_start(n_vars), sol_end(n_vars)
  real(dp) :: max_error, min_error, max_errs(n_vars)
  real(dp) :: time, dt, new_dt
  procedure(STEP_1d_type), pointer :: time_stepper

  do n = 1, n_vars
     sol_start(n) = get_solution(n, 0.0_dp)
     sol_end(n) = get_solution(n, end_time)
  end do

  print *, "Testing implementation of m_time_steppers.f90"
  print *, ""

  do i_stepper = 1, 3
     select case (i_stepper)
     case (1)
        time_stepper => STEP_forward_euler_1d
        print *, "  Forward Euler"
     case (2)
        time_stepper => STEP_rk2_1d
        print *, "  Runge-Kutta second order"
     case (3)
        time_stepper => STEP_rk4_1d
        print *, "  Runge-Kutta fourth order"
     end select

     print *, "  N_steps  Maximum difference with analytical solution"

     n_steps = 1

     do
        n_steps  = n_steps * 2
        if (n_steps > max_num_steps) exit
        dt       = end_time / n_steps
        time     = 0.0_dp
        var_list = sol_start

        do n = 1, n_steps
           call time_stepper(var_list, time, dt, set_derivatives)
        end do

        write(*, "(I10, E16.8)") n_steps, maxval(abs(var_list - sol_end))
     end do

     print *, ""
  end do

  ! Separately test adaptive stepsize algorithm
  print *, "  Runge-Kutta fourth order adaptive stepsize"
  print *, "  N_steps  Maximum difference with analytical solution"
  max_error = 1.0_dp
  min_error = epsilon(1.0_dp)

  do n = 1, max_adap_steps
     max_error = max_error * 0.5_dp
     time      = 0.0_dp
     var_list  = sol_start
     n_steps   = 0
     dt        = 1.0e-2_dp * max_error**0.25_dp

     do
        max_errs = max_error * maxval(abs(var_list)) + min_error
        call STEP_rk4a_1d(var_list, max_errs, time, dt, max_dt, new_dt, set_derivatives)
        dt = min(new_dt, end_time - time)
        n_steps = n_steps + 1
        if (time > end_time - epsilon(1.0_dp)) exit
     end do

     write(*, "(I10, E16.8)") n_steps, maxval(abs(var_list - sol_end))
  end do
  print *, ""

contains

  real(dp) function get_solution(eq_ix, time)
    integer, intent(in)  :: eq_ix
    real(dp), intent(in) :: time

    get_solution = sin(eq_ix * time)
  end function get_solution

  subroutine set_derivatives(vars, time, derivs)
    real(dp), intent(in)  :: vars(:), time
    real(dp), intent(out) :: derivs(:)
    integer :: eq_ix

    do eq_ix = 1, size(vars)
       derivs(eq_ix) = eq_ix * cos(eq_ix * time)
    end do

  end subroutine set_derivatives

end program test_m_time_steppers
