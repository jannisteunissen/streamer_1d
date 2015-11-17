module m_transport_schemes

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  interface
     subroutine TS_dd_1d_type(dens_c, vel_f, dif_f, flux_f)
       import
       real(dp), intent(in)    :: dens_c(:), vel_f(:), dif_f(:)
       real(dp), intent(inout) :: flux_f(:)
     end subroutine TS_dd_1d_type
  end interface

  public :: TS_dd_1d_type
  public :: TS_dd_1d
  public :: TS_dd_1d_up
  public :: TS_dd_1d_up_fl
  public :: TS_dd_1d_up_fl_test

contains

  !> Second order centered differences scheme. Densities are given at cell
  !! centers, velocities, diff coeffs and fluxes at cell faces. flux_f(i) is the
  !! flux between i and i+1. The diffusion coefficients have to be given scaled
  !! by 1 / dx.
  subroutine TS_dd_1d(dens_c, vel_f, dif_f, flux_f)
    real(dp), intent(in)  :: dens_c(:), vel_f(:), dif_f(:)
    real(dp), intent(inout) :: flux_f(:)
    integer :: n, n_cc

    n_cc = size(dens_c)
    ! if (size(vel_f) /= n_cc-1 .or. size(dif_f) /= n_cc-1 .or. size(flux_f) /= n_cc-1) &
    !      call ERR_show("TS_dd_1d: arguments have incompatible sizes")

    do n = 1, n_cc-1
       flux_f(n) = vel_f(n) * 0.5_dp * (dens_c(n)+ dens_c(n+1))

       ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
       flux_f(n)    = flux_f(n) - dif_f(n) * (dens_c(n+1) - dens_c(n))
    end do

  end subroutine TS_dd_1d

  !> Upwind scheme. Densities are given at cell centers, velocities, diff coeffs
  !! and fluxes at cell faces. flux_f(i) is the flux between i and i+1.
  !! The diffusion coefficients have to be given scaled by 1 / dx.
  subroutine TS_dd_1d_up(dens_c, vel_f, dif_f, flux_f)
    real(dp), intent(in)  :: dens_c(:), vel_f(:), dif_f(:)
    real(dp), intent(inout) :: flux_f(:)
    integer :: n, n_cc

    n_cc = size(dens_c)
    ! if (size(vel_f) /= n_cc-1 .or. size(dif_f) /= n_cc-1 .or. size(flux_f) /= n_cc-1) &
    !      call ERR_show("TS_dd_1d_up: arguments have incompatible sizes")

    do n = 1, n_cc-1

       if (vel_f(n) < 0.0_dp) then
          flux_f(n) = vel_f(n) * dens_c(n+1)
       else
          flux_f(n) = vel_f(n) * dens_c(n)
       end if

       ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
       flux_f(n)    = flux_f(n) - dif_f(n) * (dens_c(n+1) - dens_c(n))
    end do

  end subroutine TS_dd_1d_up

  !> Upwind scheme with Koren flux limiter. Densities are given at cell centers,
  !! velocities, diff coeffs and fluxes at cell faces. flux_f(i) is the flux between i and i+1.
  !! The diffusion coefficients have to be given scaled by 1 / dx.
  !! See http://oai.cwi.nl/oai/asset/2269/2269A.pdf
  subroutine TS_dd_1d_up_fl(dens_c, vel_f, dif_f, flux_f)
    real(dp), intent(in)  :: dens_c(:), vel_f(:), dif_f(:)
    real(dp), intent(inout) :: flux_f(:)

    integer :: n, n_cc
    real(dp) :: theta, grad(size(dens_c)-1)

    n_cc = size(dens_c)
    ! if (size(vel_f) /= n_cc-1 .or. size(dif_f) /= n_cc-1 .or. size(flux_f) /= n_cc-1) &
    !      call ERR_show("TS_dd_1d_3ufl error: arguments have incompatible sizes")

    do n = 1, n_cc-1
       grad(n) = dens_c(n+1) - dens_c(n)
    end do

    ! Approximation at left and right boundary with first order upwind
    do n = 1, n_cc-1, n_cc-2
       if (vel_f(n) < 0.0_dp) then
          flux_f(n) = vel_f(n) * dens_c(n+1) - dif_f(n) * grad(n)
       else
          flux_f(n) = vel_f(n) * dens_c(n) - dif_f(n) * grad(n)
       end if

       ! This is a centered difference approximation
       ! flux_f(n) = vel_f(n) * 0.5_dp * (dens_c(n) + dens_c(n+1)) - dif_f(n) * grad(n)
    end do

    do n = 2, n_cc - 2

       ! Advective part with flux limiter
       if (vel_f(n) < 0.0_dp) then
          theta     = limiter_grad_ratio(grad(n), grad(n+1))
          flux_f(n) = vel_f(n) * (dens_c(n+1) - limiter_koren(theta) * grad(n+1))
       else
          theta     = limiter_grad_ratio(grad(n), grad(n-1))
          flux_f(n) = vel_f(n) * (dens_c(n) + limiter_koren(theta) * grad(n-1))
       end if

       ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
       flux_f(n) = flux_f(n) - dif_f(n) * grad(n)
    end do

  end subroutine TS_dd_1d_up_fl

  subroutine TS_dd_1d_up_fl_test(dens_c, vel_f, dif_f, flux_f)
     real(dp), intent(in)  :: dens_c(:), vel_f(:), dif_f(:)
     real(dp), intent(inout) :: flux_f(:)

     integer :: n, n_cc
     real(dp) :: theta, grad(size(dens_c)-1)

     n_cc = size(dens_c)

     do n = 1, n_cc-1
        grad(n) = dens_c(n+1) - dens_c(n)
     end do

     ! Approximation at left and right boundary with first order upwind
     do n = 1, n_cc-1, n_cc-2
        if (vel_f(n) < 0.0_dp) then
           flux_f(n) = vel_f(n) * dens_c(n+1) - dif_f(n) * grad(n)
        else
           flux_f(n) = vel_f(n) * dens_c(n) - dif_f(n) * grad(n)
        end if
     end do

     do n = 2, n_cc - 2
        ! Advective part with flux limiter
        if (vel_f(n) < 0.0_dp) then
           theta     = limiter_grad_ratio(grad(n+1), grad(n))
           flux_f(n) = vel_f(n) * (dens_c(n+1) - limiter_koren(theta) * grad(n))
        else
           theta     = limiter_grad_ratio(grad(n-1), grad(n))
           flux_f(n) = vel_f(n) * (dens_c(n) + limiter_koren(theta) * grad(n))
        end if

        ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
        flux_f(n) = flux_f(n) - dif_f(n) * grad(n)
     end do

  end subroutine TS_dd_1d_up_fl_test

  elemental function limiter_koren(theta)
    real(dp), intent(in) :: theta
    real(dp)             :: limiter_koren
    real(dp), parameter  :: one_sixth = 1.0_dp / 6.0_dp
    limiter_koren = max(0.0d0, min(1.0_dp, theta, (1.0_dp + 2.0_dp * theta) * one_sixth))
  end function limiter_koren

  elemental function limiter_grad_ratio(numerator, denominator)
    real(dp), intent(in) :: numerator, denominator
    real(dp)             :: limiter_grad_ratio
    real(dp), parameter  :: eps = epsilon(1.0d0)
    ! Avoid division by zero, and ensure that at zero gradients we have a ratio of 1
    limiter_grad_ratio = (sign(eps, numerator) + numerator) / (denominator + sign(eps, denominator))
  end function limiter_grad_ratio

end module m_transport_schemes
