!> Module to compute electric fields in 1D
! In 1D the electric field is simply given by
! E_i = E_0 - Sum of charge before current point / epsilon0,

module m_efield_1d

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer               :: EF_grid_size
  real(dp)              :: EF_delta_x, EF_inv_delta_x
  real(dp)              :: EF_applied_field
  real(dp), allocatable :: EF_values(:)
  logical               :: EF_is_constant

  public :: EF_initialize
  public :: EF_compute
  public :: EF_compute_and_get
  public :: EF_compute_and_get_st
  public :: EF_get_at_pos
  public :: EF_get_values
  public :: EF_get_values_st
  public :: EF_get_min_field

contains

  subroutine EF_initialize(cfg)
    use m_config
    type(CFG_t), intent(in) :: cfg
    
    call CFG_get(cfg, "grid_num_points", EF_grid_size)
    call CFG_get(cfg, "grid_delta_x", EF_delta_x)
    EF_inv_delta_x = 1.0_dp / EF_delta_x
    call CFG_get(cfg, "sim_applied_efield", EF_applied_field)
    call CFG_get(cfg, "sim_constant_efield", EF_is_constant)

    ! The electric field is defined at cell faces, so it includes one extra point.
    ! e.g., if the charge density is defined at 0, dx, 2*dx, then the electric field is
    ! defined at -0.5*dx, 0.5*dx, 1.5*dx, 2.5*dx, with the first value equal to the applied field.
    ! These extra values are mostly useful as a buffer for EF_get_at_pos
    allocate(EF_values(EF_grid_size+1))
    EF_values = EF_applied_field

  end subroutine EF_initialize

  subroutine EF_compute(net_charge)
    use m_units_constants
    real(dp), intent(in) :: net_charge(:)
    real(dp)             :: conv_fac
    integer              :: iz

    if (size(net_charge) /= EF_grid_size) then
       print *, "EF_compute: argument has wrong size"
       stop
    end if

    conv_fac = EF_delta_x / UC_eps0
    if (.not. EF_is_constant) then
       do iz = 2, EF_grid_size+1
          EF_values(iz) = EF_values(iz-1) + net_charge(iz-1) * conv_fac
       end do
    end if

  end subroutine EF_compute

  subroutine EF_compute_and_get_st(net_charge, out_efield)
    real(dp), intent(in) :: net_charge(:)
    real(dp), intent(out) :: out_efield(:)
    call EF_compute(net_charge)
    call EF_get_values_st(out_efield)
  end subroutine EF_compute_and_get_st

  subroutine EF_compute_and_get(net_charge, out_efield)
    real(dp), intent(in) :: net_charge(:)
    real(dp), intent(out) :: out_efield(:)
    call EF_compute(net_charge)
    call EF_get_values(out_efield)
  end subroutine EF_compute_and_get

  !> Get the electric field at a position in the domain (useful for the particle model)
  real(dp) function EF_get_at_pos(pos)
    real(dp), intent(in) :: pos
    real(dp) :: Efield_pos, temp
    integer :: lowIx

    ! EF_values(1) is defined at -0.5 * EF_delta_x
    lowIx = nint(pos * EF_inv_delta_x) + 1
    lowIx = min(EF_grid_size, max(1, lowIx))

    Efield_pos = (lowIx - 1.5_dp) * EF_delta_x
    temp = (pos - Efield_pos) * EF_inv_delta_x

    ! Do linear interpolation between lowIx and lowIx + 1 in the Efield array, given the position
    EF_get_at_pos = (1.0_dp - temp) * EF_values(lowIx) + temp * EF_values(lowIx+1)

  end function EF_get_at_pos

  !> Get a copy of the electric field at cell centers
  subroutine EF_get_values(out_efield)
    real(dp), intent(out) :: out_efield(:)
    out_efield(:) = 0.5_dp * (EF_values(1:EF_grid_size) + EF_values(2:EF_grid_size+1))
  end subroutine EF_get_values

  real(dp) function EF_get_min_field()
     EF_get_min_field = minval(abs(EF_values))
  end function EF_get_min_field

  !> Get a copy of the electric field at cell faces (interior ones)
  subroutine EF_get_values_st(out_efield)
    real(dp), intent(out) :: out_efield(:)
    out_efield(:) = EF_values(2:EF_grid_size) ! Return only the interior points
  end subroutine EF_get_values_st

end module m_efield_1d
