!> Module to compute electric fields in 1D
! In 1D the electric field is simply given by
! E_i = E_0 - Sum of charge before current point / epsilon0,

module m_efield_1d
  use m_phys_domain
  use m_init_cond_1d
  use m_units_constants
  
  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  real(dp)              :: pot_left ! potential 0 on the right
  real(dp), allocatable :: EF_values(:)
  logical               :: EF_is_constant
  integer               :: iz_d

  public :: EF_initialize
  public :: EF_compute
  public :: EF_compute_and_get
  public :: EF_compute_and_get_st
  public :: EF_get_at_pos
  public :: EF_get_values
  public :: EF_get_values_st
  public :: EF_get_min_field
  public :: pot_left

contains

  subroutine EF_initialize(cfg)
    use m_config
    type(CFG_t), intent(in) :: cfg

    call CFG_get(cfg, "pot_left", pot_left)
    call CFG_get(cfg, "sim_constant_efield", EF_is_constant)

    ! The electric field is defined at cell faces, so it includes one extra point.
    ! e.g., if the charge density is defined at 0, dx, 2*dx, then the electric field is
    ! defined at -0.5*dx, 0.5*dx, 1.5*dx, 2.5*dx, with the first value equal to the applied field.
    ! These extra values are mostly useful as a buffer for EF_get_at_pos
    allocate(EF_values(PD_grid_size+1))
    EF_values = 0.0_dp
    iz_d     = int(INIT_DI/PD_dx+1) 
  end subroutine EF_initialize

  subroutine EF_compute(net_charge, surface_charge)
    use m_units_constants
    real(dp), intent(in) :: net_charge(:), surface_charge
    real(dp)             :: conv_fac, E_corr_gas, E_corr_eps, E_corr_homog
    integer              :: iz

    if (size(net_charge) /= PD_grid_size) then
       print *, "EF_compute: argument has wrong size"
       stop
    end if


    EF_values = 0.0_dp
    conv_fac = PD_dx / UC_eps0
    do iz = 2, iz_d+1
      EF_values(iz) = EF_values(iz-1) + net_charge(iz-1) * conv_fac
    end do
    EF_values(iz_d+1) = (EF_values(iz_d+1) + surface_charge/UC_eps0)/eps_DI
    do iz = iz_d+2, PD_grid_size+1
      EF_values(iz) = EF_values(iz-1)
    end do
    
    E_corr_homog = (pot_left - sum(EF_values(:))*PD_dx)/PD_length 
    
    E_corr_gas = E_corr_homog / (INIT_DI + (1-INIT_DI)/eps_DI)
    E_corr_eps = E_corr_gas / eps_DI    
    
    EF_values(1:iz_d) = EF_values(1:iz_d) + E_corr_gas
    EF_values(iz_d+1:PD_grid_size+1) = EF_values(iz_d+1:PD_grid_size+1) + E_corr_eps
    

    

  end subroutine EF_compute

  subroutine EF_compute_and_get_st(net_charge, out_efield, surface_charge)
    real(dp), intent(in) :: net_charge(:), surface_charge
    real(dp), intent(out) :: out_efield(:)
    call EF_compute(net_charge, surface_charge)
    call EF_get_values_st(out_efield)
  end subroutine EF_compute_and_get_st

  subroutine EF_compute_and_get(net_charge, out_efield, surface_charge)
    real(dp), intent(in) :: net_charge(:), surface_charge
    real(dp), intent(out) :: out_efield(:)
    call EF_compute(net_charge, surface_charge)
    call EF_get_values(out_efield, surface_charge)
  end subroutine EF_compute_and_get

  !> Get the electric field at a position in the domain (useful for the particle model)
  real(dp) function EF_get_at_pos(pos)
    real(dp), intent(in) :: pos
    real(dp) :: Efield_pos, temp
    integer :: lowIx

    ! EF_values(1) is defined at -0.5 * PD_dx
    lowIx = nint(pos * PD_inv_dx) + 1
    lowIx = min(PD_grid_size, max(1, lowIx))

    Efield_pos = (lowIx - 1.5_dp) * PD_dx
    temp = (pos - Efield_pos) * PD_inv_dx

    ! Do linear interpolation between lowIx and lowIx + 1 in the Efield array, given the position
    EF_get_at_pos = (1.0_dp - temp) * EF_values(lowIx) + temp * EF_values(lowIx+1)

  end function EF_get_at_pos

  !> Get a copy of the electric field at cell centers
  subroutine EF_get_values(out_efield, surface_charge)
    real(dp), intent(in)  :: surface_charge
    real(dp), intent(out) :: out_efield(:)
    out_efield(:) = 0.5_dp * (EF_values(1:PD_grid_size) + EF_values(2:PD_grid_size+1))
    out_efield(iz_d) = 0.5_dp * EF_values(iz_d) + 0.5_dp * (eps_DI * EF_values(iz_d+1) - surface_charge/UC_eps0)
  end subroutine EF_get_values

  real(dp) function EF_get_min_field()
     EF_get_min_field = minval(abs(EF_values))
  end function EF_get_min_field

  !> Get a copy of the electric field at cell faces (interior ones)
  subroutine EF_get_values_st(out_efield)
    real(dp), intent(out) :: out_efield(:)
    out_efield(:) = EF_values(2:PD_grid_size) ! Return only the interior points
  end subroutine EF_get_values_st

end module m_efield_1d
