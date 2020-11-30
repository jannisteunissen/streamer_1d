!> Module to set the initial conditions
module m_init_cond

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter :: init_none_type        = 0
  integer, parameter :: init_gaussian_type    = 1
  integer, parameter :: init_block_type       = 2
  integer, parameter :: init_exponential_type = 3

  !> Type of initial condition
  integer :: init_cond_type = init_block_type

  !> Initial condition location (as coordinate)
  real(dp) :: init_location = 1.0e-3_dp

  !> Initial condition width
  real(dp) :: init_width = 1.5e-4_dp

  !> Initial condition charge type (1: positive ions, -1: electrons, 0: neutral)
  real(dp) :: init_charge_type = 0

  !> Initial condition density
  real(dp) :: init_density = 1.0e15_dp

  !> Initial background density
  real(dp) :: init_background_density = 0.0e9_dp

  !> Initial energy for particles
  real(dp), protected, public :: init_energy = 0.0_dp

  public :: init_initialize
  public :: init_pos_ion_dens
  public :: init_elec_dens
  public :: init_energy_dens

contains

  subroutine init_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg
    character(len=40) :: init_type_name

    init_type_name = "block"
    call CFG_add_get(cfg, "init%type", init_type_name, &
         "Type of initial condition (none, block, gaussian, exp)")

    select case (init_type_name)
    case ("none")
       init_cond_type = init_none_type
    case ("block")
       init_cond_type = init_block_type
    case ("gaussian", "Gaussian")
       init_cond_type = init_gaussian_type
    case ("exp", "exponential")
       init_cond_type = init_exponential_type
    case default
       print *, "Wrong init%type:", trim(init_type_name)
       print *, "Options are: none, block, gaussian, exp"
       error stop
    end select

    call CFG_add_get(cfg, "init%width", init_width, &
         "Initial condition width (m)")
    call CFG_add_get(cfg, "init%charge_type", init_charge_type, &
         "Initial condition charge type (1: p.ions, -1: e-, 0: neutral)")
    call CFG_add_get(cfg, "init%density", init_density, &
         "Initial condition density (1/m^3)")
    call CFG_add_get(cfg, "init%location", init_location, &
         "Initial condition location (m)")
    call CFG_add_get(cfg, "init%energy_eV", init_energy, &
         "Initial particle energy (eV)")
    call CFG_add_get(cfg, "init%background_density", &
         init_background_density, "Initial background density (1/m^3)")
  end subroutine init_initialize

  elemental function get_dens(x) result(dens)
    real(dp), intent(in) :: x
    real(dp)             :: dens

    select case (init_cond_type)
    case (init_gaussian_type)
       dens = init_density * &
            exp(-(x - init_location)**2 / (2 * init_width**2))
    case (init_block_type)
       if (abs(x - init_location) < init_width) then
          dens = init_density
       else
          dens = 0.0_dp
       end if
    case (init_exponential_type)
       dens = init_density * exp(-abs(x - init_location)/init_width)
    case default
       dens = 0.0_dp
    end select

    dens = dens + init_background_density
  end function get_dens

  elemental function init_elec_dens(x) result(elec_dens)
    real(dp), intent(in) :: x
    real(dp)             :: elec_dens

    if (init_charge_type < 1) then
       elec_dens = get_dens(x)
    else
       elec_dens = 0.0_dp
    end if
  end function init_elec_dens

  elemental function init_pos_ion_dens(x) result(pos_ion_dens)
    real(dp), intent(in) :: x
    real(dp)             :: pos_ion_dens

    if (init_charge_type > -1) then
       pos_ion_dens = get_dens(x)
    else
       pos_ion_dens = 0.0_dp
    end if
  end function init_pos_ion_dens

  elemental function init_energy_dens(x) result(energy_dens)
    real(dp), intent(in) :: x
    real(dp)             :: energy_dens

    if (init_charge_type < 1) then
       energy_dens = get_dens(x) * init_energy
    else
       energy_dens = 0.0_dp
    end if
  end function init_energy_dens

end module m_init_cond
