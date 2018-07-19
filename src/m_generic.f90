!> Module with generic functionality for particle and fluid models, in
!> particular a description of the domain and routines to compute the electric
!> field.
module m_generic

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter, public :: model_fluid    = 1
  integer, parameter, public :: model_particle = 2
  integer, protected, public :: model_type     = 1

  real(dp), protected, public :: domain_length = 2.0e-3_dp
  integer, protected, public  :: domain_nx  = 500
  real(dp), protected, public :: domain_dx     = 2e-6_dp
  real(dp), protected, public :: domain_inv_dx = 5e5_dp

  !> The applied voltage
  real(dp), protected :: applied_voltage = 1e3_dp

  !> Whether a dielectric is present on the left and right
  logical, protected, public :: dielectric_present(2)

  !> Width of the dielectric(s) (outside the computational domain)
  real(dp), protected, public :: dielectric_width = 1.0e-3_dp

  !> Relative permittivity of the dielectric(s)
  real(dp), protected, public :: dielectric_eps = 2.0_dp

  !> Electric field at cell centers, so field_cc(1) is centered at x = 0.5 dx
  real(dp), allocatable, protected, public :: field_cc(:)

  !> Electric field at cell faces, so field_fc(1) is centered at x = 0.0 dx
  real(dp), allocatable, protected, public :: field_fc(:)

  integer, parameter :: init_gaussian_type = 1
  integer, parameter :: init_block_type    = 2

  !> Type of initial condition (init_gaussian_type, init_block_type)
  integer :: init_cond_type = init_block_type

  !> Initial condition location (as coordinate)
  real(dp)            :: init_location = 1.0e-3_dp

  !> Initial condition width
  real(dp)            :: init_width = 1.0e-4_dp

  !> Initial condition charge type (1: positive ions, -1: electrons, 0: neutral)
  real(dp)            :: init_charge_type = 0

  !> Initial condition density
  real(dp)            :: init_density = 1.0e9_dp

  !> Initial condition density
  real(dp)            :: init_energy = 1.0_dp

  !> Name of gas mixture
  character(len=40), protected, public :: gas_name = "N2"

  !> Gas pressure in bar
  real(dp), protected, public :: gas_pressure = 1.0_dp

  !> Gas temperature in Kelvin
  real(dp), protected, public :: gas_temperature = 300.0_dp

  public :: generic_initialize
  public :: compute_field
  public :: get_field_at

  public :: init_pos_ion_dens
  public :: init_elec_dens
  public :: init_energy_dens
  ! public :: field_emission
  ! public :: pot_left
  ! public :: work_fun

contains

  subroutine generic_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg

    character(len=40) :: init_type_name
    character(len=40) :: model_type_name
    character(len=40) :: dielectric_side

    call CFG_add_get(cfg, "gas%name", gas_name, &
         "Name of gas mixture")
    call CFG_add_get(cfg, "gas%pressure", gas_pressure, &
         "Gas pressure in bar")
    call CFG_add_get(cfg, "gas%temperature", gas_temperature, &
         "Gas temperature in Kelvin")

    call CFG_add_get(cfg, "domain%length", domain_length, &
         "Length of the domain")
    call CFG_add_get(cfg, "domain%n_cell", domain_nx, &
         "Number of cells in the domain")

    domain_dx        = domain_length / domain_nx
    domain_inv_dx    = 1/domain_dx

    model_type_name = "fluid"
    call CFG_add_get(cfg, "model_type", model_type_name, &
         "Type of model. Can be particle or fluid")

    select case (model_type_name)
    case ("part", "particle")
       model_type = model_particle
    case ("fluid")
       model_type = model_fluid
    case default
       print *, "Invalid simulation type given: ", MODEL_type_name
       print *, "Supported are: particle, fluid"
       error stop
    end select

    allocate(field_cc(domain_nx))
    field_cc(:) = 0.0_dp

    ! Field at cell faces includes one extra point
    allocate(field_fc(domain_nx+1))
    field_fc(:) = 0.0_dp

    call CFG_add_get(cfg, "field%voltage", applied_voltage, &
         "Voltage difference (V) over domain (including dielectrics)")

    dielectric_present(:) = [.false., .false.]
    call CFG_add_get(cfg, "dielectric%present", dielectric_present, &
         "Whether a dielectric is present on the left and right")
    call CFG_add_get(cfg, "dielectric%width", dielectric_width, &
         "Width of the dielectric(s) in (m)")
    call CFG_add_get(cfg, "dielectric%eps", dielectric_eps, &
         "Relative permittivity of the dielectric(s)")

    init_type_name = "block"
    call CFG_add_get(cfg, "init%type", init_type_name, &
         "Type of initial condition (block or gaussian)")

    select case (init_type_name)
    case ("block")
       init_cond_type = init_block_type
    case ("gaussian", "Gaussian")
       init_cond_type = init_gaussian_type
    case default
       print *, "Wrong init%type:", trim(init_type_name)
       print *, "Options are: block, gaussian"
       error stop
    end select
  end subroutine generic_initialize

  subroutine compute_field(net_charge, surface_charge, time)
    use m_units_constants
    real(dp), intent(in) :: net_charge(:)
    real(dp), intent(in) :: surface_charge(2)
    real(dp), intent(in) :: time
    real(dp)             :: conv_fac, E_corr_gas, E_corr_eps
    real(dp)             :: pot_diff, pot_correct
    integer              :: n

    if (size(net_charge) /= domain_nx) &
         error stop "compute_field: argument has wrong size"

    conv_fac = domain_dx / UC_eps0

    ! Starting guess for electric field on the left
    if (dielectric_present(1)) then
       field_fc(1) = surface_charge(1) / UC_eps0
    else
       field_fc(1) = 0.0_dp
    end if

    ! Handle the interior region
    do n = 2, domain_nx+1
       field_fc(n) = field_fc(n-1) + net_charge(n-1) * conv_fac
    end do

    ! Compute total potential difference in the domain, with a weight of 0.5 for
    ! the boundary field values
    pot_diff = -domain_dx * (sum(field_fc(2:domain_nx)) + &
         0.5_dp * (field_fc(1) + field_fc(domain_nx+1)))

    if (dielectric_present(2)) then
       ! Add contribution of field in right dielectric
       pot_diff = pot_diff - dielectric_width / dielectric_eps * &
            (field_fc(domain_nx+1) + surface_charge(2))
    end if

    ! Have to correct the potential by this amount
    pot_correct = get_potential(time) - pot_diff

    ! Determine the correction between the dielectrics
    E_corr_gas = -pot_correct / (domain_length + &
         count(dielectric_present) * dielectric_width/dielectric_eps)

    field_fc(:) = field_fc(:) + E_corr_gas

    ! Average to get cell-centered field
    do n = 1, domain_nx
       field_cc(n) = 0.5_dp * (field_fc(n) + field_fc(n+1))
    end do

    ! print *, "total pot", -sum(field_cc) * domain_dx
  end subroutine compute_field

  function get_potential(time) result(pot)
    real(dp), intent(in) :: time
    real(dp)             :: pot

    pot = applied_voltage
  end function get_potential

  !> Get the electric field at a position in the domain (useful for the particle
  !> model)
  function get_field_at(x) result(field)
    real(dp), intent(in) :: x
    real(dp)             :: field, tmp
    integer              :: low_ix

    ! Do linear interpolation, note that field_fc(1) is defined at 0.0
    tmp    = x * domain_inv_dx
    low_ix = floor(tmp) + 1
    ! Coefficient for upper point between 0 and 1
    tmp    = low_ix - tmp

    if (low_ix < 1 .or. low_ix > domain_nx) then
       print *, "get_field_at invalid position x = ", x
       error stop
    end if

    field = (1.0_dp - tmp) * field_fc(low_ix) + tmp * field_fc(low_ix+1)
  end function get_field_at

  ! real(dp) function field_emission(fld)
  ! use m_units_constants
  !   real(dp), intent(in) :: fld
  !   real(dp)             :: W_ev, A, T

  !   W_ev = sqrt(UC_elem_charge * fld/(4*UC_pi*UC_eps0))
  !   A = 120 * (1.0_dp/UC_elem_charge) * 10000
  !   T = 270_dp

  !   field_emission = A * T**2 * exp(-UC_elem_charge*(work_fun-W_ev)/(UC_boltzmann_const*T))
  ! end function field_emission

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
    case default
       ! This should never occur, as input is checked before
       dens = 0.0_dp
    end select
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

end module m_generic
