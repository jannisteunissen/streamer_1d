!> Module with generic functionality for particle and fluid models, in
!> particular a description of the domain and routines to compute the electric
!> field.
module m_generic
  use m_init_cond
  use m_flux_scheme

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  integer, parameter :: model_fluid    = 1
  integer, parameter :: model_particle = 2
  integer, protected :: model_type     = 1

  !> Length of domain
  real(dp), protected :: domain_length = 2.0e-3_dp

  !> Number of grid cells
  integer, protected  :: domain_nx     = 500

  !> Grid spacing
  real(dp), protected :: domain_dx

  !> Inverse grid spacing
  real(dp), protected :: domain_inv_dx

  !> Apply E = -voltage/domain_length through Neumann boundary condition
  logical, protected :: voltage_neumann_bc = .false.

  !> Voltage (V) on right domain boundary (left is grounded)
  real(dp), protected :: voltage_v0       = -2e4_dp

  !> Voltage rise time (s)
  real(dp), protected :: voltage_risetime = 0.0_dp

  !> Amplitude of sinusoidal voltage (V)
  real(dp), protected :: voltage_sin_v0   = 0.0_dp

  !> Frequency (1/s) of sinusoidal voltage
  real(dp), protected :: voltage_sin_freq = 1.0e8_dp

  !> Whether a dielectric is present on the left and right
  logical, protected :: dielectric_present(2) = [.true., .true.]

  !> Width of the dielectric(s) (outside the computational domain)
  real(dp), protected :: dielectric_width = 1.0e-3_dp

  !> Relative permittivity of the dielectric(s)
  real(dp), protected :: dielectric_eps = 2.0_dp

  !> Electric field at cell centers, so field_cc(1) is centered at x = 0.5 dx
  real(dp), allocatable, protected :: field_cc(:)

  !> Electric field at cell faces, so field_fc(1) is centered at x = 0.0 dx
  real(dp), allocatable, protected :: field_fc(:)

  !> Electric potential at cell centers
  real(dp), allocatable, protected :: potential(:)

  !> Gas pressure in bar
  real(dp), protected :: gas_pressure = 1.0_dp

  !> Gas temperature in Kelvin
  real(dp), protected :: gas_temperature = 300.0_dp

  !> Gas number density in m^-3
  real(dp), protected :: gas_number_dens

  !> Whether ion transport is included
  logical, protected :: ion_transport = .true.

  !> Mobility for ions
  real(dp), protected :: ion_mobility = 2e-4_dp

  !> Diffusion coefficient for ions
  real(dp), protected :: ion_diffusion = 2e-5_dp

  !> Secondary emission yield for ion impact
  real(dp), protected :: ion_secondary_emission_yield = 1.0e-2_dp

  !> Number of secondary-emission generating photons per ionization
  real(dp), protected :: se_photons_per_ionization = 0.0_dp

contains

  !> Initialize the parameters
  subroutine generic_initialize(cfg)
    use m_config
    use m_units_constants
    type(CFG_t), intent(inout) :: cfg

    character(len=40) :: model_name

    call init_initialize(cfg)

    call CFG_add_get(cfg, "gas%pressure", gas_pressure, &
         "Gas pressure in bar")
    call CFG_add_get(cfg, "gas%temperature", gas_temperature, &
         "Gas temperature in Kelvin")

    gas_number_dens = 1.0e5_dp * gas_pressure / &
         (UC_boltzmann_const * gas_temperature)

    call CFG_add_get(cfg, "domain%length", domain_length, &
         "Length of the domain")
    call CFG_add_get(cfg, "domain%nx", domain_nx, &
         "Number of cells in the domain")

    domain_dx        = domain_length / domain_nx
    domain_inv_dx    = 1/domain_dx

    model_name = "fluid"
    call CFG_add_get(cfg, "model", model_name, &
         "Type of model. Can be particle or fluid")

    select case (model_name)
    case ("part", "particle")
       model_type = model_particle
    case ("fluid")
       model_type = model_fluid
    case default
       print *, "Invalid simulation type given: ", model_name
       print *, "Supported are: particle, fluid"
       error stop
    end select

    allocate(field_cc(domain_nx))
    field_cc(:) = 0.0_dp

    allocate(potential(domain_nx))
    potential(:) = 0.0_dp

    ! Field at cell faces includes one extra point
    allocate(field_fc(domain_nx+1))
    field_fc(:) = 0.0_dp

    call CFG_add_get(cfg, "voltage%neumann_bc", voltage_neumann_bc, &
         "Apply E = -voltage/domain_length through Neumann boundary condition")
    call CFG_add_get(cfg, "voltage%v0", voltage_v0, &
         "Voltage amplitude (V) over domain (including dielectrics)")
    call CFG_add_get(cfg, "voltage%risetime", voltage_risetime, &
         "Voltage rise time (s)")
    call CFG_add_get(cfg, "voltage%sin_v0", voltage_sin_v0, &
         "Sinusoidal voltage amplitude (V)")
    call CFG_add_get(cfg, "voltage%sin_freq", voltage_sin_freq, &
         "Sinusoidal voltage frequency (1/s)")
    call CFG_add_get(cfg, "dielectric%present", dielectric_present, &
         "Whether a dielectric is present on the left and right")
    call CFG_add_get(cfg, "dielectric%width", dielectric_width, &
         "Width of the dielectric(s) in (m)")
    call CFG_add_get(cfg, "dielectric%eps", dielectric_eps, &
         "Relative permittivity of the dielectric(s)")

    call CFG_add_get(cfg, "ion%transport", ion_transport, &
         "Whether ion transport is included")
    call CFG_add_get(cfg, "ion%mobility", ion_mobility, &
         "Mobility for ions (m^2/Vs)")
    call CFG_add_get(cfg, "ion%diffusion", ion_diffusion, &
         "Diffusion coefficient for ions (m^2/s)")
    call CFG_add_get(cfg, "ion%se_yield", ion_secondary_emission_yield, &
         "Secondary emission yield for ions")

    call CFG_add_get(cfg, "photons%se_per_ionization", &
         se_photons_per_ionization, &
         "Number of secondary-emission generating photons per ionization")

  end subroutine generic_initialize

  !> Compute the electric field and store it
  subroutine compute_field(net_charge, surface_charge, time)
    use m_units_constants
    real(dp), intent(in) :: net_charge(:)
    real(dp), intent(in) :: surface_charge(2)
    real(dp), intent(in) :: time
    real(dp)             :: conv_fac, E_applied, E_corr_gas
    real(dp)             :: pot_diff, pot_correct
    integer              :: n, nx

    nx = domain_nx
    conv_fac = domain_dx / UC_eps0

    if (size(net_charge) /= domain_nx) &
         error stop "compute_field: argument has wrong size"

    ! Starting guess for electric field on the left
    if (dielectric_present(1)) then
       field_fc(1) = surface_charge(1) / UC_eps0
    else
       field_fc(1) = 0.0_dp
    end if

    ! Handle the interior region
    do n = 2, nx+1
       field_fc(n) = field_fc(n-1) + net_charge(n-1) * conv_fac
    end do

    ! Compute total potential difference in the domain, with a weight of 0.5 for
    ! the boundary field values
    pot_diff = -domain_dx * (sum(field_fc(2:nx)) + &
         0.5_dp * (field_fc(1) + field_fc(nx+1)))

    if (dielectric_present(2)) then
       ! Add contribution of field in right dielectric
       pot_diff = pot_diff - dielectric_width / dielectric_eps * &
            (field_fc(nx+1) + surface_charge(2) / UC_eps0)
    end if

    if (voltage_neumann_bc) then
       ! Correct field
       E_applied = -get_potential(time)/domain_length
       E_corr_gas = E_applied - field_fc(nx+1)
    else
       ! Have to correct the potential by this amount
       pot_correct = get_potential(time) - pot_diff

       ! Determine the correction between the dielectrics
       E_corr_gas = -pot_correct / (domain_length + &
            count(dielectric_present) * dielectric_width/dielectric_eps)
    end if

    field_fc(:) = field_fc(:) + E_corr_gas

    ! Average to get cell-centered field
    do n = 1, nx
       field_cc(n) = 0.5_dp * (field_fc(n) + field_fc(n+1))
    end do

    ! Compute cell-centered potential
    potential(1) = -0.5_dp * domain_dx * field_fc(1)
    do n = 2, nx - 1
       potential(n) = potential(n-1) - domain_dx * field_fc(n)
    end do
    potential(nx) = potential(nx-1) - 0.5_dp * domain_dx * field_fc(nx)

  end subroutine compute_field

  !> Get the potential at the current time
  function get_potential(time) result(pot)
    real(dp), intent(in) :: time
    real(dp)             :: pot, rise_factor
    real(dp), parameter  :: pi = acos(-1.0_dp)

    if (voltage_risetime > 0) then
       rise_factor = 1 - exp(-time/voltage_risetime)
    else
       rise_factor = 1
    end if

    pot = voltage_v0 * rise_factor + voltage_sin_v0 * &
         sin(2 * pi * voltage_sin_freq * time)
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

end module m_generic
