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
  integer, protected, public  :: domain_nx     = 500
  real(dp), protected, public :: domain_dx     = 2e-6_dp
  real(dp), protected, public :: domain_inv_dx = 5e5_dp

  !> The applied voltage
  real(dp), protected :: applied_voltage = -2e4_dp

  !> Whether a dielectric is present on the left and right
  logical, protected, public :: dielectric_present(2) = [.true., .true.]

  !> Width of the dielectric(s) (outside the computational domain)
  real(dp), protected, public :: dielectric_width = 1.0e-3_dp

  !> Relative permittivity of the dielectric(s)
  real(dp), protected, public :: dielectric_eps = 2.0_dp

  !> Electric field at cell centers, so field_cc(1) is centered at x = 0.5 dx
  real(dp), allocatable, protected, public :: field_cc(:)

  !> Electric field at cell faces, so field_fc(1) is centered at x = 0.0 dx
  real(dp), allocatable, protected, public :: field_fc(:)

  !> Electric potential at cell centers
  real(dp), allocatable, protected, public :: potential(:)

  integer, parameter :: init_none_type     = 0
  integer, parameter :: init_gaussian_type = 1
  integer, parameter :: init_block_type    = 2

  !> Type of initial condition (init_gaussian_type, init_block_type)
  integer :: init_cond_type = init_block_type

  !> Initial condition location (as coordinate)
  real(dp)            :: init_location = 1.0e-3_dp

  !> Initial condition width
  real(dp)            :: init_width = 3.0e-4_dp

  !> Initial condition charge type (1: positive ions, -1: electrons, 0: neutral)
  real(dp)            :: init_charge_type = 0

  !> Initial condition density
  real(dp)            :: init_density = 1.0e15_dp

  !> Initial background density
  real(dp)            :: init_background_density = 0.0e9_dp

  !> Initial energy for particles
  real(dp)            :: init_energy = 1.0_dp

  !> Name of gas mixture
  character(len=40), protected, public :: gas_name = "N2"

  !> Gas pressure in bar
  real(dp), protected, public :: gas_pressure = 1.0_dp

  !> Gas temperature in Kelvin
  real(dp), protected, public :: gas_temperature = 300.0_dp

  !> Gas number density in m^-3
  real(dp), protected, public :: gas_number_dens

  !> Whether ion transport is included
  logical, protected, public :: ion_transport = .true.

  !> Mobility for ions
  real(dp), protected, public :: ion_mobility = 5e-4_dp

  !> Diffusion coefficient for ions
  real(dp), protected, public :: ion_diffusion = 1e-4_dp

  !> Secondary emission yield for ion impact
  real(dp), protected, public :: ion_secondary_emission_yield = 1.0e-2_dp

  !> Number of secondary-emission generating photons per ionization
  real(dp), protected, public :: se_photons_per_ionization = 0.0_dp

  public :: generic_initialize
  public :: compute_field
  public :: get_field_at

  public :: init_pos_ion_dens
  public :: init_elec_dens
  public :: init_energy_dens

  public :: get_flux_1d

contains

  subroutine generic_initialize(cfg)
    use m_config
    use m_units_constants
    type(CFG_t), intent(inout) :: cfg

    character(len=40) :: init_type_name
    character(len=40) :: model_name

    call CFG_add_get(cfg, "gas%name", gas_name, &
         "Name of gas mixture")
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

    call CFG_add_get(cfg, "field%voltage", applied_voltage, &
         "Voltage difference (V) over domain (including dielectrics)")
    call CFG_add_get(cfg, "dielectric%present", dielectric_present, &
         "Whether a dielectric is present on the left and right")
    call CFG_add_get(cfg, "dielectric%width", dielectric_width, &
         "Width of the dielectric(s) in (m)")
    call CFG_add_get(cfg, "dielectric%eps", dielectric_eps, &
         "Relative permittivity of the dielectric(s)")

    init_type_name = "block"
    call CFG_add_get(cfg, "init%type", init_type_name, &
         "Type of initial condition (none, block, gaussian)")

    select case (init_type_name)
    case ("none")
       init_cond_type = init_none_type
    case ("block")
       init_cond_type = init_block_type
    case ("gaussian", "Gaussian")
       init_cond_type = init_gaussian_type
    case default
       print *, "Wrong init%type:", trim(init_type_name)
       print *, "Options are: none, block, gaussian"
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
    call CFG_add_get(cfg, "init%background_density", &
         init_background_density, "Initial background density (1/m^3)")

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

  subroutine compute_field(net_charge, surface_charge, time)
    use m_units_constants
    real(dp), intent(in) :: net_charge(:)
    real(dp), intent(in) :: surface_charge(2)
    real(dp), intent(in) :: time
    real(dp)             :: conv_fac, E_corr_gas
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
            (field_fc(nx+1) + surface_charge(2))
    end if

    ! Have to correct the potential by this amount
    pot_correct = get_potential(time) - pot_diff

    ! Determine the correction between the dielectrics
    E_corr_gas = -pot_correct / (domain_length + &
         count(dielectric_present) * dielectric_width/dielectric_eps)

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

  !> Compute advective and diffusive flux
  subroutine get_flux_1d(cc, v, dc, dx, flux, nc, ngc)
    integer, intent(in)   :: nc               !< Number of cells
    integer, intent(in)   :: ngc              !< Number of ghost cells
    real(dp), intent(in)  :: cc(1-ngc:nc+ngc) !< Cell-centered values
    !> Input: velocities at cell faces
    real(dp), intent(in)  :: v(1:nc+1)
    !> Input: diffusion coefficients at cell faces
    real(dp), intent(in)  :: dc(1:nc+1)
    !> Grid spacing
    real(dp), intent(in)  :: dx
    !> Output: flux at cell faces
    real(dp), intent(out) :: flux(1:nc+1)
    real(dp)              :: gradp, gradc, gradn, inv_dx
    integer               :: n

    inv_dx = 1/dx

    do n = 1, nc+1
       gradc = cc(n) - cc(n-1)  ! Current gradient
       if (v(n) < 0.0_dp) then
          gradn = cc(n+1) - cc(n) ! Next gradient
          flux(n) = v(n) * (cc(n) - koren_mlim(gradc, gradn))
       else                     ! v(n) > 0
          gradp = cc(n-1) - cc(n-2) ! Previous gradient
          flux(n) = v(n) * (cc(n-1) + koren_mlim(gradc, gradp))
       end if
       ! Add diffusive flux (central differences)
       flux(n) = flux(n) - dc(n) * gradc * inv_dx
    end do

  end subroutine get_flux_1d

  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  elemental function koren_mlim(a, b) result(bphi)
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: sixth = 1/6.0_dp
    real(dp)             :: bphi, aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       bphi = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value a/b
       bphi = a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/6
       bphi = sixth * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 1
       bphi = b
    end if
  end function koren_mlim

end module m_generic
