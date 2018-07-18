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
  integer, protected, public  :: domain_ncell  = 500
  real(dp), protected, public :: domain_dx     = 2e-6_dp
  real(dp), protected, public :: domain_inv_dx = 5e5_dp

  !> The applied voltage
  real(dp), protected :: applied_voltage = 1e3_dp

  !> Whether a dielectric is present on the left
  logical, protected, public :: left_dielectric_present = .false.

  !> The boundary coordinate of the left dielectric
  real(dp), protected, public :: left_dielectric_x = -1.0_dp

  !> Relative permittivity of the left dielectric
  real(dp), protected, public :: left_dielectric_eps = 2.0_dp

  !> Cell face index of the dielectric boundary, where 1 is at 0.0
  integer, protected, public :: left_dielectric_iface = -1

  !> Electric field at cell centers
  real(dp), allocatable, protected, public :: field_cc(:)

  !> Electric field at cell faces
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

    call CFG_add_get(cfg, "gas%name", gas_name, &
         "Name of gas mixture")
    call CFG_add_get(cfg, "gas%pressure", gas_pressure, &
         "Gas pressure in bar")
    call CFG_add_get(cfg, "gas%temperature", gas_temperature, &
         "Gas temperature in Kelvin")

    call CFG_add_get(cfg, "domain%length", domain_length, &
         "Length of the domain")
    call CFG_add_get(cfg, "domain%n_cell", domain_ncell, &
         "Number of cells in the domain")

    domain_dx        = domain_length / domain_ncell
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
       error stop "Invalid simulation type"
    end select

    allocate(field_cc(domain_ncell))
    field_cc(:) = 0.0_dp

    ! Field at cell faces includes one extra point
    allocate(field_fc(domain_ncell+1))
    field_fc(:) = 0.0_dp

    ! field_cc(1) is centered at x = 0.5 dx
    ! field_fc(1) is centered at x = 0.0 dx

    call CFG_add_get(cfg, "field%voltage", applied_voltage, &
         "Voltage (V) at right side of domain")

    call CFG_add_get(cfg, "left_dielectric%x", left_dielectric_x, &
         "Interface position of left dielectric")
    call CFG_add_get(cfg, "left_dielectric%eps", left_dielectric_eps, &
         "Relative permittivity of left dielectric")

    left_dielectric_present = (&
         left_dielectric_x >= 0.0_dp .and. &
         left_dielectric_x < domain_length .and. &
         abs(left_dielectric_eps - 1.0_dp) > 0.0_dp)

    if (left_dielectric_present) then
       left_dielectric_iface = nint(left_dielectric_x / domain_dx) + 1
       left_dielectric_x = (left_dielectric_iface-1) * domain_dx
    end if

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

  subroutine compute_field(net_charge, left_surface_charge, time)
    use m_units_constants
    real(dp), intent(in) :: net_charge(:), left_surface_charge, time
    real(dp)             :: conv_fac, E_corr_gas, E_corr_eps
    real(dp)             :: pot_diff, pot_correct
    integer              :: n, iface

    if (size(net_charge) /= domain_ncell) &
         error stop "compute_field: argument has wrong size"

    conv_fac = domain_dx / UC_eps0

    if (.not. left_dielectric_present) then
       field_fc(1) = 0.0_dp
       do n = 2, domain_ncell+1
          field_fc(n) = field_fc(n-1) + net_charge(n-1) * conv_fac
       end do

       ! Compute total potential difference. Here we assign a weight of 0.5 to the
       ! boundary values (on the dielectric and the wall)
       pot_diff = -domain_dx * (sum(field_fc(2:domain_ncell)) + &
            0.5_dp * (field_fc(1) + field_fc(domain_ncell+1)))

       ! Have to correct the potential by this amount
       pot_correct = get_potential(time) - pot_diff

       E_corr_gas = -pot_correct / domain_length
       field_fc = field_fc + E_corr_gas

       do n = 1, domain_ncell
          field_cc(n) = 0.5_dp * (field_fc(n) + field_fc(n+1))
       end do
    else
       iface = left_dielectric_iface
       ! First we start from a guess field_fc(1) = 0, which is the electric field
       ! at the left domain boundary. Inside the dielectric, the field is constant.
       field_fc(1:iface-1) = 0.0_dp

       ! Note that on the dielectric boundary, the electric field is different on
       ! both sides, but we only store the value outside the dielectric for now.
       field_fc(iface) = left_surface_charge / UC_eps0

       ! Handle the region outside the dielectric
       do n = iface+1, domain_ncell+1
          field_fc(n) = field_fc(n-1) + net_charge(n-1) * conv_fac
       end do

       ! Compute total potential difference. Here we assign a weight of 0.5 to the
       ! boundary values (on the dielectric and the wall)
       pot_diff = -domain_dx * (sum(field_fc(iface+1:domain_ncell)) + &
            0.5_dp * (field_fc(iface) + field_fc(domain_ncell+1)))

       ! Have to correct the potential by this amount
       pot_correct = get_potential(time) - pot_diff

       ! Correction in gas phase
       E_corr_eps = -pot_correct / (left_dielectric_x + &
            (domain_length-left_dielectric_x) * left_dielectric_eps)
       E_corr_gas = E_corr_eps * left_dielectric_eps

       field_fc(1:iface-1) = field_fc(1:iface-1) + E_corr_eps
       field_fc(iface:domain_ncell+1) = &
            field_fc(iface:domain_ncell+1) + E_corr_gas

       ! Determine cell-centered field
       field_cc(1:iface-1) = E_corr_eps
       do n = iface, domain_ncell
          field_cc(n) = 0.5_dp * (field_fc(n) + field_fc(n+1))
       end do
    end if

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

    if (low_ix < 1 .or. low_ix > domain_ncell) then
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

    if (x > left_dielectric_x) then
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
    else
       dens = 0.0_dp
    end if
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
