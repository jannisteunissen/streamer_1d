!> Implementation of a drift-diffusion-reaction plasma fluid model
module m_fluid_1d
  use m_generic
  use m_lookup_table

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter :: n_ghost_cells = 2

  integer, parameter :: forward_euler    = 1
  integer, parameter :: trapezoidal      = 2
  integer, parameter :: rk2              = 3
  integer, parameter :: rk4              = 4
  integer            :: time_step_method = rk2

  integer, parameter :: num_arrays = 3
  integer, parameter :: iv_elec    = 1
  integer, parameter :: iv_pion    = 2
  integer, parameter :: iv_nion    = 3

  integer, parameter :: num_scalars = 4
  integer, parameter :: i_lbound_elec = 1
  integer, parameter :: i_rbound_elec = 2
  integer, parameter :: i_lbound_pion = 3
  integer, parameter :: i_rbound_pion = 4

  logical :: fluid_limit_velocity = .false.

  integer :: if_mob, if_dif, if_src
  integer :: if_att, fluid_num_fld_coef

  type state_t
     real(dp), allocatable :: a(:, :) ! arrays
     real(dp), allocatable :: s(:)    ! scalars
  end type state_t

  type(state_t) :: fluid_state

  type(LT_t) :: fluid_lkp_fld

  !> Maximum electric field for the transport data table
  real(dp) :: fluid_max_field

  !> Set density to zero below this value
  real(dp) :: fluid_small_density

  ! Extrapolation coefficients for transport data
  real(dp) :: extrap_src_dydx
  real(dp) :: extrap_src_y0
  real(dp) :: extrap_mob_c0
  real(dp) :: extrap_mob_c1

  public :: fluid_initialize
  public :: fluid_advance
  public :: fluid_write_output

contains

  !> Initialize the fluid module
  subroutine fluid_initialize(cfg)
    use m_transport_data
    use m_config

    type(CFG_t), intent(inout) :: cfg
    integer, parameter         :: name_len = 100
    character(len=200)         :: input_file
    integer                    :: n, table_size
    real(dp)                   :: xx, x(2), y(2)
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=100)         :: data_name
    character(len=40)          :: integrator

    input_file          = "input/n2_transport_data_siglo.txt"
    table_size          = 1000
    integrator          = "rk2"
    fluid_max_field     = 2.5e7_dp
    fluid_small_density = 1.0_dp

    call CFG_add_get(cfg, "fluid%input_file", input_file, &
         "Input file with cross sections")
    call CFG_add_get(cfg, "fluid%table_size", table_size, &
         "Size of lookup table for transport coefficients")
    call CFG_add_get(cfg, "fluid%table_max_efield", fluid_max_field, &
         "Maximum electric field in transport coefficient table")
    call CFG_add_get(cfg, "fluid%small_density", fluid_small_density, &
         "Smallest allowed density in the fluid model (1/m3)")
    call CFG_add_get(cfg, "fluid%integrator", integrator, &
         "Time integrator (euler, trapezoidal, rk2, rk4)")

    call CFG_add(cfg, "fluid%fld_mob", "efield[V/m]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    call CFG_add(cfg, "fluid%fld_dif", "efield[V/m]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(cfg, "fluid%fld_alpha", "efield[V/m]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(cfg, "fluid%fld_eta", "efield[V/m]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")

    call CFG_add_get(cfg, "fluid%limit_velocity", fluid_limit_velocity, &
         "If true, keep the velocity constant for E > table_max_efield")

    ! Exit here if the fluid model is not used
    if (model_type /= model_fluid) return

    select case (integrator)
    case ("euler")
       time_step_method = forward_euler
    case ("trapezoidal")
       time_step_method = trapezoidal
    case ("rk2")
       time_step_method = rk2
    case ("rk4")
       time_step_method = rk4
    case default
       print *, "Unknown time integrator: ", trim(integrator)
       error stop
    end select

    n = n_ghost_cells
    allocate(fluid_state%a(1-n:domain_nx+n, num_arrays))
    allocate(fluid_state%s(num_scalars))

    fluid_state%a(:, :) = 0.0_dp
    fluid_state%s(:)    = 0.0_dp

    ! Create a lookup table for the model coefficients
    fluid_lkp_fld = LT_create(0.0_dp, fluid_max_field, table_size, 0)

    call CFG_get(cfg, "fluid%fld_mob", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_add_col(fluid_lkp_fld, x_data, y_data)
    if_mob = fluid_lkp_fld%n_cols

    call CFG_get(cfg, "fluid%fld_dif", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_add_col(fluid_lkp_fld, x_data, y_data)
    if_dif = fluid_lkp_fld%n_cols

    call CFG_get(cfg, "fluid%fld_alpha", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_add_col(fluid_lkp_fld, x_data, y_data)
    if_src = fluid_lkp_fld%n_cols
    call CFG_get(cfg, "fluid%fld_eta", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_add_col(fluid_lkp_fld, x_data, y_data)
    if_att = fluid_lkp_fld%n_cols

    fluid_num_fld_coef = fluid_lkp_fld%n_cols

    ! Determine extrapolation coefficients for transport data
    x = [fluid_max_field, 0.8_dp * fluid_max_field]

    ! Linear extrapolation for ionization coefficient
    y = LT_get_col(fluid_lkp_fld, if_src, x)

    extrap_src_dydx = (y(2) - y(1)) / (x(2) - x(1))
    extrap_src_y0   = y(2)

    ! Exponential decay for mobility: mu = c0 * exp(-c1 * E)
    y = LT_get_col(fluid_lkp_fld, if_mob, x)

    extrap_mob_c1 = -log(y(2)/y(1)) / (x(2) - x(1))
    extrap_mob_c0 = y(2) / exp(-extrap_mob_c1 * x(2))

    ! Initialization of electron and ion density
    do n = 1, domain_nx
       xx = (n-1) * domain_dx
       fluid_state%a(n, iv_elec) = init_elec_dens(xx)
       fluid_state%a(n, iv_pion) = init_pos_ion_dens(xx)
       fluid_state%a(n, iv_nion) = 0.0_dp
    end do

  end subroutine fluid_initialize

  subroutine set_boundary_conditions(state)
    type(state_t), intent(inout) :: state

    if (dielectric_present(1)) then
       ! Inside dielectric, density is zero
       state%a(0, iv_elec) = 0.0_dp
       state%a(-1, iv_elec) = 0.0_dp
    else
       ! Neumann boundary condition on the left
       state%a(0, iv_elec) = state%a(1, iv_elec)
       state%a(-1, iv_elec) = state%a(2, iv_elec)
    end if

    if (dielectric_present(2)) then
       state%a(domain_nx+1, iv_elec) = 0.0_dp
       state%a(domain_nx+2, iv_elec) = 0.0_dp
    else
       ! Neumann boundary condition on the right
       state%a(domain_nx+1, iv_elec) = state%a(domain_nx, iv_elec)
       state%a(domain_nx+2, iv_elec) = state%a(domain_nx-1, iv_elec)
    end if

  end subroutine set_boundary_conditions

  !> Compute the time derivatives of the fluid variables
  subroutine fluid_derivs(state, time, derivs, dt_max)
    use m_units_constants
    type(state_t), intent(inout)    :: state  !< Current state
    real(dp), intent(in)            :: time   !< Current time
    type(state_t), intent(inout)    :: derivs !< Derivatives
    !> If present, output the maximal time step
    real(dp), intent(out), optional :: dt_max

    type(LT_loc_t), allocatable :: fld_locs(:)
    integer                     :: n, nx
    real(dp)                    :: surface_charge(2), se
    real(dp), allocatable       :: mob_e(:), diff_e(:), src_e(:)
    real(dp), allocatable       :: mob_i(:), diff_i(:)
    real(dp), allocatable       :: flux(:)
    real(dp), allocatable       :: source(:)
    real(dp), allocatable       :: sigma(:)
    real(dp)                    :: dt_cfl, dt_dif, dt_drt
    real(dp), parameter         :: eps = 1e-100_dp

    nx = domain_nx

    derivs%a(:, :) = 0.0_dp
    derivs%s(:)    = 0.0_dp

    call set_boundary_conditions(state)

    allocate(flux(1:nx+1))
    allocate(mob_e(nx+1))
    allocate(mob_i(nx+1))
    allocate(diff_e(nx+1))
    allocate(diff_i(nx+1))
    allocate(src_e(nx+1))
    allocate(fld_locs(nx+1))

    ! Get electric field
    source = UC_elem_charge * (state%a(1:nx, iv_pion) - state%a(1:nx, iv_elec) &
         - state%a(1:nx, iv_nion))
    surface_charge = UC_elem_charge * (&
         state%s([i_lbound_pion, i_rbound_pion]) - &
         state%s([i_lbound_elec, i_rbound_elec]))
    call compute_field(source, surface_charge, time)

    ! Get locations in the lookup table
    fld_locs(:) = LT_get_loc(fluid_lkp_fld, abs(field_fc))

    ! Get coefficients
    mob_e  = LT_get_col_at_loc(fluid_lkp_fld, if_mob, fld_locs)
    diff_e = LT_get_col_at_loc(fluid_lkp_fld, if_dif, fld_locs)
    src_e  = LT_get_col_at_loc(fluid_lkp_fld, if_src, fld_locs)

    ! Extrapolate the ionization coefficient, assuming the energy per ionization
    ! remains constant
    where (abs(field_fc) > fluid_max_field)
       src_e = extrap_src(abs(field_fc))
       mob_e = extrap_mob(abs(field_fc))
    end where

    ! Compute electron flux
    call get_flux_1d(state%a(:, iv_elec), -mob_e * field_fc, diff_e, &
         domain_dx, flux, nx, n_ghost_cells)

    ! Compute source term per cell using the fluxes at cell faces
    call cell_face_to_center(source, src_e * abs(flux))
    derivs%a(1:nx, iv_elec) = source
    derivs%a(1:nx, iv_pion) = source

    do n = 1, nx
       derivs%a(n, iv_elec) = derivs%a(n, iv_elec) + &
            domain_inv_dx * (flux(n) - flux(n+1))
    end do

    ! Take into account boundary fluxes
    derivs%s(i_lbound_elec) = -flux(1)
    derivs%s(i_rbound_elec) = flux(nx+1)

    ! Photon fluxes on left and right sides. Assume photons are not absorbed and
    ! have a 50% chance of going left/right.
    if (se_photons_per_ionization > 0.0_dp) then
       se = 0.5_dp * sum(source) * domain_dx * se_photons_per_ionization

       ! Electrons will only be released if the field points towards the surface
       if (field_fc(1) < 0.0_dp) then
          derivs%a(1, iv_elec) = derivs%a(1, iv_elec) + se * domain_inv_dx
          derivs%s(i_lbound_elec) = derivs%s(i_lbound_elec) - se
       end if

       if (field_fc(nx+1) > 0.0_dp) then
          derivs%a(nx, iv_elec) = derivs%a(nx, iv_elec) + se * domain_inv_dx
          derivs%s(i_rbound_elec) = derivs%s(i_rbound_elec) - se
       end if
    end if

    ! Ion transport
    if (ion_transport) then
       diff_i(:) = ion_diffusion
       mob_i(:) = ion_mobility

       ! Compute ion flux
       call get_flux_1d(state%a(:, iv_pion), mob_i * field_fc, &
            diff_i, domain_dx, flux, nx, n_ghost_cells)
       do n = 1, nx
          derivs%a(n, iv_pion) = derivs%a(n, iv_pion) + &
               domain_inv_dx * (flux(n) - flux(n+1))
       end do

       ! Take into account boundary fluxes
       derivs%s(i_lbound_pion) = -flux(1)
       derivs%s(i_rbound_pion) = flux(nx+1)

       ! Secondary emission of electrons
       se = -ion_secondary_emission_yield * flux(1)
       derivs%a(1, iv_elec) = derivs%a(1, iv_elec) + se * domain_inv_dx
       derivs%s(i_lbound_elec) = derivs%s(i_lbound_elec) - se

       se = ion_secondary_emission_yield * flux(nx+1)
       derivs%a(nx, iv_elec) = derivs%a(nx, iv_elec) + se * domain_inv_dx
       derivs%s(i_rbound_elec) = derivs%s(i_rbound_elec) - se
    else
       diff_i(:) = 0.0_dp
       mob_i(:)  = 0.0_dp
    end if

    if (present(dt_max)) then
       ! Determine maximal time step

       ! CFL condition
       dt_cfl = domain_dx / max(eps, &
            maxval(abs(field_fc * max(mob_e, mob_i))))

       ! Diffusion condition
       dt_dif = 0.5_dp * domain_dx**2 / &
            max(eps, maxval(diff_e), maxval(diff_i))

       ! Determine conductivity (overestimate by taking maximum, and assume
       ! quasi-neutrality for ion contribution)
       allocate(sigma(nx+1))
       sigma(:) = (mob_e + mob_i) * &
            max(state%a(0:nx, iv_elec), state%a(1:nx+1, iv_elec))

       ! Dielectric relaxation time
       dt_drt = UC_eps0 / (UC_elem_charge * max(eps, maxval(sigma)))

       ! Take the minimum of the CFL condition with Courant number 0.5 and
       ! the combined CFL-diffusion condition with Courant number 1.0. The
       ! 0.5 is emperical, to have good accuracy (and TVD/positivity) in
       ! combination with the explicit trapezoidal rule
       dt_max = min(0.5_dp * dt_cfl, 1/(1/dt_cfl + 1/dt_dif))

       ! Use a 'safety' factor of 0.9
       dt_max = 0.9_dp * min(dt_max, dt_drt)
    end if
  end subroutine fluid_derivs

  !> Average source terms defined at cell faces to get source terms at cell
  !> centers
  subroutine cell_face_to_center(dens_c, src_f)
    real(dp), intent(inout) :: dens_c(:)
    real(dp), intent(in)    :: src_f(:)
    integer                 :: n

    do n = 1, size(dens_c)
       dens_c(n) = 0.5_dp * (src_f(n) + src_f(n+1))
    end do
  end subroutine cell_face_to_center

  !> Write output files
  subroutine fluid_write_output(base_fname, time, dt, ix)
    use m_units_constants
    character(len=*), intent(in) :: base_fname
    real(dp), intent(in)         :: time, dt
    integer, intent(in)          :: ix
    integer                      :: n, nx
    character(len=200)           :: fname
    real(dp)                     :: total_charge
    real(dp)                     :: max_field, deriv_max_field
    real(dp), save               :: prev_max_field, prev_time
    integer                      :: my_unit

    nx = domain_nx

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_fluid_", ix, ".txt"
    open(newunit=my_unit, file=trim(fname))

    write(my_unit, "(A)") "x field electron pos_ion potential"
    do n = 1, domain_nx
       write(my_unit, *) domain_dx * (n-0.5_dp), &
            field_cc(n), &
            fluid_state%a(n, iv_elec), &
            fluid_state%a(n, iv_pion), &
            potential(n)
    end do
    close(my_unit)

    write(*, "(A,E9.2,A,A)") " t = ", time, " wrote ", trim(fname)

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_fluid_scalars.txt"
    if (ix == 1) then
       open(newunit=my_unit, file=trim(fname))
       write(my_unit, *) "time dt E_max total_charge sum_elec ", &
            "sum_pion sigma_l sigma_r max_elec max_pion deriv_Emax"
    else
       open(newunit=my_unit, file=trim(fname), access='append')
    end if

    total_charge = domain_dx * sum(fluid_state%a(1:nx, iv_pion) &
         - fluid_state%a(1:nx, iv_elec)) &
         - fluid_state%s(i_lbound_elec) - fluid_state%s(i_rbound_elec) &
         + fluid_state%s(i_lbound_pion) + fluid_state%s(i_rbound_pion)

    max_field = maxval(abs(field_fc))

    if (ix == 1) then
       deriv_max_field = 0.0_dp
    else
       deriv_max_field = (max_field - prev_max_field) / (time - prev_time)
    end if

    write(my_unit, *) time, dt, max_field, total_charge, &
         domain_dx * sum(fluid_state%a(1:nx, iv_elec)), &
         domain_dx * sum(fluid_state%a(1:nx, iv_pion)), &
         fluid_state%s(i_lbound_pion) - fluid_state%s(i_lbound_elec), &
         fluid_state%s(i_rbound_pion) - fluid_state%s(i_rbound_elec), &
         maxval(fluid_state%a(1:nx, iv_elec)), &
         maxval(fluid_state%a(1:nx, iv_pion)), deriv_max_field
    close(my_unit)

    prev_max_field = max_field
    prev_time      = time
  end subroutine fluid_write_output

  subroutine fluid_advance(dt, time, max_dt)
    real(dp), intent(in)  :: dt
    real(dp), intent(in)  :: time
    real(dp), intent(out) :: max_dt
    real(dp), parameter   :: one_sixth = 1 / 6.0_dp

    type(state_t) :: derivs
    type(state_t) :: sum_derivs
    type(state_t) :: substep

    ! Allocate by assignment
    derivs = fluid_state

    select case (time_step_method)
    case (forward_euler)
       call fluid_derivs(fluid_state, time, derivs, max_dt)
       fluid_state%a = fluid_state%a + dt * derivs%a
       fluid_state%s = fluid_state%s + dt * derivs%s
    case (trapezoidal)
       substep = fluid_state

       call fluid_derivs(substep, time, derivs)
       substep%a = substep%a + dt * derivs%a
       substep%s = substep%s + dt * derivs%s

       call fluid_derivs(substep, time+dt, derivs, max_dt)
       substep%a = substep%a + dt * derivs%a
       substep%s = substep%s + dt * derivs%s

       fluid_state%a = 0.5_dp * (fluid_state%a + substep%a)
       fluid_state%s = 0.5_dp * (fluid_state%s + substep%s)
    case (rk2)
       substep = fluid_state

       ! Step 1 (at initial time)
       call fluid_derivs(substep, time, derivs)
       substep%a = substep%a + 0.5_dp * dt * derivs%a
       substep%s = substep%s + 0.5_dp * dt * derivs%s

       ! Step 2 (at initial time + dt/2)
       call fluid_derivs(substep, time + 0.5_dp * dt, derivs, max_dt)

       ! Result (at initial time + dt)
       fluid_state%a = fluid_state%a + dt * derivs%a
       fluid_state%s = fluid_state%s + dt * derivs%s
    case (rk4)
       substep = fluid_state

       ! Step 1 (at initial time)
       call fluid_derivs(substep, time, derivs)
       sum_derivs = derivs

       ! Step 2 (at initial time + dt/2)
       substep%a = substep%a + 0.5_dp * dt * derivs%a
       substep%s = substep%s + 0.5_dp * dt * derivs%s
       call fluid_derivs(substep, time + 0.5_dp * dt, derivs)

       sum_derivs%a = sum_derivs%a + 2 * derivs%a
       sum_derivs%s = sum_derivs%s + 2 * derivs%s

       ! Step 3 (at initial time + dt/2)
       substep%a = fluid_state%a + 0.5_dp * dt * derivs%a
       substep%s = fluid_state%s + 0.5_dp * dt * derivs%s
       call fluid_derivs(substep, time + 0.5_dp * dt, derivs)

       sum_derivs%a = sum_derivs%a + 2 * derivs%a
       sum_derivs%s = sum_derivs%s + 2 * derivs%s

       ! Step 4 (at initial time + dt)
       substep%a = fluid_state%a + dt * derivs%a
       substep%s = fluid_state%s + dt * derivs%s
       call fluid_derivs(substep, time + dt, derivs, max_dt)

       sum_derivs%a = sum_derivs%a + derivs%a
       sum_derivs%s = sum_derivs%s + derivs%s

       ! Combine time derivatives at steps
       fluid_state%a = fluid_state%a + dt * one_sixth * sum_derivs%a
       fluid_state%s = fluid_state%s + dt * one_sixth * sum_derivs%s
    case default
       error stop "Unknown time stepping scheme"
    end select

    where (fluid_state%a < fluid_small_density)
       fluid_state%a = 0.0_dp
    end where
  end subroutine fluid_advance

  !> Extrapolate the ionization coefficient
  elemental real(dp) function extrap_src(fld)
    real(dp), intent(in) :: fld
    ! Linear extrapolation
    extrap_src = extrap_src_y0 + extrap_src_dydx * (fld - fluid_max_field)
  end function extrap_src

  !> Extrapolate the mobility
  elemental real(dp) function extrap_mob(fld)
    real(dp), intent(in) :: fld

    if (fluid_limit_velocity) then
       ! Limit the maximal velocity (for testing purposes)
       extrap_mob = extrap_mob_c0 * exp(-extrap_mob_c1 * fluid_max_field) * fluid_max_field/fld
    else
       ! Exponential decay
       extrap_mob = extrap_mob_c0 * exp(-extrap_mob_c1 * fld)
    end if

  end function extrap_mob

end module m_fluid_1d
