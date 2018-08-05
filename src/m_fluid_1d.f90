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

  integer            :: num_arrays
  integer            :: iv_elec, iv_pion, iv_nion

  integer, parameter :: num_scalars = 4
  integer, parameter :: i_lbound_elec = 1
  integer, parameter :: i_rbound_elec = 2
  integer, parameter :: i_lbound_pion = 3
  integer, parameter :: i_rbound_pion = 4

  integer :: if_mob, if_dif, if_src
  integer :: if_att, fluid_num_fld_coef

  type state_t
     real(dp), allocatable :: a(:, :) ! arrays
     real(dp), allocatable :: s(:)    ! scalars
  end type state_t

  type(state_t) :: fluid_state

  type(LT_t) :: fluid_lkp_fld

  !> Whether ion transport is included
  logical :: ion_transport = .true.

  !> Mobility for ions
  real(dp) :: ion_mobility = 5e-4_dp

  !> Diffusion coefficient for ions
  real(dp) :: ion_diffusion = 1e-4_dp

  !> Secondary emission yield for ion impact
  real(dp) :: secondary_emission_yield = 1.0e-2_dp

  real(dp) :: photons_per_ionization = 1e-2_dp

  !> Secondary emission yield for photons
  real(dp) :: photon_emission_yield = 1.0e-5_dp

  !> Set density to zero below this value
  real(dp) :: fluid_small_density = 1.0_dp

  public :: fluid_initialize
  public :: fluid_advance
  public :: fluid_write_output

contains

  subroutine fluid_initialize(cfg)
    use m_transport_data
    use m_config

    type(CFG_t), intent(inout) :: cfg
    integer, parameter         :: name_len = 100
    character(len=200)         :: input_file
    integer                    :: n, table_size
    real(dp)                   :: max_efield, xx
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=100)         :: data_name
    character(len=40)          :: integrator

    input_file = "input/n2_transport_data_siglo.txt"
    call CFG_add_get(cfg, "fluid%input_file", input_file, &
         "Input file with cross sections")

    table_size = 1000
    call CFG_add_get(cfg, "fluid%table_size", table_size, &
         "Size of lookup table for transport coefficients")

    max_efield = 2.5e7
    call CFG_add_get(cfg, "fluid%table_max_efield", max_efield, &
         "Maximum electric field in transport coefficient table")

    integrator = "rk2"
    call CFG_add_get(cfg, "fluid%integrator", integrator, &
         "Time integrator (euler, trapezoidal, rk2, rk4)")

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

    num_arrays = 3
    iv_elec = 1
    iv_pion = 2
    iv_nion = 3

    n = n_ghost_cells
    allocate(fluid_state%a(1-n:domain_nx+n, num_arrays))
    allocate(fluid_state%s(num_scalars))

    fluid_state%a(:, :) = 0.0_dp
    fluid_state%s(:)    = 0.0_dp

    ! Create a lookup table for the model coefficients
    fluid_lkp_fld = LT_create(0.0_dp, max_efield, table_size, 0)

    call CFG_add(cfg, "fluid%fld_mob", "efield[V/m]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    call CFG_add(cfg, "fluid%fld_dif", "efield[V/m]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(cfg, "fluid%fld_alpha", "efield[V/m]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(cfg, "fluid%fld_eta", "efield[V/m]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")
    ! call CFG_add(cfg, "fluid%fld_det", "efield[V/m]_vs_det[1/s]", &
    !      "The name of the detachment rate coeff.")

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

    ! if (fluid_use_detach) then
    !    call CFG_get(cfg, "fluid_fld_det", data_name)
    !    call TD_get_td_from_file(input_file, gas_name, &
    !         trim(data_name), x_data, y_data)
    !    call LT_add_col(fluid_lkp_fld, x_data, y_data)
    !    if_det= fluid_lkp_fld%n_cols
    ! end if

    fluid_num_fld_coef = fluid_lkp_fld%n_cols

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

  subroutine fluid_derivs(state, time, derivs, dt_max)
    use m_units_constants
    type(state_t), intent(inout)    :: state
    real(dp), intent(in)            :: time
    type(state_t), intent(inout)    :: derivs
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
    mob_e = LT_get_col_at_loc(fluid_lkp_fld, if_mob, fld_locs)
    diff_e = LT_get_col_at_loc(fluid_lkp_fld, if_dif, fld_locs)
    src_e = LT_get_col_at_loc(fluid_lkp_fld, if_src, fld_locs)

    ! Compute electron flux
    call get_flux_1d(state%a(:, iv_elec), -mob_e * field_fc, diff_e, &
         domain_dx, flux, nx, n_ghost_cells)

    ! Compute source term per cell using the fluxes at cell faces
    call set_stagg_source_1d(source, src_e * abs(flux))
    derivs%a(1:nx, iv_elec) = source
    derivs%a(1:nx, iv_pion) = source

    ! ~~~ Attachment source ~~~ TODO: try different formulations

    ! Detachment process
    ! if (fluid_use_detach) then
    !    source(1) = 0
    !    source(n_cc) = 0
    !    source(2:n_cc-1) = (det_c(1:n_cc-2) + det_c(2:n_cc-1)) * 0.5_dp * state%a(2:n_cc-1, iv_nion)
    !    time_derivs(:, iv_elec) = time_derivs(:, iv_elec) + source
    !    time_derivs(:, iv_nion) = time_derivs(:, iv_nion) - source
    ! end if

    do n = 1, nx
       derivs%a(n, iv_elec) = derivs%a(n, iv_elec) + &
            domain_inv_dx * (flux(n) - flux(n+1))
    end do

    ! Take into account boundary fluxes
    derivs%s(i_lbound_elec) = -flux(1)
    derivs%s(i_rbound_elec) = flux(nx+1)

    ! Photon fluxes on left and right sides
    ! se = 0.5_dp * sum(source) * &
    !      photons_per_ionization * photon_emission_yield

    ! derivs%a(1, iv_elec) = derivs%a(1, iv_elec) + se * domain_inv_dx
    ! derivs%a(nx, iv_elec) = derivs%a(nx, iv_elec) + se * domain_inv_dx
    ! derivs%s(i_lbound_elec) = derivs%s(i_lbound_elec) - se
    ! derivs%s(i_rbound_elec) = derivs%s(i_lbound_elec) - se

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
       se = -secondary_emission_yield * flux(1)
       derivs%a(1, iv_elec) = derivs%a(1, iv_elec) + se * domain_inv_dx
       derivs%s(i_lbound_elec) = derivs%s(i_lbound_elec) - se

       se = secondary_emission_yield * flux(nx+1)
       derivs%a(nx, iv_elec) = derivs%a(nx, iv_elec) + se * domain_inv_dx
       derivs%s(i_rbound_elec) = derivs%s(i_rbound_elec) - se
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
       dt_max = min(dt_max, dt_drt)
    end if
  end subroutine fluid_derivs

  subroutine set_stagg_source_1d(dens_c, src_f)
    real(dp), intent(inout) :: dens_c(:)
    real(dp), intent(in)    :: src_f(:)
    integer                 :: n

    do n = 1, size(dens_c)
       dens_c(n) = 0.5_dp * (src_f(n) + src_f(n+1))
    end do
  end subroutine set_stagg_source_1d

  subroutine fluid_write_output(base_fname, time, ix)
    use m_units_constants
    character(len=*), intent(in) :: base_fname
    real(dp), intent(in)         :: time
    integer, intent(in)          :: ix
    integer                      :: n, nx
    character(len=200)           :: fname
    real(dp)                     :: total_charge
    integer                      :: my_unit

    nx = domain_nx

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_fluid_", ix, ".txt"
    open(newunit=my_unit, file=trim(fname))

    write(my_unit, "(A)") "x field electron pos_ion"
    do n = 1, domain_nx
       write(my_unit, "(4e13.5)") domain_dx * (n-0.5_dp), &
            field_cc(n), &
            fluid_state%a(n, iv_elec), &
            fluid_state%a(n, iv_pion)
    end do
    close(my_unit)

    write(*, "(A,E9.2,A,A)") " t = ", time, " wrote ", trim(fname)

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_fluid_scalars.txt"
    if (ix == 1) then
       open(newunit=my_unit, file=trim(fname))
       write(my_unit, *) "# Header TODO"
    else
       open(newunit=my_unit, file=trim(fname), access='append')
    end if

    total_charge = domain_dx * sum(fluid_state%a(1:nx, iv_pion) &
         - fluid_state%a(1:nx, iv_elec)) &
         - fluid_state%s(i_lbound_elec) - fluid_state%s(i_rbound_elec) &
         + fluid_state%s(i_lbound_pion) + fluid_state%s(i_rbound_pion)

    write(my_unit, *) time, fluid_state%s(i_lbound_elec), &
         fluid_state%s(i_rbound_elec), total_charge
    close(my_unit)
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

end module m_fluid_1d
