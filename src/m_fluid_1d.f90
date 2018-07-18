!> Implementation of a drift-diffusion-reaction plasma fluid model
module m_fluid_1d
  use m_generic
  use m_lookup_table

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter :: n_ghost_cells = 2

  integer, parameter :: forward_euler    = 1
  integer, parameter :: rk2              = 2
  integer, parameter :: rk4              = 4
  integer            :: time_step_method = forward_euler

  integer            :: num_arrays
  integer            :: iv_elec, iv_pion, iv_en, iv_nion

  integer, parameter :: num_scalars = 1
  integer, parameter :: i_lbound_elec = 1

  integer :: if_mob, if_dif, if_src, if_en
  integer :: if_att, if_det, fluid_num_fld_coef
  integer :: fluid_ie_mob, fluid_ie_dif, fluid_ie_src
  integer :: fluid_ie_att, fluid_ie_loss, fluid_num_en_coef
  ! logical :: fluid_use_en, fluid_use_en_mob, fluid_use_en_dif
  ! logical :: fluid_use_en_src, fluid_use_detach

  type state_t
     real(dp), allocatable :: a(:, :) ! arrays
     real(dp), allocatable :: s(:)    ! scalars
  end type state_t

  type(state_t) :: fluid_state

  real(dp) :: fluid_small_dens
  type(LT_t) :: fluid_lkp_fld, fluid_lkp_en

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
    real(dp)                   :: max_energy, max_efield, xx
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=100)         :: data_name

    input_file = "input/n2_transport_data_siglo.txt"
    call CFG_add_get(cfg, "fluid%input_file", input_file, &
         "Input file with cross sections")

    table_size = 1000
    call CFG_add_get(cfg, "fluid%table_size", table_size, &
         "Size of lookup table for transport coefficients")

    ! call CFG_get(cfg, "fluid_lkptbl_max_energy", max_energy)

    max_efield = 2.5e7
    call CFG_add_get(cfg, "fluid%table_max_efield", max_efield, &
         "Maximum electric field in transport coefficient table")
    ! call CFG_get(cfg, "fluid%small_density", fluid_small_dens, &
    !      "Small fluid density")
    ! fluid_max_energy = max_energy

    ! if (MODEL_type == MODEL_fluid_ee) then
    !    fluid_use_en     = .true.
    !    call CFG_get(cfg, "fluid_use_en_mob", fluid_use_en_mob)
    !    call CFG_get(cfg, "fluid_use_en_dif", fluid_use_en_dif)
    !    call CFG_get(cfg, "fluid_use_en_src", fluid_use_en_src)
    !    call CFG_get(cfg, "fluid_use_detach", fluid_use_detach)
    ! else
    !    fluid_use_en     = .false.
    !    fluid_use_en_mob = .false.
    !    fluid_use_en_dif = .false.
    !    fluid_use_en_src = .false.
    !    fluid_use_detach = .false.
    ! end if

    num_arrays = 3
    iv_elec = 1
    iv_pion = 2
    iv_nion = 3

    ! if (fluid_use_en) then
    !    num_arrays = num_arrays + 1
    !    iv_en = num_arrays
    ! end if

    n = n_ghost_cells
    allocate(fluid_state%a(1-n:domain_ncell+n, num_arrays))
    allocate(fluid_state%s(num_scalars))

    fluid_state%a(:, :) = 0.0_dp
    fluid_state%s(:)    = 0.0_dp

    ! Create a lookup table for the model coefficients
    fluid_lkp_fld = LT_create(0.0_dp, max_efield, table_size, 0)

    call CFG_add(cfg, "fluid%fld_mob", "efield[V/m]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    ! call CFG_add(cfg, "fluid%fld_en", "efield[V/m]_vs_energy[eV]", &
    !      "The name of the energy(fld) coefficient")
    call CFG_add(cfg, "fluid%fld_dif", "efield[V/m]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(cfg, "fluid%fld_alpha", "efield[V/m]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(cfg, "fluid%fld_eta", "efield[V/m]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")
    ! call CFG_add(cfg, "fluid%fld_loss", "efield[V/m]_vs_loss[eV/s]", &
    !      "The name of the energy loss coeff.")
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

    ! call CFG_get(cfg, "fluid_fld_en", data_name)
    ! call TD_get_td_from_file(input_file, gas_name, &
    !      trim(data_name), x_data, y_data)
    ! call LT_add_col(fluid_lkp_fld, x_data, y_data)
    ! if_en = fluid_lkp_fld%n_cols

    fluid_num_fld_coef = fluid_lkp_fld%n_cols
    fluid_num_en_coef  = fluid_lkp_en%n_cols

    ! Initialization of electron and ion density
    do n = 1, domain_ncell
       xx = (n-1) * domain_dx
       fluid_state%a(n, iv_elec) = init_elec_dens(xx)
       fluid_state%a(n, iv_pion) = init_pos_ion_dens(xx)
       fluid_state%a(n, iv_nion) = 0.0_dp
    end do

  end subroutine fluid_initialize

  subroutine set_boundary_conditions(state)
    type(state_t), intent(inout) :: state

    ! Neumann boundary condition on the right
    state%a(domain_ncell+1, iv_elec) = state%a(domain_ncell, iv_elec)
    state%a(domain_ncell+2, iv_elec) = state%a(domain_ncell-1, iv_elec)

    if (.not. left_dielectric_present) then
       ! Neumann boundary condition on the left
       state%a(0, iv_elec) = state%a(1, iv_elec)
       state%a(-1, iv_elec) = state%a(2, iv_elec)
    end if
  end subroutine set_boundary_conditions

  subroutine fluid_derivs(state, time, derivs)
    use m_units_constants
    type(state_t), intent(inout) :: state
    real(dp), intent(in)         :: time
    type(state_t), intent(inout) :: derivs

    integer                     :: n, nx
    real(dp), parameter         :: five_third = 5 / 3.0_dp
    real(dp), allocatable       :: mob_c(:), dif_c(:), src_c(:)
    real(dp), allocatable       :: flux(:)
    real(dp), allocatable       :: source(:)
    type(LT_loc_t), allocatable :: fld_locs(:)

    nx = domain_ncell

    call set_boundary_conditions(state)

    allocate(flux(1:nx+1))
    ! allocate(mob_c(domain_ncell))
    ! allocate(dif_c(domain_ncell))
    ! allocate(src_c(domain_ncell))
    ! allocate(att_c(domain_ncell))
    ! allocate(det_c(domain_ncell))

    ! Get electric field
    source = UC_elem_charge * (state%a(1:nx, iv_pion) - state%a(1:nx, iv_elec) &
         - state%a(1:nx, iv_nion))
    call compute_field(source, -UC_elem_charge * state%s(i_lbound_elec), time)

    ! Get locations in the lookup table
    fld_locs = LT_get_loc(fluid_lkp_fld, abs(field_fc))
    ! if (fluid_use_en) then
    !    ! There is a regularization: at density zero, we use the energy corresponding to the efield
    !    fld_en = LT_get_col_at_loc(fluid_lkp_fld, if_en, fld_locs)
    !    mean_en = (state%a(1:n_cc-1, iv_en) + state%a(2:n_cc, iv_en) + 2 * fluid_small_dens * fld_en) &
    !         / (state%a(1:n_cc-1, iv_elec) + state%a(2:n_cc, iv_elec) + 2 * fluid_small_dens)
    !    en_locs = LT_get_loc(fluid_lkp_en, mean_en)
    !    en_loss = LT_get_col_at_loc(fluid_lkp_en, fluid_ie_loss, en_locs)
    ! end if

    ! if (fluid_use_en_mob) then
    !    mob_c = LT_get_col_at_loc(fluid_lkp_en, fluid_ie_mob, en_locs)
    ! else
    mob_c = LT_get_col_at_loc(fluid_lkp_fld, if_mob, fld_locs)
    dif_c = LT_get_col_at_loc(fluid_lkp_fld, if_dif, fld_locs)
    src_c = LT_get_col_at_loc(fluid_lkp_fld, if_src, fld_locs)
    ! end if

    ! if (fluid_use_en_dif) then
    ! dif_c = LT_get_col_at_loc(fluid_lkp_en, fluid_ie_dif, en_locs)
    ! else

    ! end if

    ! if (fluid_use_en_src) then
    !    src_c = LT_get_col_at_loc(fluid_lkp_en, fluid_ie_src, en_locs)
    !    att_c = LT_get_col_at_loc(fluid_lkp_en, fluid_ie_att, en_locs)
    ! else
    ! att_c = LT_get_col_at_loc(fluid_lkp_fld, if_att, fld_locs)
    ! end if
    ! if (fluid_use_detach) det_c = LT_get_col_at_loc(fluid_lkp_fld, if_det, fld_locs)

    ! ~~~ Electron transport ~~~
    call get_flux_1d(state%a(:, iv_elec), -mob_c * field_fc, dif_c, &
         domain_dx, flux, domain_ncell, n_ghost_cells)

    ! ~~~ Ionization source ~~~ TODO: try different formulations
    call set_stagg_source_1d(source, src_c * abs(flux))
    derivs%a(1:nx, iv_elec) = source
    derivs%a(1:nx, iv_pion) = source

    ! ~~~ Attachment source ~~~ TODO: try different formulations
    ! call set_stagg_source_1d(source, att_c * abs(flux))
    ! time_derivs(:, iv_elec) = time_derivs(:, iv_elec) - source
    ! time_derivs(:, iv_nion) = time_derivs(:, iv_nion) + source

    ! Detachment process
    ! if (fluid_use_detach) then
    !    source(1) = 0
    !    source(n_cc) = 0
    !    source(2:n_cc-1) = (det_c(1:n_cc-2) + det_c(2:n_cc-1)) * 0.5_dp * state%a(2:n_cc-1, iv_nion)
    !    time_derivs(:, iv_elec) = time_derivs(:, iv_elec) + source
    !    time_derivs(:, iv_nion) = time_derivs(:, iv_nion) - source
    ! end if

    do n = 1, domain_ncell
       derivs%a(n, iv_elec) = derivs%a(n, iv_elec) + &
            domain_inv_dx * (flux(n) - flux(n+1))
    end do

    if (left_dielectric_present) then
       ! Store changes in interior cell as surface charge
       derivs%s(i_lbound_elec) = domain_dx * &
            derivs%a(left_dielectric_iface-1, iv_elec)
       derivs%a(left_dielectric_iface-1, iv_elec) = 0.0_dp
    end if

    ! if (fluid_use_en) then
    !    ! ~~~ Energy source ~~~
    !    source(n_cc)     = 0.0_dp
    !    source(1:n_cc-1) = -0.5_dp * (fld * flux + en_loss * state%a(1:n_cc-1, iv_elec))
    !    source(2:n_cc)   = source(2:n_cc) - 0.5_dp * (fld * flux + en_loss * &
    !         state%a(2:n_cc, iv_elec))

    !    ! ~~~ energy transport ~~~
    !    call fluid_transport_scheme(state%a(:, iv_en), -five_third * mob_c * fld, &
    !         five_third * dif_c * inv_delta_x, flux)

    !    ! Set time derivatives
    !    time_derivs(:, iv_en) = source
    !    call add_grad_flux_1d(time_derivs(:, iv_en), flux * inv_delta_x)
    ! end if

    ! Dirichlet
    ! time_derivs(1, :) = 0.0_dp
    ! time_derivs(n_cc, :) = 0.0_dp
    ! Neumann
    ! time_derivs(1, :) = time_derivs(2, :)
    ! time_derivs(n_cc, :) = time_derivs(n_cc-1, :)
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
    real(dp)                     :: total_charge, elec_free, elec_bound
    integer                      :: my_unit

    nx = domain_ncell
    total_charge = domain_dx * sum(fluid_state%a(1:nx, iv_pion) &
         - fluid_state%a(1:nx, iv_elec)) &
         - fluid_state%s(i_lbound_elec)
    ! elec_free = sum(fluid_state%a(1:nx, iv_elec)) * domain_dx
    ! elec_bound = fluid_state%s(i_lbound_elec)
    ! print *, "Q", elec_free, elec_bound, total_charge

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_fluid_", ix, ".txt"
    open(newunit=my_unit, file=trim(fname))

    do n = 1, domain_ncell
       write(my_unit, *) domain_dx * (n-0.5_dp), &
            fluid_state%a(n, iv_elec), &
            fluid_state%a(n, iv_pion), &
            field_cc(n)
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

    write(my_unit, *) time, fluid_state%s(i_lbound_elec), total_charge
    close(my_unit)
  end subroutine fluid_write_output

  ! subroutine fluid_get_coeffs(coeff_data, coeff_names, n_coeffs)
  !   real(dp), intent(out), allocatable :: coeff_data(:,:)
  !   character(len=*), intent(out), allocatable :: coeff_names(:)
  !   integer, intent(out) :: n_coeffs
  !   integer :: ix, n_rows, n_fld_coeffs, n_en_coeffs

  !   if (fluid_use_en) then
  !      n_en_coeffs = fluid_lkp_en%n_cols + 1
  !   else
  !      n_en_coeffs = 0
  !   end if
  !   n_fld_coeffs = fluid_lkp_fld%n_cols + 1

  !   n_coeffs = n_fld_coeffs + n_en_coeffs
  !   n_rows = fluid_lkp_fld%n_points
  !   allocate(coeff_data(n_coeffs, n_rows))
  !   allocate(coeff_names(n_coeffs))

  !   call LT_get_data(fluid_lkp_fld, coeff_data(1, :), coeff_data(2:n_fld_coeffs, :))
  !   coeff_names(1) = "efield (V/m)"
  !   coeff_names(1+if_en) = "fld_en (eV)"
  !   if (fluid_use_detach) coeff_names(1+if_det) = "fld_det (1/s)"
  !   if (.not. fluid_use_en_mob) coeff_names(1+if_mob) = "fld_mob (m2/Vs)"
  !   if (.not. fluid_use_en_dif) coeff_names(1+if_dif) = "fld_dif (m2/s)"
  !   if (.not. fluid_use_en_src) then
  !      coeff_names(1+if_src) = "fld_src (1/m)"
  !      coeff_names(1+if_att) = "fld_att (1/m)"
  !   end if

  !   if (fluid_use_en) then
  !      ix = n_fld_coeffs + 1
  !      call LT_get_data(fluid_lkp_en, coeff_data(ix, :), coeff_data(ix+1:, :))
  !      coeff_names(ix) = "energy (eV/s)"
  !      coeff_names(ix+fluid_ie_loss) = "en_loss (eV/s)"
  !      if (fluid_use_en_mob) coeff_names(ix+fluid_ie_mob) = "en_mob (m2/Vs)"
  !      if (fluid_use_en_dif) coeff_names(ix+fluid_ie_dif) = "en_dif (m2/s)"
  !      if (fluid_use_en_src) then
  !         coeff_names(ix+fluid_ie_src) = "en_src (1/m)"
  !         coeff_names(ix+fluid_ie_att) = "en_att (1/m)"
  !      end if
  !   end if

  ! end subroutine fluid_get_coeffs

  subroutine fluid_advance(time, dt)
    real(dp), intent(in) :: time
    real(dp), intent(in) :: dt
    real(dp), parameter  :: one_sixth = 1 / 6.0_dp

    type(state_t) :: derivs
    type(state_t) :: sum_derivs
    type(state_t) :: substep

    ! Allocate by assignment
    derivs = fluid_state

    select case (time_step_method)
    case (forward_euler)
       call fluid_derivs(fluid_state, time, derivs)
       fluid_state%a = fluid_state%a + dt * derivs%a
       fluid_state%s = fluid_state%s + dt * derivs%s
    case (rk2)
       ! Step 1 (at initial time)
       call fluid_derivs(fluid_state, time, derivs)
       substep = fluid_state
       substep%a = substep%a + 0.5_dp * dt * derivs%a
       substep%s = substep%s + 0.5_dp * dt * derivs%s

       ! Step 2 (at initial time + dt/2)
       call fluid_derivs(substep, time + 0.5_dp * dt, derivs)

       ! Result (at initial time + dt)
       fluid_state%a = fluid_state%a + 0.5_dp * dt * derivs%a
       fluid_state%s = fluid_state%s + 0.5_dp * dt * derivs%s
    case (rk4)
       ! Step 1 (at initial time)
       call fluid_derivs(fluid_state, time, derivs)
       sum_derivs = derivs

       ! Step 2 (at initial time + dt/2)
       substep = fluid_state
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
       call fluid_derivs(substep, time + dt, derivs)

       sum_derivs%a = sum_derivs%a + 2 * derivs%a
       sum_derivs%s = sum_derivs%s + 2 * derivs%s

       ! Combine time derivatives at steps
       fluid_state%a = fluid_state%a + dt * one_sixth * sum_derivs%a
       fluid_state%s = fluid_state%s + dt * one_sixth * sum_derivs%s
    case default
       error stop "Unknown time stepping scheme"
    end select
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
