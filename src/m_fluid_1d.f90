!> Implementation of a drift-diffusion-reaction plasma fluid model
module m_fluid_1d
  use m_generic
  use m_lookup_table
  use m_units_constants

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n_ghost_cells = 2
  real(dp), parameter :: huge_value = 1e100_dp

  integer, parameter :: forward_euler    = 1
  integer, parameter :: trapezoidal      = 2
  integer, parameter :: rk2              = 3
  integer, parameter :: rk4              = 4
  integer            :: time_step_method = rk2

  ! HTODO: Modify this to adjust when using LFA or LEA
  integer, parameter :: num_arrays = 7
  integer, parameter :: iv_elec    = 1
  integer, parameter :: iv_pion    = 2
  integer, parameter :: iv_nion    = 3
  integer, parameter :: iv_source    = 4
  integer, parameter :: iv_en = 5
  integer, parameter :: iv_flux_ad = 6
  integer, parameter :: iv_flux_dif = 7


  integer, parameter :: num_scalars = 4
  integer, parameter :: i_lbound_elec = 1
  integer, parameter :: i_rbound_elec = 2
  integer, parameter :: i_lbound_pion = 3
  integer, parameter :: i_rbound_pion = 4

  logical :: fluid_limit_velocity = .false.
  logical :: fluid_use_LEA = .false.

  ! Indices for the transport coefficients
  integer :: ix_mob             = -1   ! Mobility
  integer :: ix_diff            = -1  ! Diffusion
  integer :: ix_alpha           = -1 ! Ionization
  integer :: ix_eta             = -1   ! Attachment
  integer :: ix_energy          = -1   ! Electron mean Energy
  integer :: ix_loss            = -1   !Electron energy losses due to collisions
  integer :: fluid_num_fld_coef = -1
  integer :: fluid_num_en_coef  = -1

  type state_t
     real(dp), allocatable :: a(:, :) ! arrays
     real(dp), allocatable :: s(:)    ! scalars
  end type state_t

  type(state_t) :: fluid_state

  type(LT_t) :: fluid_lkp_fld, fluid_lkp_en

  !> Maximum electric field for the transport data table
  real(dp) :: fluid_max_field = 2.5e7_dp

  !> Maximum electron energy for the transport data table
  real(dp) :: fluid_max_en = 2.5e2_dp

  !> Disable diffusion parallel to fields above this threshold (V/m)
  real(dp) :: fluid_diffusion_emax = huge_value

  !> Set density to zero below this value
  real(dp) :: fluid_small_density = 1.0_dp

  !> Safety factor (< 1) for the time step
  real(dp) :: fluid_dt_factor = 0.9_dp


  ! Extrapolation coefficients for transport data
  real(dp) :: extrap_src_dydx
  real(dp) :: extrap_src_y0
  real(dp) :: extrap_mob_c0
  real(dp) :: extrap_mob_c1

  !> Limit fluxes to avoid the dielectric relaxation time
  logical :: fluid_avoid_drt = .false.

  !> Use a factor in the source to account for parallel diffusion
  !> See V R Soloviev and V M Krivtsov 2009 J. Phys. D: Appl. Phys. 42 125208
  integer, parameter :: source_centered = 0
  integer, parameter :: source_flux = 1
  integer, parameter :: source_drift_flux = 2
  integer, parameter :: source_factor = 3

  !> How to compute source terms in the fluid model
  integer :: fluid_source_method = source_flux
  !> The source term can increase compared to the normal case
  logical :: fluid_source_increase = .true.

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
    character(len=40)          :: integrator, source_method
    logical                    :: found, outside_range
    real(dp), allocatable      :: init_charge(:)

    input_file          = "input/n2_transport_data_siglo.txt"
    table_size          = 1000
    integrator          = "rk2"

    call CFG_add_get(cfg, "fluid%input_file", input_file, &
         "Input file with cross sections")
    call CFG_add_get(cfg, "fluid%table_size", table_size, &
         "Size of lookup table for transport coefficients")
    call CFG_add_get(cfg, "fluid%table_max_efield", fluid_max_field, &
         "Maximum electric field in transport coefficient table")
    call CFG_add_get(cfg, "fluid%diffusion_emax", fluid_diffusion_emax, &
         "Disable diffusion parallel to fields above this threshold (V/m)")
    call CFG_add_get(cfg, "fluid%small_density", fluid_small_density, &
         "Smallest allowed density in the fluid model (1/m3)")
    call CFG_add_get(cfg, "fluid%integrator", integrator, &
         "Time integrator (euler, trapezoidal, rk2, rk4)")

    ! CFG settings added for the LEA 
    call CFG_add_get(cfg, "fluid%use_LEA", fluid_use_LEA, &
       "Whether to use the LEA or LFA(default)")
    call CFG_add_get(cfg, "fluid%table_max_e_energy", fluid_max_en, &
         "Maximum electron energy in transport coefficient table")

    if (fluid_use_LEA) then
      call CFG_add(cfg, "fluid%en_mob", "energy[eV]_vs_mu[m2/Vs]", &
            "The name of the mobility coefficient")
      call CFG_add(cfg, "fluid%en_dif", "energy[eV]_vs_dif[m2/s]", &
            "The name of the diffusion coefficient")
      call CFG_add(cfg, "fluid%en_alpha", "energy[eV]_vs_alpha[1/m]", &
            "The name of the eff. ionization coeff.")
      call CFG_add(cfg, "fluid%en_eta", "energy[eV]_vs_eta[1/m]", &
            "The name of the eff. attachment coeff.")
      call CFG_add(cfg, "fluid%en_loss", "energy[eV]_vs_loss[eV/s]", &
            "The name of the  energy losses due to collisions")
    else
      call CFG_add(cfg, "fluid%fld_mob", "efield[V/m]_vs_mu[m2/Vs]", &
            "The name of the mobility coefficient")
      call CFG_add(cfg, "fluid%fld_dif", "efield[V/m]_vs_dif[m2/s]", &
            "The name of the diffusion coefficient")
      call CFG_add(cfg, "fluid%fld_alpha", "efield[V/m]_vs_alpha[1/m]", &
            "The name of the eff. ionization coeff.")
      call CFG_add(cfg, "fluid%fld_eta", "efield[V/m]_vs_eta[1/m]", &
            "The name of the eff. attachment coeff.")
    end if

    call CFG_add(cfg, "fluid%fld_energy", "efield[V/m]_vs_energy[eV]", &
         "The name of the electron energy")
    call CFG_add_get(cfg, "fluid%dt_factor", fluid_dt_factor, &
         "Safety factor (< 1) for the time step")
    call CFG_add_get(cfg, "fluid%limit_velocity", fluid_limit_velocity, &
         "If true, keep the velocity constant for E > table_max_efield")

    source_method = "centered"
    call CFG_add_get(cfg, "fluid%source_method", source_method, &
         "How to compute the source term (default, flux, factor)")
    call CFG_add_get(cfg, "fluid%source_increase", fluid_source_increase, &
         "The source term can increase compared to the normal case")

    select case (source_method)
    case ("centered")
       fluid_source_method = source_centered
    case ("flux")
       fluid_source_method = source_flux
    case ("drift_flux")
       fluid_source_method = source_drift_flux
    case ("factor")
       fluid_source_method = source_factor
    case default
       error stop "Unknown fluid%source_method (options: default, flux, factor)"
    end select

    call CFG_add_get(cfg, "fluid%avoid_drt", fluid_avoid_drt, &
         "Limit fluxes to avoid the dielectric relaxation time")

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
    allocate(init_charge(1:domain_nx))

    fluid_state%a(:, :) = 0.0_dp
    fluid_state%s(:)    = 0.0_dp

    ! Create a lookup table for the model coefficients
    fluid_lkp_fld = LT_create(0.0_dp, fluid_max_field, table_size, 0)
    call CFG_get(cfg, "fluid%fld_energy", data_name)
    call TD_get_td_from_file(input_file, &
          trim(data_name), x_data, y_data)
    outside_range = (maxval(x_data) < fluid_max_field)
    call LT_add_col(fluid_lkp_fld, x_data, y_data)
    ix_energy = fluid_lkp_fld%n_cols
    fluid_num_fld_coef = fluid_lkp_fld%n_cols
    ! Below is the conditional for LEA transport data or LFA transport data
    if (fluid_use_LEA) then
      fluid_lkp_en = LT_create(0.0_dp, fluid_max_en, table_size, 0)
      call CFG_get(cfg, "fluid%en_mob", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data)
      outside_range = (maxval(x_data) < fluid_max_en)
      call LT_add_col(fluid_lkp_en, x_data, y_data)
      ix_mob = fluid_lkp_en%n_cols
   
      call CFG_get(cfg, "fluid%en_dif", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data)
      outside_range = outside_range .or. (maxval(x_data) < fluid_max_en)
      call LT_add_col(fluid_lkp_en, x_data, y_data)
      ix_diff = fluid_lkp_en%n_cols
   
      call CFG_get(cfg, "fluid%en_alpha", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data)
      outside_range = outside_range .or. (maxval(x_data) < fluid_max_en)
      call LT_add_col(fluid_lkp_en, x_data, y_data)
      ix_alpha = fluid_lkp_en%n_cols

      call CFG_get(cfg, "fluid%en_loss", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data)
      outside_range = outside_range .or. (maxval(x_data) < fluid_max_en)
      call LT_add_col(fluid_lkp_en, x_data, y_data)
      ix_loss = fluid_lkp_en%n_cols

      call CFG_get(cfg, "fluid%en_eta", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data, found)
      if (found) then
         outside_range = outside_range .or. & 
            (maxval(x_data) < fluid_max_en)
         call LT_add_col(fluid_lkp_fld, x_data, y_data)
         ix_eta = fluid_lkp_en%n_cols
      end if
      fluid_num_en_coef = fluid_lkp_en%n_cols
    else
      call CFG_get(cfg, "fluid%fld_mob", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data)
      outside_range = (maxval(x_data) < fluid_max_field)
      call LT_add_col(fluid_lkp_fld, x_data, y_data)
      ix_mob = fluid_lkp_fld%n_cols
   
      call CFG_get(cfg, "fluid%fld_dif", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data)
      outside_range = outside_range .or. (maxval(x_data) < fluid_max_field)
      call LT_add_col(fluid_lkp_fld, x_data, y_data)
      ix_diff = fluid_lkp_fld%n_cols
   
      call CFG_get(cfg, "fluid%fld_alpha", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data)
      outside_range = outside_range .or. (maxval(x_data) < fluid_max_field)
      call LT_add_col(fluid_lkp_fld, x_data, y_data)
      ix_alpha = fluid_lkp_fld%n_cols


      call CFG_get(cfg, "fluid%fld_eta", data_name)
      call TD_get_td_from_file(input_file, &
            trim(data_name), x_data, y_data, found)
      if (found) then
         outside_range = outside_range .or. (maxval(x_data) < fluid_max_field)
         call LT_add_col(fluid_lkp_fld, x_data, y_data)
         ix_eta = fluid_lkp_fld%n_cols
      end if
      fluid_num_fld_coef = fluid_lkp_fld%n_cols
    end if


    if (outside_range) then
       print *, "At least some of the input data was not specified"
       print *, "up to fluid%table_max_efield", fluid_max_field
       print *, "OR"
       print *, "up to fluid%table_max_e_energy", fluid_max_en
       error stop
    end if


    if (fluid_use_LEA) then
      ! Determine extrapolation coefficients for transport data
      x = [fluid_max_en, 0.8_dp * fluid_max_en]
   
      ! Linear extrapolation for ionization coefficient
      y = LT_get_col(fluid_lkp_en, ix_alpha, x)
   
      extrap_src_dydx = (y(2) - y(1)) / (x(2) - x(1))
      extrap_src_y0   = y(2)
   
      ! Exponential decay for mobility: mu = c0 * exp(-c1 * E)
      y = LT_get_col(fluid_lkp_en, ix_mob, x)
    else
      ! Determine extrapolation coefficients for transport data
      x = [fluid_max_field, 0.8_dp * fluid_max_field]
   
      ! Linear extrapolation for ionization coefficient
      y = LT_get_col(fluid_lkp_fld, ix_alpha, x)
   
      extrap_src_dydx = (y(2) - y(1)) / (x(2) - x(1))
      extrap_src_y0   = y(2)
   
      ! Exponential decay for mobility: mu = c0 * exp(-c1 * E)
      y = LT_get_col(fluid_lkp_fld, ix_mob, x)
    end if

    extrap_mob_c1 = -log(y(2)/y(1)) / (x(2) - x(1))
    extrap_mob_c0 = y(2) / exp(-extrap_mob_c1 * x(2))


    ! Initialization of electron and ion density
    do n = 1, domain_nx
       xx = (n-0.5_dp) * domain_dx
       fluid_state%a(n, iv_elec) = init_elec_dens(xx)
       fluid_state%a(n, iv_pion) = init_pos_ion_dens(xx)
       fluid_state%a(n, iv_nion) = 0.0_dp
       fluid_state%a(n, iv_en) = 0.0_dp
    end do
    ! Computing the initial field with 0 surface charge to initialize
    ! initial electron energy when using LEA
    init_charge = UC_elem_charge * (fluid_state%a(1:domain_nx, iv_pion) - & 
                  fluid_state%a(1:domain_nx, iv_elec)- &
                  fluid_state%a(1:domain_nx, iv_nion))
    call compute_field(init_charge, (/0.0_dp, 0.0_dp/), 0.0_dp)
    if (fluid_use_LEA) then
      do n = 1, domain_nx
            xx = (n-0.5_dp) * domain_dx
            fluid_state%a(n, iv_en) = fluid_state%a(n, iv_elec)* &
               LT_get_col(fluid_lkp_fld, ix_energy, get_field_at(xx))
      end do
   end if

  end subroutine fluid_initialize

  !> Fill the ghost cells of the fluid variables
  !> TODO: implement this more flexibly
  subroutine set_boundary_conditions(state)
    type(state_t), intent(inout) :: state
    integer                      :: idx, iter(2)

    ! iv_en will just be a 0 array in LFA, so this need not be changed
    iter = (/iv_elec, iv_en/)
    do idx = 1,2
      if (dielectric_present(1)) then
         ! Inside dielectric, density is zero
         state%a(0, iter(idx)) = 0.0_dp
         state%a(-1, iter(idx)) = 0.0_dp
      else
         ! Neumann boundary condition on the left
         state%a(0, iter(idx)) = state%a(1, iter(idx))
         state%a(-1, iter(idx)) = state%a(2, iter(idx))
      end if
   
      if (dielectric_present(2)) then
         state%a(domain_nx+1, iter(idx)) = 0.0_dp
         state%a(domain_nx+2, iter(idx)) = 0.0_dp
      else
         ! Neumann boundary condition on the right
         state%a(domain_nx+1, iter(idx)) = state%a(domain_nx, iter(idx))
         state%a(domain_nx+2, iter(idx)) = state%a(domain_nx-1, iter(idx))
      end if
   end do
   ! Ion motion boundary conditions- only works for positive ions

   if (ion_transport) then
      if (dielectric_present(1)) then
         ! Inside dielectric, density is zero
         state%a(0, iv_pion) = 0.0_dp
         state%a(-1, iv_pion) = 0.0_dp
      else
         ! Neumann boundary condition on the left
         state%a(0, iv_pion) = state%a(1, iv_pion)
         state%a(-1, iv_pion) = state%a(2, iv_pion)
      end if
   
      if (dielectric_present(2)) then
         state%a(domain_nx+1, iv_pion) = 0.0_dp
         state%a(domain_nx+2, iv_pion) = 0.0_dp
      else
         ! Neumann boundary condition on the right
         state%a(domain_nx+1, iv_pion) = state%a(domain_nx, iv_pion)
         state%a(domain_nx+2, iv_pion) = state%a(domain_nx-1, iv_pion)
      end if
   end if
      


  end subroutine set_boundary_conditions

  !> Compute the time derivatives of the fluid variables
  subroutine fluid_derivs(state, time, derivs, dt, dt_max)
    use m_units_constants
    type(state_t), intent(inout)    :: state  !< Current state
    real(dp), intent(in)            :: time   !< Current time
    type(state_t), intent(inout)    :: derivs !< Derivatives
    real(dp), intent(in)            :: dt     !< Current time step
    !> If present, output the maximal time step
    real(dp), intent(out), optional :: dt_max

    type(LT_loc_t), allocatable :: fld_locs(:), fld_locs_cc(:)
    type(LT_loc_t), allocatable :: en_locs(:), en_locs_cc(:)
    integer                     :: n, nx
    real(dp)                    :: surface_charge(2), se
    real(dp), allocatable       :: mob_e(:), diff_e(:), src_e(:), loss_e(:), loss_e_cc(:)
    real(dp), allocatable       :: mean_en(:), elec_fc(:), fld_en(:)
    real(dp), allocatable       :: fld_en_cc(:), mean_en_cc(:)
    real(dp), allocatable       :: mob_i(:), diff_i(:), mob_cc(:)
    real(dp), allocatable       :: flux(:), flux_d(:), max_flux(:)
    real(dp), allocatable       :: flux_ene(:), flux_ene_d(:)
    real(dp), allocatable       :: source(:), tmp_vec(:), src_fac(:)
    real(dp), allocatable       :: flux_ad_cc(:), flux_dif_cc(:), source_fc(:)
    real(dp), allocatable       :: source_ene(:) 
    real(dp), allocatable       :: sigma(:)
    real(dp)                    :: dt_cfl, dt_dif, dt_drt
    real(dp)                    :: drt_fac, tmp
    real(dp), parameter         :: eps = 1e-100_dp, five_third = 5/3.0_dp
    real(dp), parameter         :: four_third = 4/3.0_dp

    nx = domain_nx

    derivs%a(:, :) = 0.0_dp
    derivs%s(:)    = 0.0_dp

    call set_boundary_conditions(state)

    allocate(source(1:nx))
    allocate(flux_ad_cc(1:nx))
    allocate(flux_dif_cc(1:nx))
    allocate(source_fc(1:nx+1))
    allocate(flux(1:nx+1))
    allocate(flux_d(1:nx+1))
    allocate(mob_i(nx+1))
    allocate(diff_i(nx+1))

    if (fluid_use_LEA) then
      allocate(flux_ene(1:nx+1))
      allocate(flux_ene_d(1:nx+1))
      allocate(source_ene(1:nx))
      allocate(loss_e_cc(1:nx))
      allocate(mean_en(1:nx+1))
      allocate(elec_fc(1:nx+1))
      allocate(fld_en(1:nx+1))
    end if

    ! Get electric field
    source = UC_elem_charge * (state%a(1:nx, iv_pion) - state%a(1:nx, iv_elec) &
         - state%a(1:nx, iv_nion))
    surface_charge = UC_elem_charge * (&
         state%s([i_lbound_pion, i_rbound_pion]) - &
         state%s([i_lbound_elec, i_rbound_elec]))
    call compute_field(source, surface_charge, time)

    ! Get locations in the lookup table
    fld_locs = LT_get_loc(fluid_lkp_fld, abs(field_fc))
    if (fluid_use_LEA) then
      ! Obtain energy values corresponding to the field values as fc
      fld_en = LT_get_col_at_loc(fluid_lkp_fld, ix_energy, fld_locs)
      call cell_center_to_face(state%a(0:nx+1, iv_en), mean_en)
      call cell_center_to_face(state%a(0:nx+1, iv_elec), elec_fc)
      mean_en = (mean_en + fluid_small_density*fld_en) / &
                (elec_fc + fluid_small_density)
      ! Below is an alternative approach to the above method
      ! without using cc-fc conversion
      !mean_en = (state%a(0:nx, iv_en) + state%a(1:nx+1, iv_en) + &
      !          2*fluid_small_density *fld_en) / &
      !          (state%a(0:nx, iv_elec) + state%a(1:nx+1, iv_elec) + &
      !          2 * fluid_small_density)
      en_locs = LT_get_loc(fluid_lkp_en, mean_en)
      ! Get LEA transport coefficients
      mob_e  = LT_get_col_at_loc(fluid_lkp_en, ix_mob, en_locs)
      diff_e = LT_get_col_at_loc(fluid_lkp_en, ix_diff, en_locs)
      ! Extrapolate coefficients
      where (abs(mean_en) > fluid_max_en)
         mob_e = extrap_mob(abs(mean_en))
      end where
    else
      ! Get LFA transport coefficients
      mob_e  = LT_get_col_at_loc(fluid_lkp_fld, ix_mob, fld_locs)
      diff_e = LT_get_col_at_loc(fluid_lkp_fld, ix_diff, fld_locs)

      ! Extrapolate coefficients
      where (abs(field_fc) > fluid_max_field)
         mob_e = extrap_mob(abs(field_fc))
      end where
    end if



    ! Compute diffusive electron flux
    flux_d = 0
    call add_diff_flux_1d(nx, n_ghost_cells, state%a(:, iv_elec), &
         diff_e, domain_dx, flux_d)


    if (fluid_diffusion_emax < huge_value) then
       ! If the diffusive flux is parallel to the field, and the field above a
       ! threshold, set the flux to zero
       where (flux_d * field_fc > 0 .and. abs(field_fc) > fluid_diffusion_emax)
          flux_d = 0
       end where
    end if

    ! Compute advective electron flux
    flux = 0
    call add_drift_flux_1d(nx, n_ghost_cells, state%a(:, iv_elec), &
         -mob_e * field_fc, flux)
    

    ! Computing advective and diffusive energy fluxes
    if (fluid_use_LEA) then
      flux_ene_d = 0
      call add_diff_flux_1d(nx, n_ghost_cells, state%a(:, iv_en), &
            diff_e, domain_dx, flux_ene_d)
   
      flux_ene = 0
      call add_drift_flux_1d(nx, n_ghost_cells, state%a(:, iv_en), &
            -mob_e * field_fc, flux_ene)
      flux_ene = flux_ene + flux_ene_d
    end if

   ! Computing the source factor correction terms
    if (fluid_source_method == source_factor) then
       allocate(src_fac(nx))
       src_fac = 1 - (flux_d(1:nx) * sign(1.0_dp, field_fc(1:nx)) + &
            flux_d(2:nx+1) * sign(1.0_dp, field_fc(2:nx+1))) / &
            max(abs(flux(1:nx) + flux(2:nx+1)), 1e-10_dp)
       src_fac = max(0.0_dp, src_fac)
       if (.not. fluid_source_increase) src_fac = min(1.0_dp, src_fac)
    end if

    ! Add electron fluxes

    !call cell_face_to_center(fluid_state%a(1:nx, iv_flux_ad), abs(flux))
    fluid_state%a(1:nx, iv_flux_ad) = flux(1:nx)
    call cell_face_to_center(fluid_state%a(1:nx, iv_flux_dif), flux_d)
    flux = flux + flux_d
    source_fc = mob_e*abs(field_fc)

    if (fluid_avoid_drt) then
       drt_fac = UC_eps0 / max(1e-100_dp, UC_elem_charge * dt)
       allocate(max_flux(nx+1))
       ! max_flux = drt_fac * max(abs(field_fc), 1.0e6_dp)
       do n = 1, nx+1
          tmp = abs(state%a(n, iv_elec) - state%a(n-1, iv_elec)) / &
               max(state%a(n, iv_elec), state%a(n-1, iv_elec), 1e-10_dp)
          tmp = tmp * diff_e(n) * domain_inv_dx / mob_e(n)
          max_flux(n) = drt_fac * max(abs(field_fc(n)), tmp)
       end do

       where (abs(flux) > max_flux)
          flux = sign(max_flux, flux)
       end where
    end if

    select case (fluid_source_method)
    case (source_flux, source_drift_flux)
       if (fluid_use_LEA) then
         ! Compute source terms using the fluxes at cell faces
         src_e = LT_get_col_at_loc(fluid_lkp_en, ix_alpha, en_locs)
         loss_e = LT_get_col_at_loc(fluid_lkp_en, ix_loss, en_locs)
   
         ! Remove attachment coefficient from source
         if (ix_eta > 0) then
            src_e = src_e - LT_get_col_at_loc(fluid_lkp_en, ix_eta, en_locs)
         end if
       else
         ! Compute source terms using the fluxes at cell faces
         src_e = LT_get_col_at_loc(fluid_lkp_fld, ix_alpha, fld_locs)
   
         ! Remove attachment coefficient from source
         if (ix_eta > 0) then
            src_e = src_e - LT_get_col_at_loc(fluid_lkp_fld, ix_eta, fld_locs)
         end if
       end if
       !call cell_face_to_center(fluid_state%a(1:nx, iv_source), src_e*source_fc)
       if (fluid_source_method == source_drift_flux) then
          call cell_face_to_center(source, src_e * abs(flux-flux_d))
       else
          call cell_face_to_center(source, src_e * abs(flux))

          ! Ensure the source term is not bigger than before
          if (.not. fluid_source_increase) then
             allocate(tmp_vec(nx))
             call cell_face_to_center(tmp_vec, src_e * abs(flux-flux_d))
             where (abs(tmp_vec) < abs(source))
                source = tmp_vec
             end where
          end if
       end if
    case (source_centered, source_factor)
       ! Compute source term at cell centers
       fld_locs_cc = LT_get_loc(fluid_lkp_fld, abs(field_cc(1:nx)))
       if (fluid_use_LEA) then
         error stop "Evaluating source at cell centers for LEA! This doesnt work, evaluate at face centers(fluid%source_term= flux)"

         !fld_en_cc = LT_get_col_at_loc(fluid_lkp_fld, ix_energy, fld_locs_cc)
         !mean_en_cc = (state%a(1:nx, iv_en) + &
         !         fluid_small_density *fld_en_cc) / &
         !         (state%a(1:nx, iv_elec) + &
         !         fluid_small_density)
         !en_locs_cc = LT_get_loc(fluid_lkp_fld, mean_en_cc)
         !src_e       = LT_get_col_at_loc(fluid_lkp_en, ix_alpha, en_locs_cc)
         !mob_cc      = LT_get_col_at_loc(fluid_lkp_en, ix_mob, en_locs_cc)
         !loss_e      = LT_get_col_at_loc(fluid_lkp_en, ix_loss, en_locs_cc)
         !source_ene = -(field_cc(1:nx)*mob_cc*field_cc(1:nx)) - loss_e*state%a(1:nx, iv_elec)
   
         ! Remove attachment coefficient from source
         !if (ix_eta > 0) then
         !   src_e = src_e - & 
         !      LT_get_col_at_loc(fluid_lkp_en, ix_eta, en_locs_cc)
         !end if
       else
         src_e       = LT_get_col_at_loc(fluid_lkp_fld, ix_alpha, fld_locs_cc)
         mob_cc      = LT_get_col_at_loc(fluid_lkp_fld, ix_mob, fld_locs_cc)
   
         ! Remove attachment coefficient from source
         if (ix_eta > 0) then
            src_e = src_e - LT_get_col_at_loc(fluid_lkp_fld, ix_eta, fld_locs_cc)
         end if

       end if
   

       source = src_e * abs(mob_cc * field_cc(1:nx) * state%a(1:nx, iv_elec))

       if (fluid_source_method == source_factor) then
          source = source * src_fac
       end if
    end select

    
    ! Computing the energy source term at cell centers
    if (fluid_use_LEA) then
      call cell_face_to_center(source_ene, -1*(flux*field_fc))
      call cell_face_to_center(loss_e_cc, loss_e)
      source_ene = source_ene - state%a(1:nx, iv_elec)*loss_e_cc
      derivs%a(1:nx, iv_en) = source_ene
    else
      derivs%a(1:nx, iv_en) = 0.0_dp
    end if
    derivs%a(1:nx, iv_elec) = source
    derivs%a(1:nx, iv_pion) = source
    fluid_state%a(1:nx, iv_source) = source

    do n = 1, nx
       derivs%a(n, iv_elec) = derivs%a(n, iv_elec) + &
            domain_inv_dx * (flux(n) - flux(n+1))
      if (fluid_use_LEA) then
         derivs%a(n, iv_en) = derivs%a(n, iv_en) + &
               (5.0_dp/3.0_dp)*domain_inv_dx * (flux_ene(n) - flux_ene(n+1))
      end if
    end do

    ! Take into account boundary fluxes
    derivs%s(i_lbound_elec) = -flux(1)
    derivs%s(i_rbound_elec) = flux(nx+1)

    ! Photon fluxes on left and right sides. Assume photons are not absorbed and
    ! have a 50% chance of going left/right.
    ! HTODO: Account for electron energy change due to photoemission
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
       call get_flux_1d(nx, n_ghost_cells, state%a(:, iv_pion), &
            mob_i * field_fc, diff_i, domain_dx, flux)
       do n = 1, nx
          derivs%a(n, iv_pion) = derivs%a(n, iv_pion) + &
               domain_inv_dx * (flux(n) - flux(n+1))
       end do

       ! Take into account boundary fluxes
       derivs%s(i_lbound_pion) = -flux(1)
       derivs%s(i_rbound_pion) = flux(nx+1)
       !if (flux(1) < 0.0_dp) then
       !print *, "Left ion flux is neg"
       !end if

       ! HTODO: Check the correctness of the electron energy BCs
       ! Secondary emission of electrons
       se = max(0.0_dp, -ion_secondary_emission_yield * flux(1))
       derivs%a(1, iv_elec) = derivs%a(1, iv_elec) + se * domain_inv_dx
       derivs%s(i_lbound_elec) = derivs%s(i_lbound_elec) - se
       ! HTODO:Adding the energy to the system due to the se electrons
       if (fluid_use_LEA) then
         derivs%a(1, iv_en) = derivs%a(1, iv_en) + &
               se * domain_inv_dx * &
               five_third*UC_boltzmann_const*(derivs%a(1, iv_en)/ &
               (fluid_small_density + derivs%a(1, iv_elec)))
       end if

       se = max(0.0_dp, ion_secondary_emission_yield * flux(nx+1))
       !if (se > 0.0_dp) then
       !print *, "Adding ion SE right"
       !end if
       derivs%a(nx, iv_elec) = derivs%a(nx, iv_elec) + se * domain_inv_dx
       derivs%s(i_rbound_elec) = derivs%s(i_rbound_elec) - se
       ! HTODO:Adding the energy to the system due to the se electrons
       if (fluid_use_LEA) then
         derivs%a(nx, iv_en) = derivs%a(nx, iv_en) + &
               se * domain_inv_dx * &
               five_third*UC_boltzmann_const*(derivs%a(nx, iv_en)/ &
               (fluid_small_density + derivs%a(nx, iv_elec)))
       end if
      !call cell_face_to_center(flux_cc, flux_d)
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
       if (fluid_avoid_drt) then
          dt_drt = 1e100_dp
       else
          dt_drt = UC_eps0 / (UC_elem_charge * max(eps, maxval(sigma)))
       end if

       ! Take the minimum of the CFL condition with Courant number 0.5 and
       ! the combined CFL-diffusion condition with Courant number 1.0. The
       ! 0.5 is emperical, to have good accuracy (and TVD/positivity) in
       ! combination with the explicit trapezoidal rule
       dt_max = min(0.5_dp * dt_cfl, 1/(1/dt_cfl + 1/dt_dif))

       ! Use a 'safety' factor
       dt_max = fluid_dt_factor * min(dt_max, dt_drt)
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

  !>  Terms defined at cell centers with appropriately filled 
  !> ghost cells will be averaged to cell face values
  subroutine cell_center_to_face(cc_vals, fc_vals)
    real(dp), intent(in) :: cc_vals(:)
    real(dp), intent(inout)    :: fc_vals(:)
    integer                 :: n

    do n = 1, size(cc_vals) - 1
      fc_vals(n) = 0.5_dp*(cc_vals(n) + cc_vals(n+1))
    end do
  end subroutine cell_center_to_face

  !> Write output files
  subroutine fluid_write_output(base_fname, time, dt, ix)
    use m_units_constants
    character(len=*), intent(in) :: base_fname
    real(dp), intent(in)         :: time, dt
    integer, intent(in)          :: ix
    integer                      :: n, nx
    character(len=200)           :: fname
    type(LT_loc_t), allocatable  :: fld_locs_cc(:)
    real(dp), allocatable        :: energy_e(:)
    real(dp)                     :: total_charge
    real(dp)                     :: max_field, deriv_max_field
    real(dp)                     :: velocity, max_ne_loc
    real(dp), save               :: prev_max_field, prev_time, prev_ne_max_loc
    integer                      :: my_unit

    nx = domain_nx

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_fluid_", ix, ".txt"
    open(newunit=my_unit, file=trim(fname))

    ! Get locations in the lookup table
    fld_locs_cc = LT_get_loc(fluid_lkp_fld, abs(field_cc))

    energy_e = LT_get_col_at_loc(fluid_lkp_fld, ix_energy, fld_locs_cc)
    write(my_unit, "(A)") "x field electron pos_ion potential electron_source energy_density mean_energy flux_ad flux_diff"
    do n = 1, domain_nx
       write(my_unit, *) domain_dx * (n-0.5_dp), &
            field_cc(n), &
            fluid_state%a(n, iv_elec), &
            fluid_state%a(n, iv_pion), &
            potential(n), &
            fluid_state%a(n, iv_source), &
            fluid_state%a(n, iv_en), &
            (fluid_state%a(n, iv_en) + fluid_small_density*energy_e(n))/ &
            (fluid_state%a(n, iv_elec) + fluid_small_density), &
            fluid_state%a(n, iv_flux_ad), &
            fluid_state%a(n, iv_flux_dif)
    end do
    close(my_unit)

    write(*, "(A,E9.2,A,A)") " t = ", time, " wrote ", trim(fname)

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_fluid_scalars.txt"
    if (ix == 1) then
       open(newunit=my_unit, file=trim(fname))
       write(my_unit, *) "time dt E_max total_charge sum_elec ", &
            "sum_pion sigmapos_l sigmaneg_l sigmapos_r sigmaneg_r max_elec max_pion deriv_Emax velocity"
    else
       open(newunit=my_unit, file=trim(fname), access='append')
    end if

    total_charge = domain_dx * sum(fluid_state%a(1:nx, iv_pion) &
         - fluid_state%a(1:nx, iv_elec)) &
         - fluid_state%s(i_lbound_elec) - fluid_state%s(i_rbound_elec) &
         + fluid_state%s(i_lbound_pion) + fluid_state%s(i_rbound_pion)

    max_field = maxval(abs(field_fc))


    max_ne_loc = domain_dx*(findloc(fluid_state%a(1:nx, iv_elec), & 
                 maxval(fluid_state%a(1:nx, iv_elec)), dim=1)-0.5_dp)


    if (ix == 1) then
       deriv_max_field = 0.0_dp
       velocity = 0.0_dp
    else
       deriv_max_field = (max_field - prev_max_field) / (time - prev_time)
       velocity = (max_ne_loc - prev_ne_max_loc) / (time - prev_time)
    end if

    write(my_unit, *) time, dt, max_field, total_charge, &
         domain_dx * sum(fluid_state%a(1:nx, iv_elec)), &
         domain_dx * sum(fluid_state%a(1:nx, iv_pion)), &
         fluid_state%s(i_lbound_pion) , fluid_state%s(i_lbound_elec), &
         fluid_state%s(i_rbound_pion) , fluid_state%s(i_rbound_elec), &
         maxval(fluid_state%a(1:nx, iv_elec)), &
         maxval(fluid_state%a(1:nx, iv_pion)), deriv_max_field, &
         velocity
    close(my_unit)

    prev_max_field = max_field
    prev_ne_max_loc = max_ne_loc
    prev_time      = time
  end subroutine fluid_write_output

  !> Advance the fluid state over dt
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
       call fluid_derivs(fluid_state, time, derivs, dt, max_dt)
       fluid_state%a = fluid_state%a + dt * derivs%a
       fluid_state%s = fluid_state%s + dt * derivs%s
    case (trapezoidal)
       substep = fluid_state

       call fluid_derivs(substep, time, derivs, dt)
       substep%a = substep%a + dt * derivs%a
       substep%s = substep%s + dt * derivs%s

       call fluid_derivs(substep, time+dt, derivs, dt, max_dt)
       substep%a = substep%a + dt * derivs%a
       substep%s = substep%s + dt * derivs%s

       fluid_state%a = 0.5_dp * (fluid_state%a + substep%a)
       fluid_state%s = 0.5_dp * (fluid_state%s + substep%s)
    case (rk2)
       substep = fluid_state

       ! Step 1 (at initial time)
       call fluid_derivs(substep, time, derivs, dt)
       substep%a = substep%a + 0.5_dp * dt * derivs%a
       substep%s = substep%s + 0.5_dp * dt * derivs%s

       ! Step 2 (at initial time + dt/2)
       call fluid_derivs(substep, time + 0.5_dp * dt, derivs, dt, max_dt)

       ! Result (at initial time + dt)
       fluid_state%a = fluid_state%a + dt * derivs%a
       fluid_state%s = fluid_state%s + dt * derivs%s
    case (rk4)
       substep = fluid_state

       ! Step 1 (at initial time)
       call fluid_derivs(substep, time, derivs, dt)
       sum_derivs = derivs

       ! Step 2 (at initial time + dt/2)
       substep%a = substep%a + 0.5_dp * dt * derivs%a
       substep%s = substep%s + 0.5_dp * dt * derivs%s
       call fluid_derivs(substep, time + 0.5_dp * dt, derivs, dt)

       sum_derivs%a = sum_derivs%a + 2 * derivs%a
       sum_derivs%s = sum_derivs%s + 2 * derivs%s

       ! Step 3 (at initial time + dt/2)
       substep%a = fluid_state%a + 0.5_dp * dt * derivs%a
       substep%s = fluid_state%s + 0.5_dp * dt * derivs%s
       call fluid_derivs(substep, time + 0.5_dp * dt, derivs, dt)

       sum_derivs%a = sum_derivs%a + 2 * derivs%a
       sum_derivs%s = sum_derivs%s + 2 * derivs%s

       ! Step 4 (at initial time + dt)
       substep%a = fluid_state%a + dt * derivs%a
       substep%s = fluid_state%s + dt * derivs%s
       call fluid_derivs(substep, time + dt, derivs, dt, max_dt)

       sum_derivs%a = sum_derivs%a + derivs%a
       sum_derivs%s = sum_derivs%s + derivs%s

       ! Combine time derivatives at steps
       fluid_state%a = fluid_state%a + dt * one_sixth * sum_derivs%a
       fluid_state%s = fluid_state%s + dt * one_sixth * sum_derivs%s
    case default
       error stop "Unknown time stepping scheme"
    end select

    !where (fluid_state%a < fluid_small_density)
    !   fluid_state%a = 0.0_dp
    !end where
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
