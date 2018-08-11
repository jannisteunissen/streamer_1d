!> Module for 1D particle simulations, making use of the particle_core routines
module m_particle_1d
  use m_generic
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp          = kind(0.0d0)

  integer, parameter :: PM_num_vars = 4
  integer, parameter :: iv_elec     = 1
  integer, parameter :: iv_pion     = 2
  integer, parameter :: iv_en       = 3
  integer, parameter :: iv_ppc      = 4

  integer, parameter :: num_scalars = 4
  integer, parameter :: i_lbound_elec = 1
  integer, parameter :: i_rbound_elec = 2
  integer, parameter :: i_lbound_pion = 3
  integer, parameter :: i_rbound_pion = 4

  integer, parameter :: n_ghost_cells = 2

  real(dp), allocatable :: PM_vars(:, :)
  real(dp), allocatable :: PM_scalars(:)

  !> Transverse area of the 1D domain (automatically determined)
  real(dp) :: PM_transverse_area
  !> Volume of a grid cell
  real(dp) :: PM_cell_volume
  !> Inverse volume of a grid cell
  real(dp) :: PM_inv_cell_volume

  !> Goal number of particles per cell
  real(dp) :: PM_part_per_cell = 1.0e2_dp

  !> Weight factor for velocity when merging particles
  real(dp) :: PM_v_merge_weight = 1.0e-12_dp

  !> Adapt weights when particle number changes by this factor
  real(dp) :: PM_merge_factor = 1.1_dp

  !> Maximal mobility to estimate dielectric relaxation time
  real(dp) :: PM_max_mobility_drt = 0.2_dp

  real(dp)                       :: merge_prev_npart
  character(len=10), allocatable :: gas_components(:)
  real(dp), allocatable          :: gas_fractions(:)
  type(PC_t)                     :: pc
  type(PC_events_t)              :: events

  ! Public routines
  public :: particle_initialize
  public :: particle_advance
  public :: particle_write_output

contains

  !> Initialize particle module
  subroutine particle_initialize(cfg)
    use m_config
    use m_random
    use m_units_constants
    use m_cross_sec

    type(CFG_t), intent(inout) :: cfg

    character(len=200)      :: cs_file
    type(CS_t), allocatable :: cross_secs(:)
    integer                 :: n_gas_comp
    integer                 :: n_part_init, n_part_max
    integer                 :: n, i, num_part, tbl_size
    real(dp)                :: x, sum_elec_dens, max_elec_dens
    real(dp)                :: elec_dens, ion_dens
    real(dp)                :: pos(3), vel(3), accel(3), max_eV

    cs_file             = "input/cross_sections_siglo.txt"
    max_eV              = 1e3_dp
    n_part_init         = 10*1000
    n_part_max          = 10*1000*1000
    tbl_size            = 1000

    call CFG_add(cfg, "gas%components", ["N2"], &
         "Gas components", .true.)
    call CFG_add(cfg, "gas%fractions", [1.0_dp], &
         "Gas components fractions", .true.)
    call CFG_get_size(cfg, "gas%components", n_gas_comp)
    allocate(gas_components(n_gas_comp))
    allocate(gas_fractions(n_gas_comp))

    call CFG_get(cfg, "gas%components", gas_components)
    call CFG_get(cfg, "gas%fractions", gas_fractions)

    call CFG_add_get(cfg, "gas%cross_sections", cs_file, &
         "File with cross section data")
    call CFG_add_get(cfg, "particle%max_eV", max_eV, &
         "Maximum energy (eV) for particles")
    call CFG_add_get(cfg, "particle%initial_number", n_part_init, &
         "Min. number of initial simulation particles")
    call CFG_add_get(cfg, "particle%max_number", n_part_max, &
         "The maximum number of particles")
    call CFG_add_get(cfg, "particle%part_per_cell", PM_part_per_cell, &
         "The desired number of particles per cell")
    call CFG_add_get(cfg, "particle%v_merge_weight", PM_v_merge_weight, &
         "Velocity weight factor (s) for particle merging")
    call CFG_add_get(cfg, "particle%merge_increase_factor", &
         PM_merge_factor, &
         "Adapt weights if particle count increases by this factor")
    call CFG_add_get(cfg, "particle%table_size", tbl_size, &
         "Size of the lookup table for collision rates")
    call CFG_add_get(cfg, "particle%max_mobility_drt", PM_max_mobility_drt, &
         "Maximum mobility to estimate dielectric relaxation time")

    ! Exit here if we don't use the particle model
    if (model_type /= model_particle) return

    do n = 1, n_gas_comp
       call CS_add_from_file(trim(cs_file), trim(gas_components(n)), &
            gas_fractions(n) * gas_number_dens, max_eV, cross_secs)
    end do

    if (size(cross_secs) < 1) &
         error stop "No cross sections given to particle module!"

    sum_elec_dens = 0.0_dp
    max_elec_dens = 0.0_dp
    do n = 1, domain_nx
       x             = (n - 0.5_dp) * domain_dx
       elec_dens     = init_elec_dens(x)
       sum_elec_dens = sum_elec_dens + elec_dens
       max_elec_dens = max(max_elec_dens, elec_dens)
    end do

    if (sum_elec_dens < epsilon(1.0_dp)) &
         error stop "Initial electron density seems to be zero!"

    PM_cell_volume     = max(n_part_init / sum_elec_dens, &
         PM_part_per_cell / max_elec_dens)
    PM_transverse_area = PM_cell_volume / domain_dx
    PM_inv_cell_volume = 1 / PM_cell_volume

    allocate(PM_vars(-1:domain_nx+2, PM_num_vars))
    PM_vars(:, :) = 0.0_dp
    allocate(PM_scalars(num_scalars))
    PM_scalars(:) = 0.0_dp

    pc%accel_function => accel_func
    call pc%initialize(UC_elec_mass, cross_secs, tbl_size, max_eV, n_part_max)

    pc%outside_check => outside_check

    where (pc%colls(:)%type == CS_ionize_t .or. &
         pc%colls(:)%type == CS_attach_t)
       pc%coll_is_event(:) = .true.
    end where

    ! Initialization of electron and ion density
    do n = 1, domain_nx
       x         = (n - 1) * domain_dx
       elec_dens = init_elec_dens(x)
       ion_dens  = init_pos_ion_dens(x)

       num_part = nint(elec_dens * PM_cell_volume)

       ! Uniform position in the cell
       pos(2:3) = [0.0_dp, 0.0_dp]
       vel(:)   = 0
       ! Accel is set after computing electric field
       accel(:) = 0.0_dp

       do i = 1, num_part
          pos(1) = x + (i-0.5_dp) * domain_dx / num_part
          call pc%create_part(pos, vel, accel, 1.0_dp, 0.0_dp)
       end do

       pm_vars(n, iv_pion) = num_part / PM_cell_volume
    end do

    call PM_update_efield(0.0_dp)
    call pc%set_accel()

    merge_prev_npart = pc%get_num_sim_part()
  end subroutine particle_initialize

  subroutine set_elec_density()
    PM_vars(:, iv_elec) = 0.0_dp
    call pc%loop_ipart(add_elec_to_dens)
  end subroutine set_elec_density

  subroutine set_elec_en_density()
    PM_vars(:, iv_en) = 0.0_dp
    call pc%loop_ipart(add_elec_en_to_dens)
  end subroutine set_elec_en_density

  subroutine set_particle_per_cell()
    PM_vars(:, iv_ppc) = 0.0_dp
    call pc%loop_ipart(add_particle_to_count)
  end subroutine set_particle_per_cell

  subroutine add_to_dens_cic(amount, x, iv)
    real(dp), intent(in) :: amount, x
    integer, intent(in)  :: iv
    integer              :: low_ix
    real(dp)             :: tmp, dens_dif

    ! Do linear interpolation, note that dens(1) is defined at 0.5 * dx
    tmp    = x * domain_inv_dx - 0.5_dp
    low_ix = floor(tmp) + 1
    ! Coefficient for upper point between 0 and 1
    tmp    = low_ix - tmp

    dens_dif = amount * PM_inv_cell_volume

    if (low_ix < 0 .or. low_ix > domain_nx) then
       error stop "a particle is outside the computational domain"
    end if

    if (low_ix == 0) then
       PM_vars(1, iv) = PM_vars(1, iv) + dens_dif
    else if (low_ix == domain_nx) then
       PM_vars(domain_nx, iv) = PM_vars(domain_nx, iv) + dens_dif
    else
       PM_vars(low_ix, iv)   = PM_vars(low_ix, iv) + (1.0_dp - tmp) * dens_dif
       PM_vars(low_ix+1, iv) = PM_vars(low_ix+1, iv) + tmp * dens_dif
    end if
  end subroutine add_to_dens_cic

  subroutine add_elec_to_dens(my_part)
    type(PC_part_t), intent(in) :: my_part
    call add_to_dens_cic(my_part%w, my_part%x(1), iv_elec)
  end subroutine add_elec_to_dens

  subroutine add_particle_to_count(my_part)
    type(PC_part_t), intent(in) :: my_part
    call add_to_dens_cic(PM_cell_volume, my_part%x(1), iv_ppc)
  end subroutine add_particle_to_count

  subroutine add_elec_en_to_dens(my_part)
    use m_units_constants
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: energy
    energy = 0.5_dp * my_part%w * UC_elec_mass * sum(my_part%v**2)
    call add_to_dens_cic(energy, my_part%x(1), iv_en)
  end subroutine add_elec_en_to_dens

  subroutine PM_update_efield(time)
    use m_units_constants
    real(dp), intent(in) :: time
    real(dp) :: source(1:domain_nx)
    real(dp) :: surface_charge(2)

    call set_elec_density()

    source = UC_elem_charge * (PM_vars(1:domain_nx, iv_pion) - &
         PM_vars(1:domain_nx, iv_elec))
    surface_charge = UC_elem_charge * (&
         PM_scalars([i_lbound_pion, i_rbound_pion]) - &
         PM_scalars([i_lbound_elec, i_rbound_elec]))

    call compute_field(source, surface_charge, time)
  end subroutine PM_update_efield

  function accel_func(my_part) result(accel)
    use m_units_constants
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: accel(3)

    accel = get_accel(my_part%x(1))
  end function accel_func

  function get_accel(x) result(accel)
    use m_units_constants
    real(dp), intent(in) :: x
    real(dp)             :: accel(3)
    real(dp), parameter  :: accel_fac = UC_elec_charge / UC_elec_mass

    ! Only acceleration in the z-direction
    accel(1) = accel_fac * get_field_at(x)
    accel(2:3) = 0.0_dp
  end function get_accel

  subroutine particle_advance(dt, time, dt_limit)
    use m_units_constants
    real(dp), intent(in)  :: dt
    real(dp), intent(in)  :: time
    real(dp), intent(out) :: dt_limit
    integer               :: n_in, n_out, ix(1)
    real(dp)              :: field_change, vmax, sigma
    real(dp)              :: dt_cfl, dt_drt
    real(dp), allocatable :: tmp(:)
    real(dp), parameter   :: eps = 1e-100_dp

    n_in = pc%get_num_sim_part()
    call pc%advance_openmp(dt, events)

    if (ion_transport) then
       call update_ion_density(dt)
    end if

    tmp = field_fc
    call PM_update_efield(time)

    ! Update acceleration and correct positions and velocities
    call pc%after_mover(dt)

    ! Handle ionization events etc.
    call handle_events(events)
    n_out = pc%get_num_sim_part()

    ! Determine relative deviation in electric field
    tmp = abs(field_fc - tmp) / max(maxval(abs(field_fc)), eps)
    ix = maxloc(tmp)
    field_change = tmp(ix(1))

    if (n_out > PM_merge_factor * merge_prev_npart) then
       call pc%merge_and_split([.true., .false., .false.], &
         PM_v_merge_weight, .true., get_desired_weight, &
         domain_dx, PC_merge_part_rxv, PC_split_part)
       merge_prev_npart = pc%get_num_sim_part()
       call PM_update_efield(time)
       call pc%set_accel()
    end if

    ! Determine CFL-like time step
    vmax   = maxval(abs(pc%particles(1:pc%n_part)%v(1)))
    dt_cfl = domain_dx / max(vmax, eps)

    ! Dielectric relaxation time estimate
    sigma = PM_max_mobility_drt * UC_elem_charge * &
         max(eps, maxval(PM_vars(1:domain_nx, iv_elec)))
    dt_drt = UC_eps0 / sigma

    dt_limit = 0.9_dp * min(dt_cfl, dt_drt)
    ! if (field_change > 1e-1_dp) then
    !    dt_limit = min(dt_limit, dt * 1e-1_dp / field_change)
    ! end if

    ! print *, ix(1), maxval(abs(field_fc)), field_fc(ix(1))
    ! print *, dt, dt_cfl, dt_drt
  end subroutine particle_advance

  subroutine update_ion_density(dt)
    real(dp), intent(in)  :: dt
    integer               :: n, nx, n_elec
    real(dp), allocatable :: mob_i(:), diff_i(:), flux(:)
    real(dp)              :: x(3), v(3), a(3), w, se

    nx = domain_nx
    allocate(flux(1:nx+1))
    allocate(mob_i(nx+1))
    allocate(diff_i(nx+1))

    diff_i(:) = ion_diffusion
    mob_i(:)  = ion_mobility

    ! Compute ion flux
    call get_flux_1d(PM_vars(:, iv_pion), mob_i * field_fc, &
         diff_i, domain_dx, flux, nx, n_ghost_cells)

    do n = 1, nx
       PM_vars(n, iv_pion) = PM_vars(n, iv_pion) + &
            domain_inv_dx * (flux(n) - flux(n+1)) * dt
    end do

    ! Take into account boundary fluxes
    PM_scalars(i_lbound_pion) = PM_scalars(i_lbound_pion) &
         - flux(1) * dt
    PM_scalars(i_rbound_pion) = PM_scalars(i_rbound_pion) &
         + flux(nx+1) * dt

    ! Secondary emission of electrons on lower boundary
    se = -ion_secondary_emission_yield * flux(1) * dt
    n_elec = random_round(se * PM_transverse_area)

    ! Create at most PM_part_per_cell electrons
    if (n_elec > PM_part_per_cell) then
       w      = n_elec / PM_part_per_cell
       n_elec = nint(PM_part_per_cell)
    else
       w = 1.0_dp
    end if

    x(1)   = 0.1_dp * domain_dx
    x(2:3) = 0.0_dp
    v(:)   = 0.0_dp
    a(:)   = get_accel(x(1))

    ! Create electrons
    do n = 1, n_elec
       call pc%create_part(x, v, a, w, 0.0_dp)
    end do
    PM_scalars(i_lbound_elec) = PM_scalars(i_lbound_elec) &
         - n_elec * w / PM_transverse_area

    ! Secondary emission of electrons on upper boundary
    se = ion_secondary_emission_yield * flux(nx+1) * dt
    n_elec = random_round(se * PM_transverse_area)

    if (n_elec > PM_part_per_cell) then
       w      = n_elec / PM_part_per_cell
       n_elec = nint(PM_part_per_cell)
    else
       w = 1.0_dp
    end if

    x(1) = domain_length - 0.1_dp * domain_dx
    a(:) = get_accel(x(1))

    ! Create electrons
    do n = 1, n_elec
       call pc%create_part(x, v, a, w, 0.0_dp)
    end do

    PM_scalars(i_rbound_elec) = PM_scalars(i_rbound_elec) &
         - n_elec * w / PM_transverse_area
  end subroutine update_ion_density

  integer function random_round(x)
    real(dp), intent(in) :: x

    random_round = floor(x)
    if (pc%rng%unif_01() < x - random_round) then
       random_round = random_round + 1
    end if
  end function random_round

  subroutine handle_events(events)
    use m_cross_sec
    type(PC_events_t), intent(inout) :: events
    integer                          :: n, n_se_photons
    real(dp)                         :: w, sum_ionizations
    real(dp)                         :: x(3), v(3), a(3), tmp

    sum_ionizations = 0.0_dp

    do n = 1, events%n_stored
       select case (events%list(n)%ctype)
       case (CS_ionize_t)
          w = events%list(n)%part%w
          call add_to_dens_cic(w, events%list(n)%part%x(1), iv_pion)
          sum_ionizations = sum_ionizations + w

       case (CS_attach_t)
          call add_to_dens_cic(-events%list(n)%part%w, &
               events%list(n)%part%x(1), iv_pion)

       case (PC_particle_went_out)
          if (events%list(n)%part%x(1) < 0.0_dp) then
             PM_scalars(i_lbound_elec) = PM_scalars(i_lbound_elec) + &
                  events%list(n)%part%w / PM_transverse_area
          else
             PM_scalars(i_rbound_elec) = PM_scalars(i_rbound_elec) + &
                  events%list(n)%part%w / PM_transverse_area
          end if

       case default
          print *, events%list(n)%ctype
          error stop "Unknown event"
       end select
    end do

    events%n_stored = 0

    ! Create electrons at walls due to secondary emission. Create at most
    ! PM_part_per_cell electrons at each boundary cell.
    if (se_photons_per_ionization > 0.0_dp) then
       ! Determine number of photons per side
       tmp = 0.5_dp * sum_ionizations * se_photons_per_ionization
       n_se_photons = random_round(tmp)

       if (n_se_photons > PM_part_per_cell) then
          w            = n_se_photons / PM_part_per_cell
          n_se_photons = nint(PM_part_per_cell)
       else
          w = 1.0_dp
       end if

       v(:)   = 0.0_dp
       x(1)   = 0.1_dp * domain_dx
       x(2:3) = 0.0_dp
       a(:)   = get_accel(x(1))

       ! Electrons will only be released if the field points towards the surface
       if (a(1) > 0.0_dp) then
          do n = 1, n_se_photons
             ! Create on lower boundary
             PM_scalars(i_lbound_elec) = PM_scalars(i_lbound_elec) - &
                  w / PM_transverse_area
             call pc%create_part(x, v, a, w, 0.0_dp)
          end do
       end if

       x(1) = domain_length - 0.1_dp * domain_dx
       a(:)   = get_accel(x(1))

       if (a(1) < 0.0_dp) then
          do n = 1, n_se_photons
             ! Create on upper boundary
             PM_scalars(i_rbound_elec) = PM_scalars(i_rbound_elec) - &
                  w / PM_transverse_area
             call pc%create_part(x, v, a, w, 0.0_dp)
          end do
       end if
    end if

  end subroutine handle_events

  real(dp) function get_desired_weight(my_part)
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: elec_dens
    integer                     :: ix

    ix = floor(my_part%x(1) * domain_inv_dx) + 1
    elec_dens = PM_vars(ix, iv_elec)
    get_desired_weight = elec_dens * PM_cell_volume / PM_part_per_cell
    get_desired_weight = max(1.0_dp, get_desired_weight)
  end function get_desired_weight

  subroutine particle_write_output(base_fname, time, dt, ix)
    use m_units_constants
    character(len=*), intent(in) :: base_fname
    real(dp), intent(in)         :: time, dt
    integer, intent(in)          :: ix
    integer                      :: n, nx, my_unit
    character(len=200)           :: fname
    real(dp)                     :: total_charge
    real(dp), allocatable        :: eedf(:, :)

    nx = domain_nx
    call set_elec_density()
    call set_elec_en_density()
    call set_particle_per_cell()

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_particle_", ix, ".txt"
    open(newunit=my_unit, file=trim(fname))

    write(my_unit, "(A)") "x field electron pos_ion energy potential num_part"
    do n = 1, domain_nx
       write(my_unit, *) domain_dx * (n-0.5_dp), &
            field_cc(n), &
            PM_vars(n, iv_elec), &
            PM_vars(n, iv_pion), &
            PM_vars(n, iv_en), &
            potential(n), &
            PM_vars(n, iv_ppc)
    end do
    close(my_unit)

    write(*, "(A,E9.2,A,A)") " t = ", time, " wrote ", trim(fname)

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_particle_scalars.txt"
    if (ix == 1) then
       open(newunit=my_unit, file=trim(fname))
       write(my_unit, *) "time dt max_field total_charge ", &
            "sum_elec sum_pion sigma_l sigma_r max_elec max_pion"
    else
       open(newunit=my_unit, file=trim(fname), access='append')
    end if

    total_charge = domain_dx * sum(PM_vars(1:nx, iv_pion) &
         - PM_vars(1:nx, iv_elec)) &
         - PM_scalars(i_lbound_elec) - PM_scalars(i_rbound_elec) &
         + PM_scalars(i_lbound_pion) + PM_scalars(i_rbound_pion)

    write(my_unit, *) time, dt, maxval(abs(field_fc)), total_charge, &
         domain_dx * sum(PM_vars(1:nx, iv_elec)), &
         domain_dx * sum(PM_vars(1:nx, iv_pion)), &
         PM_scalars(i_lbound_pion) - PM_scalars(i_lbound_elec), &
         PM_scalars(i_rbound_pion) - PM_scalars(i_rbound_elec), &
         maxval(PM_vars(1:nx, iv_elec)), &
         maxval(PM_vars(1:nx, iv_pion))
    close(my_unit)

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_particle_eedf_", &
         ix, ".txt"
    open(newunit=my_unit, file=trim(fname))
    write(my_unit, *) "energy(eV) count"
    call get_eedf(eedf, [0.0_dp, 1e2_dp], [0.0_dp, 1e8_dp], 100)
    do n = 1, size(eedf, 2)
       write(my_unit, *) eedf(:, n)
    end do
    close(my_unit)

  end subroutine particle_write_output

  subroutine get_eedf(eedf, energy_range, field_range, n_bins)
    use m_units_constants
    real(dp), intent(inout), allocatable :: eedf(:,:)
    real(dp), intent(in)                 :: energy_range(2)
    real(dp), intent(in)                 :: field_range(2)
    integer, intent(in)                  :: n_bins
    integer                              :: ix
    real(dp)                             :: fac, delta_en
    real(dp)                             :: accel_range(2)

    allocate(eedf(2, n_bins))

    ! Set centers of the bins
    do ix = 1, n_bins
       fac = (ix - 0.5_dp) / n_bins
       eedf(1, ix) = energy_range(1) + &
            fac * (energy_range(2) - energy_range(1))
    end do

    ! Convert to range of velocity**2
    eedf(1, :) = eedf(1, :) * UC_elec_volt / (0.5_dp * pc%get_mass())
    accel_range = field_range * UC_elem_charge / UC_elec_mass

    call pc%histogram(get_vel_squared, select_part_by_accel, &
         accel_range, eedf(1,:), eedf(2,:))

    ! Normalize
    delta_en = (energy_range(2) - energy_range(1)) / n_bins
    eedf(2, :) = eedf(2, :) / (sum(eedf(2, :)) * delta_en)

    ! Convert to energy range
    eedf(1, :) = eedf(1, :) * 0.5_dp * pc%get_mass() / UC_elec_volt
  end subroutine get_eedf

  real(dp) function get_vel_squared(my_part)
    type(PC_part_t), intent(in) :: my_part
    get_vel_squared = sum(my_part%v**2)
  end function get_vel_squared

  pure logical function select_part_by_accel(my_part, accel_range)
    type(PC_part_t), intent(in) :: my_part
    real(dp), intent(in)        :: accel_range(:)
    real(dp)                    :: accel
    accel                = abs(my_part%a(1))
    select_part_by_accel = (accel >= accel_range(1) .and. &
         accel <= accel_range(2))
  end function select_part_by_accel

  pure integer function outside_check(my_part)
    type(PC_part_t), intent(in) :: my_part

    if (my_part%x(1) < 0.0_dp .or. my_part%x(1) > domain_length) then
       outside_check = 1
    else
       outside_check = 0
    end if
  end function outside_check

end module m_particle_1d
