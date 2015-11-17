!> \mainpage
!! This webpage contains documentation for a 1d streamer simulation simulation code.
!! It can simulate these 'ionization waves' or 1d streamer using either a particle or various fluid models.
!!
!! \section secLabel Usage overview
!! Usage with the default configuration is easy: just type "make" to build the program and then execute
!! "streamer_1d" with no arguments. If you want to execute the program with a configuration file,
!! then use "streamer_1d [file]".
!! Input files are located in the "input" directory. Output files are located in the "output" directory.
!! The input file for the particle model contains cross sections.
!! The input file for the fluid models contains transport data (mobility, diffusion, ionization rate etc.).
!! The output files have the name and type of the simulation prepended (e.g., simName_type_000001.txt for the first one).
!! They contain the following columns: position, electric field, electron density, ion density, mean energy and possibly more.


program streamer_1d
  use m_model_choice
  use m_config
  use m_efield_1d
  use m_output_1d
  use m_fluid_dd_1d
  use m_particle_1d
  use m_phys_domain

  implicit none

  integer, parameter   :: dp    = kind(0.0d0)
  integer, parameter   :: s_len = 100
  character(len=s_len) :: sim_name

  integer              :: steps, n_steps_apm
  integer              :: output_cntr
  integer              :: info_cntr
  integer              :: prev_apm_part

  real(dp)             :: small_dens
  real(dp)             :: sim_time, sim_end_time
  real(dp)             :: dt, max_dt
  real(dp)             :: output_dt, min_field
  real(dp)             :: time_now, time_start
  real(dp)             :: apm_increase

  logical              :: do_apm
  logical              :: stop_dens, stop_fld, stop_time, stop_sim
  type(CFG_t)          :: cfg

  call initialize_all(cfg) ! Initialize all necessary modules

  ! Initialize variables
  prev_apm_part = PM_get_num_sim_part()
  sim_time      = 0.0_dp
  output_cntr   = 0
  info_cntr     = 0
  steps         = 0
  stop_sim      = .false.
  stop_dens     = .false.
  min_field     = CFG_fget_real(cfg, "sim_min_field")
  dt            = CFG_fget_real(cfg, "sim_initial_dt")
  max_dt        = CFG_fget_real(cfg, "sim_max_dt")
  output_dt     = CFG_fget_real(cfg, "output_interval")
  sim_end_time  = CFG_fget_real(cfg, "sim_end_time")
  n_steps_apm   = CFG_fget_int(cfg, "apm_steps_between")
  small_dens    = CFG_fget_real(cfg, "fluid_small_density")
  apm_increase  = CFG_fget_real(cfg, "apm_increase_factor")

  call OUT_write_coeffs(sim_name)
  call cpu_time(time_start)

  ! Here the simulation starts
  do
     call cpu_time(time_now)
     if (time_now - time_start > 10 * info_cntr .or. stop_sim) then
        if (MODEL_type == MODEL_part) then
           print *, "Number of physical particles", PM_get_num_real_part()
           write(*, "(I8,A,E8.2,A,E8.2,A,E8.2,A,E8.2)") steps, &
                " -- t = ", sim_time, ", dt = ", dt, &
                ", wct = ", time_now - time_start, &
                ", n_part = ", real(PM_get_num_sim_part(), dp)
        else
           write(*, "(I8,A,E8.2,A,E8.2,A,E8.2)") steps, &
                " -- t = ", sim_time, ", dt = ", dt, &
                ", wct = ", time_now - time_start
        end if
        info_cntr = info_cntr + 1
     end if

     if (sim_time >= output_cntr * output_dt .or. stop_sim) then
        output_cntr = output_cntr + 1
        call OUT_write_vars(sim_name, output_cntr, sim_time)
     end if

     if (stop_sim) exit ! Exit after output

     select case (MODEL_type)
     case (MODEL_part)
        do_apm = (n_steps_apm > 0 .and. mod(steps, n_steps_apm) == 0) &
             .or. PM_get_num_sim_part() > prev_apm_part * apm_increase
        call PM_advance(dt, do_apm)
        if (do_apm) prev_apm_part = PM_get_num_sim_part()

        sim_time = sim_time + dt
        stop_dens = (PM_max_edens_at_boundary() > small_dens)
     case (MODEL_fluid_lfa, MODEL_fluid_ee)
        call FL_advance(sim_time, dt)
        stop_dens = (FL_max_edens_at_boundary() > small_dens)
     end select

     stop_fld  = EF_get_min_field() < min_field
     stop_time = sim_time > sim_end_time

     if (stop_dens) print *, "Stopping because discharge reached boundary"
     if (stop_fld)  print *, "Stopping because field is getting too low"
     if (stop_time) print *, "Stopping because end time was reached"
     stop_sim = stop_dens .or. stop_fld .or. stop_time

     steps = steps + 1
  end do

  write(*, "(A,E10.4,A)") "Simulation ended after ", 1.0d9 * sim_time, " ns"

contains

  !> Initializes everything needed for the simulation
  subroutine initialize_all(cfg)
    use m_gas
    use m_init_cond_1d
    use m_transport_data
    use m_cross_sec

    type(CFG_t), intent(inout)        :: cfg
    integer                           :: ix
    integer                           :: nn, n_gas_comp
    integer                           :: n_cs_files
    character(len=s_len)              :: MODEL_type_name
    character(len=s_len)              :: cfg_name, tmp_name, prev_name
    character(len=s_len), allocatable :: comp_names(:)
    character(len=s_len), allocatable :: cs_files(:)
    real(dp), allocatable             :: comp_fracs(:)
    type(CS_t), allocatable           :: cross_secs(:)

    call create_sim_config(cfg)      ! Create default parameters for the simulation
    call CFG_sort(cfg)

    sim_name = ""
    prev_name = ""
    do ix = 1, command_argument_count()
       call get_command_argument(ix, cfg_name)
       call CFG_read_file(cfg, trim(cfg_name))

       call CFG_get(cfg, "sim_name", tmp_name)
       if (sim_name == "") then
          sim_name = tmp_name
       else if (tmp_name /= "" .and. tmp_name /= prev_name) then
          sim_name = trim(sim_name) // "_" // trim(tmp_name)
       end if
       prev_name = tmp_name
    end do

    call CFG_get(cfg, "sim_type", MODEL_type_name)
    call CFG_write(cfg, "output/" // trim(sim_name) // "_config.txt")

    call MODEL_initialize(MODEL_type_name)
    call PD_set(CFG_fget_real(cfg, "grid_dx"), CFG_fget_int(cfg, "grid_size"))

    n_gas_comp = CFG_fget_size(cfg, "gas_comp_names")
    if (n_gas_comp /= CFG_fget_size(cfg, "gas_comp_fracs")) &
         stop "streamer_1d: gas_comp_names/fracs have unequal size"

    allocate(comp_names(n_gas_comp))
    allocate(comp_fracs(n_gas_comp))

    call CFG_get(cfg, "gas_comp_names", comp_names)
    call CFG_get(cfg, "gas_comp_fracs", comp_fracs)

    ! Initialize gas and electric field module
    call GAS_initialize(comp_names, comp_fracs, &
         CFG_fget_real(cfg, "gas_pressure"), &
         CFG_fget_real(cfg, "gas_temperature"))
    call EF_initialize(cfg)
    call INIT_init(cfg)
    call OUT_init(cfg)

    select case (MODEL_type)
    case (MODEL_part)
       call CFG_get_size(cfg, "gas_cs_files", n_cs_files)
       allocate(cs_files(n_cs_files))
       call CFG_get(cfg, "gas_cs_files", cs_files)
       if (n_cs_files /= n_gas_comp) then
          print *, "streamer_1d: variables gas_cs_files and", &
               " gas_component_fracs have unequal size"
          stop
       end if

       do nn = 1, n_gas_comp
          call CS_add_from_file("input/" // trim(cs_files(nn)), &
               trim(GAS_comp_names(nn)), GAS_comp_fracs(nn) * GAS_number_dens, &
               CFG_fget_real(cfg, "part_max_ev"), cross_secs)
       end do

       call CS_write_summary(cross_secs, &
            "output/" // trim(sim_name) // "_cs_summary.txt")

       call PM_initialize(cfg, cross_secs, CFG_fget_int(cfg, "init_num_part"), &
            CFG_fget_int(cfg, "part_max_num"))

    case (MODEL_fluid_lfa, MODEL_fluid_ee)
       call FL_init_cfg(cfg)
    end select

  end subroutine initialize_all

  !> Create the parameters and their default values for the simulation
  subroutine create_sim_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    ! General simulation parameters
    call CFG_add(cfg, "sim_type", "fluid_lfa", &
         "The type of simulation to run: part, fluid_lfa, fluid_ee")
    call CFG_add(cfg, "sim_end_time", 3.0D-9, &
         "The desired endtime in seconds of the simulation")
    call CFG_add(cfg, "sim_initial_dt", 3.0D-13, &
         "The initial/fixed timestep in seconds")
    call CFG_add(cfg, "sim_max_dt", 3.0D-13, &
         "The maximal timestep in seconds")
    call CFG_add(cfg, "sim_min_field", -1.0d99, &
         "The minimum required electric field")
    call CFG_add(cfg, "sim_name", "sim", &
         "The name of the simulation")

    ! Grid parameters
    call CFG_add(cfg, "grid_size", 1000, &
         "The number of grid cells")
    call CFG_add(cfg, "grid_dx", 4.0d-6, &
         "The length of a grid cell")

    ! Gas parameters
    call CFG_add(cfg, "gas_pressure", 1.0D0, &
         "The gas pressure (bar)")
    call CFG_add(cfg, "gas_temperature", 293.0D0, &
         "The gas temperature (Kelvin)")
    call CFG_add(cfg, "gas_mixture_name", "N2", &
         "The name of the gas mixture used")
    call CFG_add(cfg, "gas_comp_names", (/"N2"/), &
         "The names of the gases used in the simulation", .true.)
    call CFG_add(cfg, "gas_cs_files", (/"cross_sections_nitrogen.txt"/), &
         "The files in which to find cross section data for each gas", .true.)
    call CFG_add(cfg, "gas_comp_fracs", (/1.0_dp /), &
         "The partial pressure of the gases", .true.)

    ! Electric field parameters
    call CFG_add(cfg, "sim_applied_efield", 1.0D7, &
         "The initial electric field")
    call CFG_add(cfg, "sim_constant_efield", .false., &
         "Whether the electric field is kept constant")

    ! Initial conditions
    call CFG_add(cfg, "init_cond_name", "gaussian", &
         "The type of initial condition")
    call CFG_add(cfg, "init_use_neg_ion", .false., &
         "Whether to use neg. ions initially")
    call CFG_add(cfg, "init_dens", 1.0d15 , &
         "The number of initial ion pairs")
    call CFG_add(cfg, "init_num_part", 1000 , &
         "The number of initial simulation particles")
    call CFG_add(cfg, "init_elec_energy", 1.0D0 , &
         "The initial energy of the electrons in eV")
    call CFG_add(cfg, "init_rel_pos", 0.5D0, &
         "The relative position of the initial seed")
    call CFG_add(cfg, "init_width", 25.0d-6, &
         "The standard deviation used for Gaussian initial profiles")
    call CFG_add(cfg, "init_background_density", 0.0D0, &
         "The background ion and electron density in 1/m^3")

    ! Output parameters
    call CFG_add(cfg, "output_interval", 1.0D-10, &
         "The timestep for writing output")
    call CFG_add(cfg, "output_head_density", 2.0D18, &
         "The density level tracked as streamer head")
    call CFG_add(cfg, "output_num_bins", 200, &
         "The number of bins to use for histograms")
    call CFG_add(cfg, "output_eedf_min_max_fields", &
         (/0.0d0, 1.0d10, 0.0d0, 1.0d6/), &
         "The field range(s) over which to get the EEDF")
    call CFG_add(cfg, "output_eedf_eV_range", (/0.0d0, 5.0d1/), &
         "The energy range over which to get the EEDF")

    call CFG_add(cfg, "apm_part_per_cell", 2.0d2, &
         "The desired number of particles per cell")
    call CFG_add(cfg, "apm_vel_rel_weight", 1.0d-12, &
         "Relative weight of vel. in the k-d tree compared to position")
    call CFG_add(cfg, "apm_steps_between", 100, &
         "Adapt weight every apm_steps_between steps")
    call CFG_add(cfg, "apm_increase_factor", 1.2_dp, &
         "Adapt weight if #part increases by this factor")

    ! Particle model related parameters
    call CFG_add(cfg, "part_lkptbl_size", 20*1000, &
         "The size of the lookup table for the collision rates")
    call CFG_add(cfg, "part_max_num", 25*1000*1000, &
         "The maximum number of particles allowed per task")
    call CFG_add(cfg, "part_max_ev", 1000.0D0, &
         "The maximum energy in eV for particles in the simulation")

    ! General fluid model parameters
    call CFG_add(cfg, "fluid_use_en_mob", .false., &
         "Whether to use energy dependent mobility")
    call CFG_add(cfg, "fluid_use_en_dif", .false., &
         "Whether to use energy dependent diffusion coefficient")
    call CFG_add(cfg, "fluid_use_en_src", .false., &
         "Whether to use energy dependent source term")
    call CFG_add(cfg, "fluid_use_detach", .false., &
         "Whether to use electron detachment")
    call CFG_add(cfg, "fluid_input_file", "transport_data_nitrogen.txt" , &
         "The input file for the fluid models")
    call CFG_add(cfg, "fluid_lkptbl_size", 1000, &
         "The transport data table size in the fluid model")
    call CFG_add(cfg, "fluid_lkptbl_max_efield", 3.0d7, &
         "The maximum electric field in the fluid model coefficients")
    call CFG_add(cfg, "fluid_lkptbl_max_energy", 1.0d2, &
         "The maximum mean energy in eV in the fluid model coefficients")
    call CFG_add(cfg, "fluid_small_density", 1.0d0, &
         "Regularization density to compute a mean energy")

    call CFG_add(cfg, "fluid_en_fld", "energy[eV]_vs_efield[V/m]", &
         "The name of the energy vs efield list")
    call CFG_add(cfg, "fluid_en_mob", "energy[eV]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    call CFG_add(cfg, "fluid_en_dif", "energy[eV]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(cfg, "fluid_en_alpha", "energy[eV]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(cfg, "fluid_en_eta", "energy[eV]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")
    call CFG_add(cfg, "fluid_en_loss", "energy[eV]_vs_loss[eV/s]", &
         "The name of the energy loss coeff.")

    call CFG_add(cfg, "fluid_fld_mob", "efield[V/m]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    call CFG_add(cfg, "fluid_fld_en", "efield[V/m]_vs_energy[eV]", &
         "The name of the energy(fld) coefficient")
    call CFG_add(cfg, "fluid_fld_dif", "efield[V/m]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(cfg, "fluid_fld_alpha", "efield[V/m]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(cfg, "fluid_fld_eta", "efield[V/m]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")
    call CFG_add(cfg, "fluid_fld_loss", "efield[V/m]_vs_loss[eV/s]", &
         "The name of the energy loss coeff.")
    call CFG_add(cfg, "fluid_fld_det", "efield[V/m]_vs_det[1/s]", &
         "The name of the detachment rate coeff.")
  end subroutine create_sim_config

end program streamer_1d
