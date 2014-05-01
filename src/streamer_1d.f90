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

   implicit none

   integer, parameter :: dp = kind(0.0d0)
   character(len=80) :: sim_name

   integer           :: sim_type
   integer           :: steps, n_steps_apm
   integer           :: output_cntr
   integer           :: info_cntr
   integer :: prev_apm_part

   real(dp)          :: err_goal, small_dens
   real(dp)          :: sim_time, sim_end_time
   real(dp)          :: dt, max_dt, new_dt
   real(dp)          :: output_dt, min_field
   real(dp)          :: time_now, time_start
   real(dp) :: apm_increase

   logical           :: do_apm, stop_sim

   call initialize_all() ! Initialize all necessary modules (routine contained in program)

   ! Initialize variables
   prev_apm_part = PC_get_num_sim_part()
   sim_time      = 0.0_dp
   output_cntr   = 0
   info_cntr     = 0
   steps         = 0
   stop_sim      = .false.
   min_field     = CFG_get_real("sim_min_field")
   dt            = CFG_get_real("sim_initial_dt")
   max_dt        = CFG_get_real("sim_max_dt")
   output_dt     = CFG_get_real("output_interval")
   sim_end_time  = CFG_get_real("sim_end_time")
   err_goal      = CFG_get_real("sim_rel_error_goal")
   n_steps_apm   = CFG_get_int("apm_steps_between")
   small_dens    = CFG_get_real("fluid_small_density")
   apm_increase  = CFG_get_real("apm_increase_factor")

   call OUT_write_coeffs(sim_name, sim_type)
   call cpu_time(time_start)

   ! Here the simulation starts
   do
      call cpu_time(time_now)
      if (time_now - time_start > 10 * info_cntr .or. stop_sim) then
         if (stop_sim) print *, "Stopping simulation"
         if (sim_type == MODEL_part) then
            print *, "Number of physical particles", PC_get_num_real_part()
            write(*, "(I8,A,E8.2,A,E8.2,A,E8.2,A,E8.2)") steps, " -- t = ", sim_time, ", dt = ", dt, &
                 ", wct = ", time_now - time_start, ", n_part = ", real(PC_get_num_sim_part(), dp)
         else
            write(*, "(I8,A,E8.2,A,E8.2,A,E8.2)") steps, " -- t = ", sim_time, ", dt = ", dt, &
                 ", wct = ", time_now - time_start
         end if
         info_cntr = info_cntr + 1
      end if

      if (sim_time >= output_cntr * output_dt .or. stop_sim) then
         output_cntr = output_cntr + 1
         call OUT_write_vars(sim_name, sim_type, output_cntr, sim_time)
      end if

      if (stop_sim) exit ! Exit after output

      select case (sim_type)
      case (MODEL_part)
         do_apm = (n_steps_apm > 0 .and. mod(steps, n_steps_apm) == 0) &
              .or. PC_get_num_sim_part() > prev_apm_part * apm_increase
         call PM_advance(dt, do_apm)
         if (do_apm) prev_apm_part = PC_get_num_sim_part()

         sim_time = sim_time + dt
         stop_sim = (PM_max_edens_at_boundary() > small_dens)
      case (MODEL_fluid_lfa, MODEL_fluid_ee)
         call FL_advance(sim_time, dt, max_dt, new_dt, err_goal)
         dt = new_dt
         stop_sim = (FL_max_edens_at_boundary() > small_dens)
      end select

      stop_sim = stop_sim .or. sim_time > sim_end_time &
           .or. (EF_get_min_field() < min_field)
      steps = steps + 1
   end do

   write(*, "(A,E10.4,A)") "The simulation has ended after ", 1.0d9 * sim_time, " ns"
   print *, CFG_get_real("sim_applied_efield"), sim_time, "LOG"
contains

   !> Initializes everything needed for the simulation
   subroutine initialize_all()
      use m_gas
      use m_init_cond_1d
      use m_transport_data
      use m_cross_sec

      integer :: nn, n_gas_comp
      character(len=100) :: sim_type_name, cs_file, cfg_name
      character(len=20), allocatable :: gas_comp_names(:)
      real(dp), allocatable :: gas_comp_fracs(:)
      type(CS_type), allocatable :: cross_secs(:)

      call create_sim_config()      ! Create default parameters for the simulation
      do nn = 1, command_argument_count()
         call get_command_argument(nn, cfg_name)
         call CFG_read_file(trim(cfg_name))
      end do
      call CFG_sort()

      call CFG_get("sim_name", sim_name)
      call CFG_get("sim_type", sim_type_name)
      call CFG_write("output/" // trim(sim_name) // "_config.txt")

      call MODEL_initialize(sim_type_name)
      sim_type = MODEL_get_type()    ! Get the type of simulation model

      n_gas_comp = CFG_get_size("gas_component_names")
      if (n_gas_comp /= CFG_get_size("gas_component_fractions")) then
         print *, "streamer_1d: variables gas_component_names and gas_component_fracs have unequal size"
         stop
      end if

      allocate(gas_comp_names(n_gas_comp))
      allocate(gas_comp_fracs(n_gas_comp))

      call CFG_get("gas_component_names", gas_comp_names)
      call CFG_get("gas_component_fractions", gas_comp_fracs)

      ! Initialize gas and electric field module
      call GAS_initialize(gas_comp_names, gas_comp_fracs, &
           CFG_get_real("gas_pressure"), CFG_get_real("gas_temperature"))
      call EF_initialize()
      call INIT_init()
      call OUT_init()

      select case (sim_type)
      case (MODEL_part)
         if (CFG_get_size("gas_crosssec_files") /= n_gas_comp) then
            print *, "streamer_1d: variables gas_crosssec_files and gas_component_fracs have unequal size"
            stop
         end if

         do nn = 1, n_gas_comp
            cs_file = CFG_get_string("gas_crosssec_files", nn)
            call CS_read_file("input/" // trim(cs_file), trim(gas_comp_names(nn)), 1.0_dp, &
                 gas_comp_fracs(nn) * GAS_get_number_dens(), CFG_get_real("part_max_energy_eV"))
         end do

         call CS_get_cross_secs(cross_secs)
         call CS_write_summary("output/" // trim(sim_name) // "_cs_summary.txt")
         call CS_write_all("output/" // trim(sim_name) // "_cs.txt")
         call PM_initialize(cross_secs, CFG_get_int("init_num_part"))
      case (MODEL_fluid_lfa, MODEL_fluid_ee)
         call FL_init_cfg(sim_type, "input/" // CFG_get_string("fluid_input_file"), &
              CFG_get_string("gas_mixture_name"))
      end select

   end subroutine initialize_all

   !> Create the parameters and their default values for the simulation
   subroutine create_sim_config()

      ! General simulation parameters
      call CFG_add("sim_type", "fluid", "The type of simulation to run, options are particle, fluid_min, fluid_ee")
      call CFG_add("sim_end_time", 1.0D-9, "The desired endtime in seconds of the simulation")
      call CFG_add("sim_initial_dt", 3.0D-13, "The initial/fixed timestep in seconds")
      call CFG_add("sim_max_dt", 3.0D-13, "The maximal timestep in seconds")
      call CFG_add("sim_min_field", 1.0d99, "The minimum required electric field")
      call CFG_add("sim_rel_error_goal", 1.0D-4, "Desired relative error in the solution between consecutive steps")
      call CFG_add("sim_name", "my_sim", "The name of the simulation")

      ! Grid parameters
      call CFG_add("grid_num_points", 2000, "The number of grid cells")
      call CFG_add("grid_delta_x", 2.0d-6, "The length of a grid cell")

      ! Gas parameters
      call CFG_add("gas_pressure", 1.0D0, "The gas pressure (bar)")
      call CFG_add("gas_temperature", 293.0D0, "The gas temperature (Kelvin)")
      call CFG_add("gas_mixture_name", "N2", "The name of the gas mixture used")
      call CFG_add("gas_component_names", (/"N2"/), "The names of the gases used in the simulation", .true.)
      call CFG_add("gas_crosssec_files", (/"cross_sections_nitrogen.txt"/), &
           & "The files in which to find cross section data for each gas", .true.)
      call CFG_add("gas_component_fractions", (/1.0_dp /), &
           & "The partial pressure of the gases (as if they were ideal gases)", .true.)
      ! call CFG_add("gas_crosssec_scaling", 1.0D0, "Scale factor for the cross sections in the input data", .true.)

      ! Electric field parameters
      call CFG_add("sim_applied_efield", 1.0D7, "The initial electric field")
      call CFG_add("sim_constant_efield", .false., "Whether the electric field is kept constant")

      ! Initial conditions
      call CFG_add("init_cond_name", "gaussian", "The type of initial condition")
      call CFG_add("init_use_neg_ion", .false., "Whether to use neg. ions initially")
      call CFG_add("init_dens", 1.0d17 , "The number of initial ion pairs")
      call CFG_add("init_num_part", 1000 , "The number of initial simulation particles")
      call CFG_add("init_elec_energy", 1.0D0 , "The initial energy of the electrons in eV")
      call CFG_add("init_rel_pos", 0.5D0, "The relative position of the initial seed")
      call CFG_add("init_width", 25.0d-6, "The standard deviation used for Gaussian initial profiles")
      call CFG_add("init_background_density", 0.0D0, "The background ion and electron density in 1/m^3")
      call CFG_add("init_cutoff_density", 1.0D15, "The cutoff density for the initial condition")

      ! Output parameters
      call CFG_add("output_interval", 1.0D-10, "The timestep for writing output")
      call CFG_add("output_head_density", 2.0D18, "The density level tracked as streamer head")
      call CFG_add("output_num_bins", 200, "The number of bins to use for histograms")
      call CFG_add("output_eedf_min_max_fields", (/0.0d0, 1.0d10, 0.0d0, 1.0d6/), "The field range(s) over which to get the EEDF")
      call CFG_add("output_eedf_eV_range", (/0.0d0, 5.0d1/), "The energy range over which to get the EEDF")

      call CFG_add("apm_part_per_cell", 1.0d2, "The desired number of particles per cell")
      call CFG_add("apm_vel_rel_weight", 1.0d-6, "Relative weight of vel. in the k-d tree compared to position")
      call CFG_add("apm_steps_between", 100, "Adapt weight every apm_steps_between steps")
      call CFG_add("apm_increase_factor", 1.2_dp, "Adapt weight if #part increases by this factor")

      ! Particle model related parameters
      call CFG_add("part_lkptbl_size", 20*1000, "The size of the lookup table for the collision rates")
      call CFG_add("part_max_number_of", 25*1000*1000, "The maximum number of particles allowed per task")
      call CFG_add("part_max_energy_eV", 1000.0D0, "The maximum energy in eV for particles in the simulation")

      ! General fluid model parameters
      call CFG_add("fluid_use_en_mob", .false., "Whether to use energy dependent mobility")
      call CFG_add("fluid_use_en_dif", .false., "Whether to use energy dependent diffusion coefficient")
      call CFG_add("fluid_use_en_src", .false., "Whether to use energy dependent source term")
      call CFG_add("fluid_use_detach", .false., "Whether to use electron detachment")
      call CFG_add("fluid_input_file", "transport_data_nitrogen.txt" , "The input file for the fluid models")
      call CFG_add("fluid_lkptbl_size", 1000, "The transport data table size in the fluid model")
      call CFG_add("fluid_lkptbl_max_efield", 3.0d7, "The maximum electric field in the fluid model coefficients")
      call CFG_add("fluid_lkptbl_max_energy", 1.0d2, "The maximum mean energy in eV in the fluid model coefficients")
      call CFG_add("fluid_small_density", 1.0d0, "Regularization density to compute a mean energy")

      call CFG_add("fluid_en_fld", "energy[eV]_vs_efield[V/m]", "The name of the energy vs efield list")
      call CFG_add("fluid_en_mob", "energy[eV]_vs_mu[m2/Vs]", "The name of the mobility coefficient")
      call CFG_add("fluid_en_dif", "energy[eV]_vs_dif[m2/s]", "The name of the diffusion coefficient")
      call CFG_add("fluid_en_alpha", "energy[eV]_vs_alpha[1/m]", "The name of the eff. ionization coeff.")
      call CFG_add("fluid_en_eta", "energy[eV]_vs_eta[1/m]", "The name of the eff. attachment coeff.")
      call CFG_add("fluid_en_loss", "energy[eV]_vs_loss[eV/s]", "The name of the energy loss coeff.")

      call CFG_add("fluid_fld_mob", "efield[V/m]_vs_mu[m2/Vs]", "The name of the mobility coefficient")
      call CFG_add("fluid_fld_en", "efield[V/m]_vs_energy[eV]", "The name of the energy(fld) coefficient")
      call CFG_add("fluid_fld_dif", "efield[V/m]_vs_dif[m2/s]", "The name of the diffusion coefficient")
      call CFG_add("fluid_fld_alpha", "efield[V/m]_vs_alpha[1/m]", "The name of the eff. ionization coeff.")
      call CFG_add("fluid_fld_eta", "efield[V/m]_vs_eta[1/m]", "The name of the eff. attachment coeff.")
      call CFG_add("fluid_fld_loss", "efield[V/m]_vs_loss[eV/s]", "The name of the energy loss coeff.")
      call CFG_add("fluid_fld_det", "efield[V/m]_vs_det[1/s]", "The name of the detachment rate coeff.")
   end subroutine create_sim_config

end program streamer_1d
