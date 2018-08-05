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
  use m_config
  use m_generic
  use m_fluid_1d
  use m_particle_1d

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  character(len=200) :: output_name

  integer :: it, end_iteration
  integer :: output_ix
  integer :: info_cntr

  real(dp) :: time, end_time
  real(dp) :: dt, dt_next, dt_max, dt_min, dt_limit
  real(dp) :: dt_output
  real(dp) :: time_elapsed
  integer  :: time_now, time_start, count_rate

  logical     :: write_output
  type(CFG_t) :: cfg

  call CFG_update_from_arguments(cfg)

  end_time = 3.0e-9_dp
  call CFG_add_get(cfg, "end%time", end_time, "End time (s)")

  end_iteration = huge(1)
  call CFG_add_get(cfg, "end%iteration", end_iteration, "End iteration")

  dt_next = 1e-12_dp
  call CFG_add_get(cfg, "dt%initial", dt_next, "Initial time step (s)")

  dt_max = 1e-11_dp
  call CFG_add_get(cfg, "dt%max", dt_max, "Maximal time step (s)")

  dt_min = 1e-13_dp
  call CFG_add_get(cfg, "dt%min", dt_max, "Minimal time step (s)")

  output_name = "output/my_sim"
  call CFG_add_get(cfg, "output%filename", output_name, &
       "Base file name for output")

  dt_output = 1.0e-10_dp
  call CFG_add_get(cfg, "output%dt", dt_output, &
         "The time step for writing output")

  print *, "Initialize generic"
  call generic_initialize(cfg)
  print *, "Initialize fluid"
  call fluid_initialize(cfg)
  print *, "Initialize particle"
  call particle_initialize(cfg)
  print *, "done"

  call check_output_folder(trim(output_name))
  call CFG_write(cfg, trim(output_name) // ".cfg", .true.)

  ! Initialize variables
  time      = 0.0_dp
  output_ix = 0
  info_cntr = 1
  it        = 0

  call system_clock(time_start, count_rate)

  ! Here the simulation starts
  do while (time < end_time .and. it < end_iteration)
     call system_clock(time_now)
     time_elapsed = (time_now - time_start) / real(count_rate, dp)

     if (time_elapsed > 10 * info_cntr) then
        info_cntr = info_cntr + 1
        write(*, "(F6.2,A)") 100.0_dp * time/end_time, "% complete"
     end if

     dt = dt_next
     write_output = (time + dt >= output_ix * dt_output)

     if (write_output) then
        dt          = output_ix * dt_output - time
        output_ix = output_ix + 1
     end if

     if (model_type == model_fluid) then
        call fluid_advance(dt, time, dt_limit)

        dt_next = get_new_dt(dt_next, dt_limit)

        if (write_output) then
           call fluid_write_output(trim(output_name), time, output_ix)
        end if
     else
        call particle_advance(dt, time)

        if (write_output) then
           call particle_write_output(trim(output_name), time, output_ix)
        end if
     end if

     it   = it + 1
     time = time + dt
  end do

  write(*, "(A,E10.4,A)") "Simulation ended after ", 1.0d9 * time, " ns"

contains

  !> Initializes everything needed for the simulation
  ! subroutine initialize_all(cfg)
  !   use m_gas
  !   use m_transport_data
  !   use m_cross_sec

  !   type(CFG_t), intent(inout)     :: cfg
  !   integer                        :: nn, n_gas_comp
  !   integer                        :: n_cs_files
  !   character(len=40)              :: model_type_name
  !   character(len=40), allocatable :: comp_names(:)
  !   ! character(len=s_len), allocatable :: cs_files(:)
  !   real(dp), allocatable          :: comp_fracs(:)
  !   ! type(CS_t), allocatable           :: cross_secs(:)
  !   real(dp)                       :: pressure, temperature
  !   ! real(dp) :: max_ev
  !   integer                        :: num_part, max_num_part

  !   call create_sim_config(cfg)      ! Create default parameters for the simulation
  !   call CFG_sort(cfg)

  !   ! call CFG_write(cfg, "output/" // trim(sim_name) // "_config.txt")

  !   allocate(comp_names(n_gas_comp))
  !   allocate(comp_fracs(n_gas_comp))

  !   call CFG_get(cfg, "gas_comp_names", comp_names)
  !   call CFG_get(cfg, "gas_comp_fracs", comp_fracs)
  !   call CFG_get(cfg, "gas_pressure", pressure)
  !   call CFG_get(cfg, "gas_temperature", temperature)

  !   ! Initialize gas and electric field module
  !   call GAS_initialize(comp_names, comp_fracs, pressure, temperature)

  !   select case (MODEL_type)
  !   case (model_fluid)
  !      call fluid_init_cfg(cfg)
  !   end select

  ! end subroutine initialize_all

  !> Create the parameters and their default values for the simulation
  ! subroutine create_sim_config(cfg)
  !   type(CFG_t), intent(inout) :: cfg

  !   ! General simulation parameters

  !   ! Grid parameters
  !   call CFG_add(cfg, "grid_size", 1000, &
  !        "The number of grid cells")
  !   call CFG_add(cfg, "grid_dx", 4.0d-6, &
  !        "The length of a grid cell")

  !   ! Gas parameters
  !   call CFG_add(cfg, "gas_pressure", 1.0D0, &
  !        "The gas pressure (bar)")
  !   call CFG_add(cfg, "gas_temperature", 293.0D0, &
  !        "The gas temperature (Kelvin)")
  !   call CFG_add(cfg, "gas_mixture_name", "N2", &
  !        "The name of the gas mixture used")
  !   call CFG_add(cfg, "gas_comp_names", (/"N2"/), &
  !        "The names of the gases used in the simulation", .true.)

  !   call CFG_add(cfg, "gas_comp_fracs", (/1.0_dp /), &
  !        "The partial pressure of the gases", .true.)

  !   ! Electric field parameters
  !   call CFG_add(cfg, "sim_applied_efield", 1.0D7, &
  !        "The initial electric field")
  !   call CFG_add(cfg, "sim_constant_efield", .false., &
  !        "Whether the electric field is kept constant")

  !   ! Initial conditions
  !   call CFG_add(cfg, "init_cond_name", "gaussian", &
  !        "The type of initial condition")
  !   call CFG_add(cfg, "init_use_neg_ion", .false., &
  !        "Whether to use neg. ions initially")
  !   call CFG_add(cfg, "init_dens", 1.0d15 , &
  !        "The number of initial ion pairs")
  !   call CFG_add(cfg, "init_num_part", 1000 , &
  !        "The number of initial simulation particles")
  !   call CFG_add(cfg, "init_elec_energy", 1.0D0 , &
  !        "The initial energy of the electrons in eV")
  !   call CFG_add(cfg, "init_rel_pos", 0.5D0, &
  !        "The relative position of the initial seed")
  !   call CFG_add(cfg, "init_width", 25.0d-6, &
  !        "The standard deviation used for Gaussian initial profiles")
  !   call CFG_add(cfg, "init_background_density", 0.0D0, &
  !        "The background ion and electron density in 1/m^3")

  !   ! Output parameters
  !   call CFG_add(cfg, "output_num_bins", 200, &
  !        "The number of bins to use for histograms")
  !   call CFG_add(cfg, "output_eedf_min_max_fields", &
  !        (/0.0d0, 1.0d10, 0.0d0, 1.0d6/), &
  !        "The field range(s) over which to get the EEDF")
  !   call CFG_add(cfg, "output_eedf_eV_range", (/0.0d0, 5.0d1/), &
  !        "The energy range over which to get the EEDF")

  !   call CFG_add(cfg, "apm_vel_rel_weight", 1.0d-12, &
  !        "Relative weight of vel. in the k-d tree compared to position")
  !   call CFG_add(cfg, "apm_steps_between", 100, &
  !        "Adapt weight every apm_steps_between steps")
  !   call CFG_add(cfg, "apm_increase_factor", 1.2_dp, &
  !        "Adapt weight if #part increases by this factor")

  !   ! General fluid model parameters
  !   call CFG_add(cfg, "fluid_use_en_mob", .false., &
  !        "Whether to use energy dependent mobility")
  !   call CFG_add(cfg, "fluid_use_en_dif", .false., &
  !        "Whether to use energy dependent diffusion coefficient")
  !   call CFG_add(cfg, "fluid_use_en_src", .false., &
  !        "Whether to use energy dependent source term")
  !   call CFG_add(cfg, "fluid_use_detach", .false., &
  !        "Whether to use electron detachment")
  !   call CFG_add(cfg, "fluid_input_file", "transport_data_nitrogen.txt" , &
  !        "The input file for the fluid models")
  !   call CFG_add(cfg, "fluid_lkptbl_size", 1000, &
  !        "The transport data table size in the fluid model")
  !   call CFG_add(cfg, "fluid_lkptbl_max_efield", 3.0d7, &
  !        "The maximum electric field in the fluid model coefficients")
  !   call CFG_add(cfg, "fluid_lkptbl_max_energy", 1.0d2, &
  !        "The maximum mean energy in eV in the fluid model coefficients")
  !   call CFG_add(cfg, "fluid_small_density", 1.0d0, &
  !        "Regularization density to compute a mean energy")

  !   call CFG_add(cfg, "fluid_en_fld", "energy[eV]_vs_efield[V/m]", &
  !        "The name of the energy vs efield list")
  !   call CFG_add(cfg, "fluid_en_mob", "energy[eV]_vs_mu[m2/Vs]", &
  !        "The name of the mobility coefficient")
  !   call CFG_add(cfg, "fluid_en_dif", "energy[eV]_vs_dif[m2/s]", &
  !        "The name of the diffusion coefficient")
  !   call CFG_add(cfg, "fluid_en_alpha", "energy[eV]_vs_alpha[1/m]", &
  !        "The name of the eff. ionization coeff.")
  !   call CFG_add(cfg, "fluid_en_eta", "energy[eV]_vs_eta[1/m]", &
  !        "The name of the eff. attachment coeff.")
  !   call CFG_add(cfg, "fluid_en_loss", "energy[eV]_vs_loss[eV/s]", &
  !        "The name of the energy loss coeff.")

  !   call CFG_add(cfg, "fluid_fld_mob", "efield[V/m]_vs_mu[m2/Vs]", &
  !        "The name of the mobility coefficient")
  !   call CFG_add(cfg, "fluid_fld_en", "efield[V/m]_vs_energy[eV]", &
  !        "The name of the energy(fld) coefficient")
  !   call CFG_add(cfg, "fluid_fld_dif", "efield[V/m]_vs_dif[m2/s]", &
  !        "The name of the diffusion coefficient")
  !   call CFG_add(cfg, "fluid_fld_alpha", "efield[V/m]_vs_alpha[1/m]", &
  !        "The name of the eff. ionization coeff.")
  !   call CFG_add(cfg, "fluid_fld_eta", "efield[V/m]_vs_eta[1/m]", &
  !        "The name of the eff. attachment coeff.")
  !   call CFG_add(cfg, "fluid_fld_loss", "efield[V/m]_vs_loss[eV/s]", &
  !        "The name of the energy loss coeff.")
  !   call CFG_add(cfg, "fluid_fld_det", "efield[V/m]_vs_det[1/s]", &
  !        "The name of the detachment rate coeff.")
  ! end subroutine create_sim_config

  subroutine check_output_folder(filename)
    character(len=*), intent(in) :: filename
    integer                      :: my_unit, iostate

    open(newunit=my_unit, file=trim(filename)//"DUMMY", iostat=iostate)

    if (iostate /= 0) then
       print *, "Output base file name: ", trim(filename)
       error stop "Cannot write to output folder"
    else
       close(my_unit, status='delete')
    end if
  end subroutine check_output_folder

  real(dp) function get_new_dt(dt, dt_limit)
    real(dp), intent(in) :: dt, dt_limit

    if (dt > dt_limit) then
       get_new_dt = dt_limit
    else
       ! Increase time step at most by 10%
       get_new_dt = min(1.1_dp * dt, dt_limit)
    end if
  end function get_new_dt

end program
