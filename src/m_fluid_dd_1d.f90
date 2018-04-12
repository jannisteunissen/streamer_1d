!> A not so accurate implementation of a drift-diffusion-reaction plasma fluid model

module m_fluid_dd_1d
   use m_lookup_table
   use m_transport_schemes
   use m_phys_domain
   use m_photoi_1d
   use m_random

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)
   integer               :: FL_num_vars
   integer               :: FL_iv_elec, FL_iv_ion, FL_iv_en, FL_iv_nion
   integer               :: FL_if_mob, FL_if_dif, FL_if_src, FL_if_en, &
        FL_if_att, FL_if_det, FL_num_fld_coef
   integer               :: FL_ie_mob, FL_ie_dif, FL_ie_src, FL_ie_att, &
        FL_ie_loss, FL_num_en_coef
   logical               :: FL_use_en, FL_use_en_mob, FL_use_en_dif, &
        FL_use_en_src, FL_use_detach

   real(dp), allocatable :: FL_vars(:,:)
   real(dp)              :: FL_surface_charge
   real(dp)              :: FL_small_dens
   real(dp)              :: FL_max_energy
   real(dp)              :: FL_pos_ion_mob, FL_pos_ion_diff
   real(dp)              :: FL_max_E, FL_max_ne, FL_max_imp, FL_E_t_max, FL_ne_t_max, FL_imp_t_max
   real(dp), allocatable :: photo(:)
   real(dp)              :: photo_yield

   type(LT_table_t) :: FL_lkp_fld, FL_lkp_en
   procedure(TS_dd_1d_type), pointer :: FL_transport_scheme
   
   type(rng_t) :: rng

   integer :: rng_int4_seed(4) = [8123, 91234, 12399, 293434]

   public :: FL_init_cfg
   public :: FL_advance
   public :: FL_get_output
   public :: FL_get_coeffs
   public :: FL_max_edens_at_boundary
   public :: FL_surface_charge
   public :: FL_max_E, FL_max_ne, FL_max_imp
   public :: FL_E_t_max, FL_ne_t_max, FL_imp_t_max
   public :: photoi_initialize
   public :: photoi_set_src


contains

   subroutine FL_init_cfg(cfg)
      use m_transport_data
      use m_init_cond_1d
      use m_config
      use m_model_choice
      use m_efield_1d

      type(CFG_t), intent(in) :: cfg
      integer, parameter      :: name_len = 100
      character(len=name_len) :: input_file, gas_name
      integer                 :: n, table_size
      real(dp)                :: max_energy, max_efield, xx
      real(dp), allocatable   :: x_data(:), y_data(:)
      character(len=100)      :: data_name

      input_file = "input/" // CFG_fget_string(cfg, "fluid_input_file")
      gas_name = CFG_fget_string(cfg, "gas_mixture_name")

      FL_transport_scheme => TS_dd_1d_up_fl

      call CFG_get(cfg, "fluid_lkptbl_size", table_size)
      call CFG_get(cfg, "fluid_lkptbl_max_energy", max_energy)
      call CFG_get(cfg, "fluid_lkptbl_max_efield", max_efield)
      call CFG_get(cfg, "fluid_small_density", FL_small_dens)
      call CFG_get(cfg, "pos_ion_mob", FL_pos_ion_mob)
      call CFG_get(cfg, "pos_ion_diff", FL_pos_ion_diff)
      FL_max_energy = max_energy
      

      if (MODEL_type == MODEL_fluid_ee) then
         FL_use_en     = .true.
         call CFG_get(cfg, "fluid_use_en_mob", FL_use_en_mob)
         call CFG_get(cfg, "fluid_use_en_dif", FL_use_en_dif)
         call CFG_get(cfg, "fluid_use_en_src", FL_use_en_src)
         call CFG_get(cfg, "fluid_use_detach", FL_use_detach)
      else
         FL_use_en     = .false.
         FL_use_en_mob = .false.
         FL_use_en_dif = .false.
         FL_use_en_src = .false.
         FL_use_detach = .false.
      end if

      if (allocated(FL_vars)) deallocate(FL_vars)
      FL_num_vars = 3
      FL_iv_elec = 1
      FL_iv_ion = 2
      FL_iv_nion = 3

      if (FL_use_en) then
         FL_num_vars = FL_num_vars + 1
         FL_iv_en = FL_num_vars
      end if  

      allocate(FL_vars(PD_grid_size, FL_num_vars))
      allocate(photo(PD_grid_size))

      ! Create a lookup table for the model coefficients
      FL_lkp_fld = LT_create(0.0_dp, max_efield, table_size, 0)
      if (FL_use_en) then
         FL_lkp_en = LT_create(0.0_dp, max_energy, table_size, 0)
         call CFG_get(cfg, "fluid_en_loss", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_en, x_data, y_data)
         FL_ie_loss = LT_get_num_cols(FL_lkp_en)
      end if

      if (FL_use_en_mob) then
         call CFG_get(cfg, "fluid_en_mob", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_en, x_data, y_data)
         FL_ie_mob = LT_get_num_cols(FL_lkp_en)
      else
         call CFG_get(cfg, "fluid_fld_mob", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_fld, x_data, y_data)
         FL_if_mob = LT_get_num_cols(FL_lkp_fld)
      end if

      if (FL_use_en_dif) then
         call CFG_get(cfg, "fluid_en_dif", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_en, x_data, y_data)
         FL_ie_dif = LT_get_num_cols(FL_lkp_en)
      else
         call CFG_get(cfg, "fluid_fld_dif", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_fld, x_data, y_data)
         FL_if_dif = LT_get_num_cols(FL_lkp_fld)
      end if

      if (FL_use_en_src) then
         call CFG_get(cfg, "fluid_en_alpha", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_en, x_data, y_data)
         FL_ie_src = LT_get_num_cols(FL_lkp_en)
         call CFG_get(cfg, "fluid_en_eta", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_en, x_data, y_data)
         FL_ie_att = LT_get_num_cols(FL_lkp_en)
      else
         call CFG_get(cfg, "fluid_fld_alpha", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_fld, x_data, y_data)
         FL_if_src = LT_get_num_cols(FL_lkp_fld)
         call CFG_get(cfg, "fluid_fld_eta", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_fld, x_data, y_data)
         FL_if_att = LT_get_num_cols(FL_lkp_fld)
      end if

      if (FL_use_detach) then
         call CFG_get(cfg, "fluid_fld_det", data_name)
         call TD_get_td_from_file(input_file, gas_name, &
              trim(data_name), x_data, y_data)
         call LT_add_col(FL_lkp_fld, x_data, y_data)
         FL_if_det= LT_get_num_cols(FL_lkp_fld)
      end if

      call CFG_get(cfg, "fluid_fld_en", data_name)
      call TD_get_td_from_file(input_file, gas_name, &
           trim(data_name), x_data, y_data)
      call LT_add_col(FL_lkp_fld, x_data, y_data)
      FL_if_en = LT_get_num_cols(FL_lkp_fld)

      FL_num_fld_coef = LT_get_num_cols(FL_lkp_fld)
      FL_num_en_coef  = LT_get_num_cols(FL_lkp_en)
      
      ! Initialization of electron and ion density
      FL_surface_charge = 0.0_dp
      do n = 1, PD_grid_size
         xx = (n-1) * PD_dx
         call INIT_get_back(xx, FL_vars(n, FL_iv_elec), FL_vars(n, FL_iv_ion))
         call INIT_get_elec_dens(xx, FL_vars(n, FL_iv_elec))
         call INIT_get_ion_dens(xx, FL_vars(n, FL_iv_ion))
         call INIT_get_nion_dens(xx, FL_vars(n, FL_iv_nion))
         if (FL_use_en) then
            FL_vars(n, FL_iv_en) = FL_vars(n, FL_iv_elec) * &
                 LT_get_col(FL_lkp_fld, FL_if_en, EF_get_at_pos((n-1)*PD_dx))
         end if
      end do
      
   
      FL_max_E = 0.0_dp
      FL_max_ne = 0.0_dp
      FL_max_imp = 0.0_dp
      FL_E_t_max = 0.0_dp
      FL_ne_t_max = 0.0_dp
      FL_imp_t_max = 0.0_dp

   end subroutine FL_init_cfg

   real(dp) function FL_max_edens_at_boundary()
      FL_max_edens_at_boundary = &
           max(FL_vars(1, FL_iv_elec), FL_vars(PD_grid_size, FL_iv_elec))
   end function FL_max_edens_at_boundary

   subroutine FL_time_derivs(vars, time, time_derivs)
      use m_efield_1d
      use m_units_constants
      use m_init_cond_1d
      real(dp), intent(in)        :: vars(:,:), time
      real(dp), intent(out)       :: time_derivs(:,:)

      integer                     :: n_cc, iz_d
      real(dp), parameter         :: five_third = 5 / 3.0_dp
      real(dp)                    :: inv_delta_x
      real(dp), dimension(size(vars,1)-1) :: &
           mob_c, dif_c, src_c, att_c, det_c, flux, fld, fld_en, en_loss, mean_en, ones
      real(dp)                    :: source(size(vars,1))
      type(LT_loc_t) :: fld_locs(size(vars,1)-1), en_locs(size(vars,1)-1)

      n_cc = size(vars, 1)
      iz_d     = int(INIT_DI/PD_dx+1) 
      inv_delta_x = 1.0_dp / PD_dx
      ones = 1.0_dp

      ! Get electric field
      source = (vars(:, FL_iv_ion) - vars(:, FL_iv_elec) - vars(:, FL_iv_nion)) * UC_elem_charge
      source(iz_d+1:PD_grid_size) = 0.0_dp
      call EF_compute_and_get_st(source, fld, FL_surface_charge)

      ! Get locations in the lookup table
      fld_locs = LT_get_loc(FL_lkp_fld, abs(fld))
      if (FL_use_en) then
         ! There is a regularization: at density zero, we use the energy corresponding to the efield
         fld_en = LT_get_col_at_loc(FL_lkp_fld, FL_if_en, fld_locs)
         mean_en = (vars(1:n_cc-1, FL_iv_en) + vars(2:n_cc, FL_iv_en) + 2 * FL_small_dens * fld_en) &
              / (vars(1:n_cc-1, FL_iv_elec) + vars(2:n_cc, FL_iv_elec) + 2 * FL_small_dens)
         en_locs = LT_get_loc(FL_lkp_en, mean_en)
         en_loss = LT_get_col_at_loc(FL_lkp_en, FL_ie_loss, en_locs)
      end if

      if (FL_use_en_mob) then
         mob_c = LT_get_col_at_loc(FL_lkp_en, FL_ie_mob, en_locs)
      else
         mob_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_mob, fld_locs)
      end if

      if (FL_use_en_dif) then
         dif_c = LT_get_col_at_loc(FL_lkp_en, FL_ie_dif, en_locs)
      else
         dif_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_dif, fld_locs)
      end if

      if (FL_use_en_src) then
         src_c = LT_get_col_at_loc(FL_lkp_en, FL_ie_src, en_locs)
         att_c = LT_get_col_at_loc(FL_lkp_en, FL_ie_att, en_locs)
      else
         src_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_src, fld_locs)
         att_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_att, fld_locs)
      end if
      if (FL_use_detach) det_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_det, fld_locs)

      ! ~~~ Electron transport ~~~
      call FL_transport_scheme(vars(:, FL_iv_elec), -mob_c * fld, dif_c * inv_delta_x, flux)
      

      ! ~~~ Ionization source ~~~ TODO: try different formulations
      call set_stagg_source_1d(source, (src_c-att_c) * abs(flux))
      time_derivs(1:iz_d, FL_iv_elec) = source(1:iz_d) + photo(1:iz_d)
      time_derivs(1:iz_d, FL_iv_ion) = source(1:iz_d) + photo(1:iz_d)
      time_derivs(iz_d, FL_iv_elec) = time_derivs(iz_d, FL_iv_elec) + sum(photo(iz_d+1:n_cc)) * photo_yield
      time_derivs(iz_d+1:n_cc, FL_iv_elec) = 0.0_dp
      time_derivs(iz_d+1:n_cc, FL_iv_ion) = 0.0_dp
      time_derivs(:, FL_iv_nion) = 0

      ! ~~~ Attachment source ~~~ TODO: try different formulations
      !call set_stagg_source_1d(source, att_c * abs(flux))
      !time_derivs(:, FL_iv_elec) = time_derivs(:, FL_iv_elec) - source
      !time_derivs(:, FL_iv_nion) = time_derivs(:, FL_iv_nion) + source

      ! Detachment process
      if (FL_use_detach) then
         source(1) = 0
         source(n_cc) = 0
         source(2:n_cc-1) = (det_c(1:n_cc-2) + det_c(2:n_cc-1)) * 0.5_dp * FL_vars(2:n_cc-1, FL_iv_nion)
         time_derivs(:, FL_iv_elec) = time_derivs(:, FL_iv_elec) + source
         time_derivs(:, FL_iv_nion) = time_derivs(:, FL_iv_nion) - source
      end if

      call add_grad_flux_1d(time_derivs(:, FL_iv_elec), flux * inv_delta_x)

      if (FL_use_en) then
         ! ~~~ Energy source ~~~
         source(n_cc)     = 0.0_dp
         source(1:n_cc-1) = -0.5_dp * (fld * flux + en_loss * vars(1:n_cc-1, FL_iv_elec))
         source(2:n_cc)   = source(2:n_cc) - 0.5_dp * (fld * flux + en_loss * &
              vars(2:n_cc, FL_iv_elec))

         ! ~~~ energy transport ~~~
         call FL_transport_scheme(vars(:, FL_iv_en), -five_third * mob_c * fld, &
              five_third * dif_c * inv_delta_x, flux)

         ! Set time derivatives
         time_derivs(:, FL_iv_en) = source
         call add_grad_flux_1d(time_derivs(:, FL_iv_en), flux * inv_delta_x)
      end if
      
      flux = 0.0_dp
      ! ~~~ Positive ion transport ~~~
      call FL_transport_scheme(vars(:, FL_iv_ion), FL_pos_ion_mob * fld, FL_pos_ion_diff * inv_delta_x * ones, flux)
      call add_grad_flux_1d(time_derivs(:, FL_iv_ion), flux * inv_delta_x)

      ! Dirichlet
      ! time_derivs(1, :) = 0.0_dp
      ! time_derivs(n_cc, :) = 0.0_dp
      ! Neumann
      time_derivs(1, :) = time_derivs(2, :)
      time_derivs(n_cc, :) = time_derivs(n_cc-1, :)
   end subroutine FL_time_derivs

   subroutine set_stagg_source_1d(dens_c, src_f)
      use m_init_cond_1d
      real(dp), intent(inout) :: dens_c(:)
      real(dp), intent(in)    :: src_f(:)
      integer                 :: n_cc, iz_d

      iz_d        = int(INIT_DI/PD_dx+1) 
      n_cc             = size(dens_c)
      dens_c(n_cc)     = 0
      dens_c(1:iz_d-1) = 0.5_dp * src_f(1:iz_d-1)
      dens_c(2:iz_d)   = dens_c(2:iz_d) + dens_c(1:iz_d-1)
   end subroutine set_stagg_source_1d

   subroutine add_grad_flux_1d(dens_c, flux_f)
      real(dp), intent(inout) :: dens_c(:)
      real(dp), intent(in)    :: flux_f(:)
      integer                 :: n_cc

      n_cc             = size(dens_c)
      dens_c(1:n_cc-1) = dens_c(1:n_cc-1) - flux_f
      dens_c(2:n_cc)   = dens_c(2:n_cc) + flux_f
   end subroutine add_grad_flux_1d
   
   subroutine add_s_charge(fld, dt)
     use m_units_constants
     use m_init_cond_1d
     real(dp), intent(in)                 :: dt, fld(PD_grid_size-1)
     real(dp), dimension(PD_grid_size-1)  :: fluxi, fluxe, ones, mob_c, dif_c
     integer                              :: iz_d
     real(dp)                             :: inv_delta_x
     type(LT_loc_t)                       :: fld_locs(PD_grid_size-1)
       
     
     ones        = 1.0_dp
     iz_d        = int(INIT_DI/PD_dx+1) 
     inv_delta_x = 1.0_dp / PD_dx 
     
     fld_locs = LT_get_loc(FL_lkp_fld, abs(fld))
     mob_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_mob, fld_locs)
     dif_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_dif, fld_locs)

     call FL_transport_scheme(FL_vars(:, FL_iv_ion), FL_pos_ion_mob * fld, FL_pos_ion_diff * inv_delta_x * ones, fluxi)
     call FL_transport_scheme(FL_vars(:, FL_iv_elec), -mob_c * fld, dif_c * inv_delta_x, fluxe)
     FL_surface_charge = FL_surface_charge + (max(fluxi(iz_d), 0.0_dp) + max(fluxe(iz_d), 0.0_dp)+ &
                              sum(photo(iz_d+1:PD_grid_size)) * photo_yield * PD_dx) * dt * UC_elem_charge
                    
     
   end subroutine add_s_charge
   
   subroutine global_result(E, ne, imp, time)
      real(dp), intent(in) :: E, ne, imp, time
      integer              :: i

      if(E > FL_max_E) then 
        FL_max_E = E
        FL_E_t_max = 0.0_dp
      end if

      if(time > 1e-10_dp .and. E < FL_max_E .and. FL_E_t_max == 0.0_dp) FL_E_t_max = time 
      
      if(ne > FL_max_ne) then
        FL_max_ne = ne
        FL_ne_t_max = 0.0_dp
      end if
      if(time > 1e-10_dp .and. ne < FL_max_ne .and. FL_ne_t_max == 0.0_dp) FL_ne_t_max = time
      
      if(imp > FL_max_imp) then 
        FL_max_imp = imp
        FL_imp_t_max = 0.0_dp
      end if
      if(time > 1e-10_dp .and. imp < FL_max_imp .and. FL_imp_t_max == 0.0_dp) FL_imp_t_max = time
      
      do i = 1, 15 
        if(time >= 2e-9_dp*(16-i) .and. time < 2.0015e-9_dp*(16-i)) then
          write(*,*) FL_max_E, FL_max_ne, FL_max_imp, FL_E_t_max, FL_ne_t_max, FL_imp_t_max
          exit
        end if
      end do
   end subroutine global_result

   subroutine FL_advance(time, dt)
      use m_time_steppers
      use m_efield_1d
      use m_units_constants
      use m_init_cond_1d
      real(dp), intent(inout) :: dt
      real(dp), intent(inout) :: time
      integer                 :: iz_d
      real(dp)                :: fld(PD_grid_size-1), mob_c(PD_grid_size-1), source(PD_grid_size)
      type(LT_loc_t)          :: fld_locs(PD_grid_size-1), fld_locsc(PD_grid_size)
      real(dp), dimension(PD_grid_size) :: vel, alpha, fldc
        
      iz_d = int(INIT_DI/PD_dx+1) 

      call STEP_expl_trap_2d(FL_vars, time, dt, FL_time_derivs)
      source = (FL_vars(:, FL_iv_ion) - FL_vars(:, FL_iv_elec) - FL_vars(:, FL_iv_nion)) * UC_elem_charge 
      source(iz_d+1:PD_grid_size) = 0.0_dp
      call EF_compute_and_get_st(source, fld, FL_surface_charge)
      fld_locs = LT_get_loc(FL_lkp_fld, abs(fld))
      call add_s_charge(fld, dt)
      if (FL_use_en) then
         where (FL_vars(:, FL_iv_en) < 0) FL_vars(:, FL_iv_en) = 0
      end if
      where (FL_vars(:, FL_iv_elec) < 0) FL_vars(:, FL_iv_elec) = 0
        
      FL_vars(iz_d+1:PD_grid_size, FL_iv_elec) = 0.0_dp
      FL_vars(iz_d+1:PD_grid_size, FL_iv_ion) = 0.0_dp
      
      
        
        
      call EF_compute_and_get_st(source, fld, FL_surface_charge)
      fld_locs = LT_get_loc(FL_lkp_fld, abs(fld))
      mob_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_mob, fld_locs)
      dt = 0.7_dp * min(PD_dx / (epsilon(1.0_dp) + maxval(fld*mob_c)), &
           UC_eps0/(UC_elem_charge * maxval(mob_c * 0.5 * (FL_vars(1:PD_grid_size-1, FL_iv_elec) + &
           FL_vars(2:PD_grid_size, FL_iv_elec)))))

      if(dt > 1e-13) dt = 1e-13
      
      call EF_compute_and_get((FL_vars(:, FL_iv_ion) - FL_vars(:, FL_iv_elec) - FL_vars(:, FL_iv_nion)) &
           * UC_elem_charge, fldc, FL_surface_charge)   
      fld_locsc = LT_get_loc(FL_lkp_fld, abs(fldc))
      vel = LT_get_col_at_loc(FL_lkp_fld, FL_if_mob, fld_locsc) * abs(fldc)
      alpha = LT_get_col_at_loc(FL_lkp_fld, FL_if_src, fld_locsc)     
      
      
      call global_result(maxval(fld), maxval(FL_vars(:, FL_iv_elec)), maxval(vel*alpha * FL_vars(:, FL_iv_elec)), time) 
        
   end subroutine FL_advance

   subroutine FL_get_output(pos_data, sca_data, data_names, n_pos, n_sca, time, head_density)
      use m_efield_1d
      use m_units_constants
      real(dp), intent(out), allocatable                :: pos_data(:,:), sca_data(:)
      real(dp), intent(in)                              :: time, head_density
      character(len=*), intent(out), allocatable        :: data_names(:)
      integer, intent(out)                              :: n_pos, n_sca
      integer                                           :: n, ix
      real(dp) , dimension(PD_grid_size)                :: temp_data, fld_en, src_c, att_c, mob_c
      type(LT_loc_t)                                    :: fld_locs(PD_grid_size)

      n_pos = 7
      if(FL_use_en) n_pos = n_pos + 2
      n_sca = 2
      allocate(pos_data(PD_grid_size, n_pos))
      allocate(sca_data(n_sca))
      allocate(data_names(n_pos+n_sca))

      do n = 1, PD_grid_size
         temp_data(n) = (n-1) * PD_dx
      end do

      data_names(n_sca+1) = "position (m)"
      pos_data(:,1) = temp_data


      call EF_compute_and_get((FL_vars(:, FL_iv_ion) - FL_vars(:, FL_iv_elec) - FL_vars(:, FL_iv_nion)) &
           * UC_elem_charge, temp_data, FL_surface_charge)
      fld_en = LT_get_col(FL_lkp_fld, FL_if_en, abs(temp_data))
      fld_locs = LT_get_loc(FL_lkp_fld, abs(temp_data))
      src_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_src, fld_locs)
      att_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_att, fld_locs)
      mob_c = LT_get_col_at_loc(FL_lkp_fld, FL_if_mob, fld_locs)

      data_names(n_sca+2) = "electric field (V/m)"
      pos_data(:,2) = temp_data

      data_names(n_sca+3) = "electron density (1/m3)"
      pos_data(:,3) = FL_vars(:, FL_iv_elec)

      ! Set sca data ********
      data_names(1) = "time"
      sca_data(1) = time

      ! Find where the electron density first exceeds head_density and store as head_pos
      !data_names(2) = "head_pos"

      !ix = -1
      !if (pos_data(1,2) > 0) then ! Check sign of electric field
      !   do n = 2, PD_grid_size
      !      if (pos_data(n, 3) > head_density) then
      !         ix = n
      !         exit
      !      end if
      !   end do
      !else
      !   do n = PD_grid_size, 2, -1
      !      if (pos_data(n-1, 3) > head_density) then
      !         ix = n
      !         exit
      !      end if
      !   end do
      !end if
      

      !if (ix == -1) then
      !   sca_data(2) = 0
      !else
         ! Interpolate between points ix-1 and ix
      !   write(*,*) pos_data(ix-1, 3), ix
      !   sca_data(2) = (head_density - pos_data(ix-1, 3)) / (pos_data(ix, 3) - pos_data(ix-1, 3))
      !   sca_data(2) = (ix - 2 + sca_data(2)) * PD_dx
      !end if

      data_names(n_sca+4) = "pos ion density (1/m3)"
      pos_data(:,4) = FL_vars(:, FL_iv_ion)

      data_names(n_sca+5) = "net charge density (1/m3)"
      pos_data(:,5) = (FL_vars(:, FL_iv_ion) - FL_vars(:, FL_iv_elec) - FL_vars(:, FL_iv_nion)) &
           * UC_elem_charge
      
      data_names(n_sca+6) = "photo-elec density (1/(s m3))"
      pos_data(:,6) = photo(:)
      
      data_names(n_sca+7) = "impact ion (1/(s m3))"
      pos_data(:,7) = FL_vars(:, FL_iv_elec) * abs(temp_data) * (src_c-att_c) * mob_c   


      if (FL_use_en) then
        data_names(n_sca+8) = "energy density (1/m3)"
        data_names(n_sca+9) = "mean energy (eV/m3)"
        pos_data(:,8) = FL_vars(:, FL_iv_en)
        pos_data(:,9) = (FL_vars(:, FL_iv_en) + FL_small_dens * fld_en) &
              / (FL_small_dens + FL_vars(:, FL_iv_elec))
      !else
      !   pos_data(:,8) = FL_vars(:, FL_iv_elec) * fld_en
      !   pos_data(:,9) = fld_en
      end if

   end subroutine FL_get_output

   subroutine FL_get_coeffs(coeff_data, coeff_names, n_coeffs)
      real(dp), intent(out), allocatable :: coeff_data(:,:)
      character(len=*), intent(out), allocatable :: coeff_names(:)
      integer, intent(out) :: n_coeffs
      integer :: ix, n_rows, n_fld_coeffs, n_en_coeffs

      if (FL_use_en) then
         n_en_coeffs = LT_get_num_cols(FL_lkp_en) + 1
      else
         n_en_coeffs = 0
      end if
      n_fld_coeffs = LT_get_num_cols(FL_lkp_fld) + 1

      n_coeffs = n_fld_coeffs + n_en_coeffs
      n_rows = LT_get_num_rows(FL_lkp_fld)
      allocate(coeff_data(n_coeffs, n_rows))
      allocate(coeff_names(n_coeffs))

      call LT_get_data(FL_lkp_fld, coeff_data(1, :), coeff_data(2:n_fld_coeffs, :))
      coeff_names(1) = "efield (V/m)"
      coeff_names(1+FL_if_en) = "fld_en (eV)"
      if (FL_use_detach) coeff_names(1+FL_if_det) = "fld_det (1/s)"
      if (.not. FL_use_en_mob) coeff_names(1+FL_if_mob) = "fld_mob (m2/Vs)"
      if (.not. FL_use_en_dif) coeff_names(1+FL_if_dif) = "fld_dif (m2/s)"
      if (.not. FL_use_en_src) then
         coeff_names(1+FL_if_src) = "fld_src (1/m)"
         coeff_names(1+FL_if_att) = "fld_att (1/m)"
      end if

      if (FL_use_en) then
         ix = n_fld_coeffs + 1
         call LT_get_data(FL_lkp_en, coeff_data(ix, :), coeff_data(ix+1:, :))
         coeff_names(ix) = "energy (eV/s)"
         coeff_names(ix+FL_ie_loss) = "en_loss (eV/s)"
         if (FL_use_en_mob) coeff_names(ix+FL_ie_mob) = "en_mob (m2/Vs)"
         if (FL_use_en_dif) coeff_names(ix+FL_ie_dif) = "en_dif (m2/s)"
         if (FL_use_en_src) then
            coeff_names(ix+FL_ie_src) = "en_src (1/m)"
            coeff_names(ix+FL_ie_att) = "en_att (1/m)"
         end if
      end if

   end subroutine FL_get_coeffs

   !> Initialize photoionization parameters
   subroutine photoi_initialize(cfg, gas_pressure)
      use m_photoi_1d
      use m_config
      type(CFG_t), intent(inout) :: cfg 
      real(dp), intent(in)           :: gas_pressure
      real(dp) :: photoi_eta = 0.1_dp

      call CFG_get(cfg, "photo_yield", photo_yield)


      call phmc_initialize(cfg, gas_pressure)

   end subroutine photoi_initialize

   !> Sets the photoionization
   subroutine photoi_set_src(dt, gas_pressure, photoi_eta)
     use m_photoi_1d
     use m_random
     use m_efield_1d
     use m_units_constants
     real(dp), intent(in)              :: dt
     real(dp), intent(in)              :: gas_pressure, photoi_eta
     real(dp), parameter               :: p_quench = 40.0e-3_dp
     real(dp)                          :: quench_fac, src(PD_grid_size)
     real(dp), dimension(PD_grid_size) :: vel, alpha, fld
     real(dp)                          :: source(PD_grid_size)
     type(LT_loc_t)                    :: fld_locs(PD_grid_size)
       
     call EF_compute_and_get((FL_vars(:, FL_iv_ion) - FL_vars(:, FL_iv_elec) - FL_vars(:, FL_iv_nion)) &
          * UC_elem_charge, fld, FL_surface_charge)   
     fld_locs = LT_get_loc(FL_lkp_fld, abs(fld))
     vel = LT_get_col_at_loc(FL_lkp_fld, FL_if_mob, fld_locs) * abs(fld)
     alpha = LT_get_col_at_loc(FL_lkp_fld, FL_if_src, fld_locs)     

     quench_fac = p_quench / (gas_pressure + p_quench)
     rng_int4_seed = get_random_seed()

     src = vel * alpha * FL_vars(:, FL_iv_elec) * quench_fac * 0.1_dp
     where (src < 0) src = 0
     
     call phmc_set_src_1d(rng, src, photo, dt)

   end subroutine photoi_set_src
   
   !> Get a random seed based on the current time
   function get_random_seed() result(seed)
     integer :: seed(4)
     integer :: time, i

     call system_clock(time)
     do i = 1, 4
        seed(i) = ishftc(time, i*8)
     end do
   end function get_random_seed
   

end module m_fluid_dd_1d
