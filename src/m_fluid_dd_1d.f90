!> A not so accurate implementation of a drift-diffusion-reaction plasma fluid model

module m_fluid_dd_1d
   use m_lookup_table
   use m_transport_schemes
   use m_phys_domain

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
   real(dp)              :: FL_small_dens
   real(dp)              :: FL_max_energy

   type(LT_table_t) :: FL_lkp_fld, FL_lkp_en
   procedure(TS_dd_1d_type), pointer :: FL_transport_scheme

   public :: FL_init_cfg
   public :: FL_advance
   public :: FL_get_output
   public :: FL_get_coeffs
   public :: FL_max_edens_at_boundary

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
      do n = 1, PD_grid_size
         xx = (n-1) * PD_dx
         call INIT_get_elec_dens(xx, FL_vars(n, FL_iv_elec))
         call INIT_get_ion_dens(xx, FL_vars(n, FL_iv_ion))
         call INIT_get_nion_dens(xx, FL_vars(n, FL_iv_nion))
         if (FL_use_en) then
            FL_vars(n, FL_iv_en) = FL_vars(n, FL_iv_elec) * &
                 LT_get_col(FL_lkp_fld, FL_if_en, EF_get_at_pos((n-1)*PD_dx))
         end if
      end do

   end subroutine FL_init_cfg

   real(dp) function FL_max_edens_at_boundary()
      FL_max_edens_at_boundary = &
           max(FL_vars(1, FL_iv_elec), FL_vars(PD_grid_size, FL_iv_elec))
   end function FL_max_edens_at_boundary

   subroutine FL_time_derivs(vars, time, time_derivs)
      use m_efield_1d
      use m_units_constants
      real(dp), intent(in)        :: vars(:,:), time
      real(dp), intent(out)       :: time_derivs(:,:)

      integer                     :: n_cc
      real(dp), parameter         :: five_third = 5 / 3.0_dp
      real(dp)                    :: inv_delta_x
      real(dp), dimension(size(vars,1)-1) :: &
           mob_c, dif_c, src_c, att_c, det_c, flux, fld, fld_en, en_loss, mean_en
      real(dp)                    :: source(size(vars,1))
      type(LT_loc_t) :: fld_locs(size(vars,1)-1), en_locs(size(vars,1)-1)

      n_cc = size(vars, 1)
      inv_delta_x = 1.0_dp / PD_dx

      ! Get electric field
      source = (vars(:, FL_iv_ion) - vars(:, FL_iv_elec) - vars(:, FL_iv_nion)) * UC_elem_charge
      call EF_compute_and_get_st(source, fld)

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
      call set_stagg_source_1d(source, src_c * abs(flux))
      time_derivs(:, FL_iv_elec) = source
      time_derivs(:, FL_iv_ion) = source
      time_derivs(:, FL_iv_nion) = 0

      ! ~~~ Attachment source ~~~ TODO: try different formulations
      call set_stagg_source_1d(source, att_c * abs(flux))
      time_derivs(:, FL_iv_elec) = time_derivs(:, FL_iv_elec) - source
      time_derivs(:, FL_iv_nion) = time_derivs(:, FL_iv_nion) + source

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

      ! Dirichlet
      ! time_derivs(1, :) = 0.0_dp
      ! time_derivs(n_cc, :) = 0.0_dp
      ! Neumann
      time_derivs(1, :) = time_derivs(2, :)
      time_derivs(n_cc, :) = time_derivs(n_cc-1, :)
   end subroutine FL_time_derivs

   subroutine set_stagg_source_1d(dens_c, src_f)
      real(dp), intent(inout) :: dens_c(:)
      real(dp), intent(in)    :: src_f(:)
      integer                 :: n_cc

      n_cc             = size(dens_c)
      dens_c(n_cc)     = 0
      dens_c(1:n_cc-1) = 0.5_dp * src_f
      dens_c(2:n_cc)   = dens_c(2:n_cc) + dens_c(1:n_cc-1)
   end subroutine set_stagg_source_1d

   subroutine add_grad_flux_1d(dens_c, flux_f)
      real(dp), intent(inout) :: dens_c(:)
      real(dp), intent(in)    :: flux_f(:)
      integer                 :: n_cc

      n_cc             = size(dens_c)
      dens_c(1:n_cc-1) = dens_c(1:n_cc-1) - flux_f
      dens_c(2:n_cc)   = dens_c(2:n_cc) + flux_f
   end subroutine add_grad_flux_1d

   subroutine FL_advance(time, dt)
      use m_time_steppers
      real(dp), intent(in)    :: dt
      real(dp), intent(inout) :: time

      call STEP_expl_trap_2d(FL_vars, time, dt, FL_time_derivs)
      where (FL_vars(:, FL_iv_en) < 0) FL_vars(:, FL_iv_en) = 0
      where (FL_vars(:, FL_iv_elec) < 0) FL_vars(:, FL_iv_elec) = 0
   end subroutine FL_advance

   subroutine FL_get_output(pos_data, sca_data, data_names, n_pos, n_sca, time, head_density)
      use m_efield_1d
      use m_units_constants
      real(dp), intent(out), allocatable                :: pos_data(:,:), sca_data(:)
      real(dp), intent(in)                              :: time, head_density
      character(len=*), intent(out), allocatable :: data_names(:)
      integer, intent(out)                              :: n_pos, n_sca
      integer                                           :: n, ix
      real(dp)                                          :: temp_data(PD_grid_size), fld_en(PD_grid_size)

      n_pos = 7
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
           * UC_elem_charge, temp_data)
      fld_en = LT_get_col(FL_lkp_fld, FL_if_en, abs(temp_data))

      data_names(n_sca+2) = "electric field (V/m)"
      pos_data(:,2) = temp_data

      data_names(n_sca+3) = "electron density (1/m3)"
      pos_data(:,3) = FL_vars(:, FL_iv_elec)

      ! Set sca data ********
      data_names(1) = "time"
      sca_data(1) = time

      ! Find where the electron density first exceeds head_density and store as head_pos
      data_names(2) = "head_pos"

      ix = -1
      if (pos_data(1,2) > 0) then ! Check sign of electric field
         do n = 2, PD_grid_size
            if (pos_data(n, 3) > head_density) then
               ix = n
               exit
            end if
         end do
      else
         do n = PD_grid_size, 2, -1
            if (pos_data(n-1, 3) > head_density) then
               ix = n
               exit
            end if
         end do
      end if

      if (ix == -1) then
         sca_data(2) = 0
      else
         ! Interpolate between points ix-1 and ix
         sca_data(2) = (head_density - pos_data(ix-1, 3)) / (pos_data(ix, 3) - pos_data(ix-1, 3))
         sca_data(2) = (ix - 2 + sca_data(2)) * PD_dx
      end if

      data_names(n_sca+4) = "ion density (1/m3)"
      pos_data(:,4) = FL_vars(:, FL_iv_ion)

      data_names(n_sca+5) = "nion density (1/m3)"
      pos_data(:,5) = FL_vars(:, FL_iv_nion)

      data_names(n_sca+6) = "energy density (1/m3)"
      data_names(n_sca+7) = "mean energy (eV/m3)"

      if (FL_use_en) then
         pos_data(:,6) = FL_vars(:, FL_iv_en)
         pos_data(:,7) = (FL_vars(:, FL_iv_en) + FL_small_dens * fld_en) &
              / (FL_small_dens + FL_vars(:, FL_iv_elec))
      else
         pos_data(:,6) = FL_vars(:, FL_iv_elec) * fld_en
         pos_data(:,7) = fld_en
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

end module m_fluid_dd_1d
