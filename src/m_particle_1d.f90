module m_particle_1d
   use m_particle_core

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)
   integer, parameter    :: PM_num_vars = 3
   integer, parameter    :: PM_iv_elec  = 1, PM_iv_ion = 2, PM_iv_en = 3

   integer               :: PM_grid_size
   real(dp), allocatable :: PM_vars(:,:)
   real(dp)              :: PM_delta_x, PM_inv_delta_x
   real(dp)              :: PM_transverse_area
   real(dp)              :: PM_grid_volume, PM_inv_grid_volume
   real(dp)              :: PM_part_per_cell, PM_vel_rel_weight

   ! Public routines
   public :: PM_initialize
   public :: PM_advance
   public :: PM_get_output
   public :: PM_get_eedf
   public :: PC_get_num_sim_part
   public :: PC_get_num_real_part
   public :: PC_get_coeffs
   public :: PM_max_edens_at_boundary

contains

   subroutine PM_initialize(cross_secs, n_part_init)
      use m_config
      use m_random
      use m_init_cond_1d
      use m_units_constants
      use m_cross_sec
      type(CS_type), intent(in) :: cross_secs(:)
      integer, intent(in)       :: n_part_init

      integer  :: n, pp, num_part
      real(dp) :: xx, sum_elec_dens
      real(dp) :: elec_dens, ion_dens, en_dens, mean_en
      real(dp) :: pos(3), vel(3), accel(3)

      PM_grid_size       = CFG_get_int("grid_num_points")
      PM_delta_x         = CFG_get_real("grid_delta_x")
      PM_inv_delta_x     = 1 / PM_delta_x

      sum_elec_dens = 0.0_dp
      do n = 1, PM_grid_size
         call INIT_get_elec_dens((n-1) * PM_delta_x, elec_dens)
         sum_elec_dens = sum_elec_dens + elec_dens
      end do

      if (size(cross_secs) < 1) then
         print *, "No cross sections given to particle module!"
         stop
      end if

      PM_grid_volume     = n_part_init / sum_elec_dens
      PM_transverse_area = PM_grid_volume / PM_delta_x
      PM_inv_grid_volume = 1 / PM_grid_volume
      PM_part_per_cell   = CFG_get_real("apm_part_per_cell")
      PM_vel_rel_weight  = CFG_get_real("apm_vel_rel_weight")

      allocate(PM_vars(PM_grid_size, PM_num_vars))
      PM_vars = 0.0_dp

      call PC_initialize(UC_elec_mass, cross_secs, &
           CFG_get_int("part_lkptbl_size"), CFG_get_real("part_max_energy_eV"))
      call PC_set_pptrs(pptr_coll_callback = coll_callback)

      ! Initialization of electron and ion density
      do n = 1, PM_grid_size
         xx = (n-1) * PM_delta_x
         call INIT_get_elec_dens(xx, elec_dens)
         call INIT_get_ion_dens(xx, ion_dens)
         call INIT_get_en_dens(xx, en_dens)
         num_part = nint(elec_dens * PM_grid_volume)

         if (elec_dens > 0.0_dp) then
            mean_en = en_dens * UC_elec_volt / elec_dens
         else
            mean_en = 0.0_dp
         end if

         pos = (/0.0_dp, 0.0_dp, xx/)

         ! Accel is set after computing electric field
         accel = 0.0_dp

         do pp = 1, num_part
            ! Set a Maxwellian velocity distribution with the correct mean energy
            vel = (/RNG_normal(), RNG_normal(), RNG_normal()/)
            vel = vel * sqrt(2 * mean_en / (3 * UC_elec_mass))

            call PC_create_part(pos, vel, accel, 1.0_dp, 0.0_dp)
         end do

         num_part = nint(ion_dens * PM_grid_volume)
         pm_vars(n, PM_iv_ion) = num_part / PM_grid_volume
      end do

      call set_elec_density()
      call PM_update_efield()
      call PC_correct_new_accel(0.0_dp, PM_set_accel)

   end subroutine PM_initialize

   real(dp) function PM_max_edens_at_boundary()
      PM_max_edens_at_boundary = &
           max(PM_vars(1, PM_iv_elec), PM_vars(PM_grid_size, PM_iv_elec))
   end function PM_max_edens_at_boundary

   subroutine set_elec_density()
      PM_vars(:, PM_iv_elec) = 0.0_dp
      call PC_loop_ipart(add_elec_to_dens)
   end subroutine set_elec_density

   subroutine set_elec_en_density()
      PM_vars(:, PM_iv_en) = 0.0_dp
      call PC_loop_ipart(add_elec_en_to_dens)
   end subroutine set_elec_en_density

   subroutine add_to_dens_cic(amount, zz, dens)
      real(dp), intent(in)    :: amount, zz
      real(dp), intent(inout) :: dens(:)

      integer                 :: low_ix
      real(dp)                :: temp, weight, dens_dif

      ! Index i is at (i-1) * dx
      temp   = zz * PM_inv_delta_x
      low_ix = floor(temp) + 1
      weight = low_ix - temp

      dens_dif = amount * PM_inv_grid_volume

      if (low_ix < 1 .or. low_ix > PM_grid_size-1) then
         print *, "Particle module: a particle is outside the computational domain"
         stop
      end if
      dens(low_ix)   = dens(low_ix) + weight * dens_dif
      dens(low_ix+1) = dens(low_ix+1) + (1-weight) * dens_dif
   end subroutine add_to_dens_cic

   real(dp) function get_from_dens_lin_interp(zz, dens)
      real(dp), intent(in)    :: zz
      real(dp), intent(inout) :: dens(:)

      integer                 :: low_ix
      real(dp)                :: temp, weight

      ! Index i is at (i-1) * dx
      temp   = zz * PM_inv_delta_x
      low_ix = floor(temp) + 1
      weight = low_ix - temp
      get_from_dens_lin_interp = dens(low_ix) * weight + dens(low_ix+1) * (1-weight)
   end function get_from_dens_lin_interp

   subroutine add_elec_to_dens(my_part)
      type(PC_part_t), intent(in) :: my_part
      call add_to_dens_cic(my_part%weight, my_part%x(3), PM_vars(:, PM_iv_elec))
   end subroutine add_elec_to_dens

   subroutine add_elec_en_to_dens(my_part)
      use m_units_constants
      type(PC_part_t), intent(in) :: my_part
      real(dp) :: energy
      energy = 0.5_dp * my_part%weight * UC_elec_mass * sum(my_part%v**2)
      call add_to_dens_cic(energy, my_part%x(3), PM_vars(:, PM_iv_en))
   end subroutine add_elec_en_to_dens

   subroutine coll_callback(my_part, col_type, col_ix)
      use m_cross_sec
      type(PC_part_t), intent(in) :: my_part
      integer, intent(in) :: col_type, col_ix

      select case(col_type)
      case (CS_ionize_t)
         call add_to_dens_cic(my_part%weight, my_part%x(3), PM_vars(:, PM_iv_ion))
      case (CS_attach_t)
         call add_to_dens_cic(-my_part%weight, my_part%x(3), PM_vars(:, PM_iv_ion))
      end select
   end subroutine coll_callback

   subroutine PM_update_efield()
      use m_efield_1d
      use m_units_constants
      call EF_compute((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge)
   end subroutine PM_update_efield

   subroutine PM_set_accel(my_part, accel)
      use m_efield_1d
      use m_units_constants
      type(PC_part_t), intent(in) :: my_part
      real(dp), intent(out)       :: accel(3)
      real(dp), parameter         :: accel_fac = UC_elec_charge / UC_elec_mass

      ! Only acceleration in the z-direction
      accel(1:2) = 0.0_dp
      accel(3) = accel_fac * EF_get_at_pos(my_part%x(3))
   end subroutine PM_set_accel

   subroutine PM_advance(dt, do_apm)
      real(dp), intent(in) :: dt
      logical, intent(in) :: do_apm

      call PC_advance(dt)
      call set_elec_density()
      call PM_update_efield()
      call PC_correct_new_accel(dt, PM_set_accel)
      if (do_apm) call PM_adapt_weights()
   end subroutine PM_advance

   subroutine PM_adapt_weights()
      real(dp) :: coord_weights(6), n_start, n_end
      n_start = PC_get_num_real_part()
      coord_weights(1:3) = (/0.0_dp, 0.0_dp, 1.0_dp/) ! xyz coords
      coord_weights(4:6) = PM_vel_rel_weight * PM_delta_x / PM_part_per_cell
      call PC_merge_and_split(coord_weights, PM_delta_x, &
           get_desired_weight, PC_merge_part_rxv, PC_split_part)

      call set_elec_density()
      call PM_update_efield()
      call PC_set_accel(PM_set_accel)
      n_end = PC_get_num_real_part()

      ! print *, (n_end - n_start) / n_start
      write(*,'(A, E12.4, A, I0)') "After merging real/sim part: ", PC_get_num_real_part(), " / ", PC_get_num_sim_part()
   end subroutine PM_adapt_weights

   real(dp) function get_desired_weight(my_part)
      type(PC_part_t), intent(in) :: my_part
      real(dp)                    :: elec_dens
      elec_dens = get_from_dens_lin_interp(my_part%x(3), PM_vars(:, PM_iv_elec))
      get_desired_weight = elec_dens * PM_grid_volume / PM_part_per_cell
      get_desired_weight = max(1.0_dp, get_desired_weight)
   end function get_desired_weight

   subroutine PM_get_output(pos_data, sca_data, data_names, n_pos, n_sca, time, head_density)
      use m_efield_1d
      use m_units_constants
      real(dp), intent(out), allocatable                :: pos_data(:,:), sca_data(:)
      real(dp), intent(in) :: time, head_density
      character(len=*), intent(out), allocatable :: data_names(:)
      integer, intent(out)                              :: n_pos, n_sca
      integer                                           :: n, ix
      real(dp)                                          :: temp_data(PM_grid_size)

      n_pos = 6
      n_sca = 2
      allocate(pos_data(PM_grid_size, n_pos))
      allocate(sca_data(n_sca))
      allocate(data_names(n_pos+n_sca))

      do n = 1, PM_grid_size
         temp_data(n) = (n-1) * PM_delta_x
      end do

      data_names(n_sca+1) = "position (m)"
      pos_data(:,1) = temp_data

      call EF_compute_and_get((PM_vars(:, PM_iv_ion) - PM_vars(:, PM_iv_elec)) * UC_elem_charge, temp_data)
      data_names(n_sca+2) = "electric field (V/m)"
      pos_data(:,2) = temp_data

      call set_elec_density()
      data_names(n_sca+3) = "electron density (1/m3)"
      pos_data(:,3) = PM_vars(:, PM_iv_elec)

      ! Set sca data ********
      data_names(1) = "time"
      sca_data(1) = time

      ! Find where the electron density first exceeds head_density and store as head_pos
      data_names(2) = "head_pos"

      ix = -1
      if (pos_data(1,2) > 0) then ! Check sign of electric field
         do n = 2, PM_grid_size
            if (pos_data(n, 3) > head_density) then
               ix = n
               exit
            end if
         end do
      else
         do n = PM_grid_size, 2, -1
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
         sca_data(2) = (ix - 2 + sca_data(2)) * PM_delta_x
      end if

      data_names(n_sca+4) = "ion density (1/m3)"
      pos_data(:,4) = PM_vars(:, PM_iv_ion)

      call set_elec_en_density()
      data_names(n_sca+5) = "energy density (eV/m3)"
      pos_data(:,5) = PM_vars(:, PM_iv_en) / UC_elec_volt

      data_names(n_sca+6) = "mean energy density (1/m3)"
      where (PM_vars(:, PM_iv_elec) > 0)
         pos_data(:,6) = PM_vars(:, PM_iv_en) / (UC_elec_volt * PM_vars(:, PM_iv_elec))
      elsewhere
         pos_data(:,6) = 0.0_dp
      end where
   end subroutine PM_get_output

   subroutine PM_get_eedf(eedf, energy_range, field_range, n_bins)
      use m_units_constants
      real(dp), intent(inout), allocatable :: eedf(:,:)
      real(dp), intent(in) :: energy_range(2), field_range(2)
      integer, intent(in) :: n_bins

      integer :: ix
      real(dp) :: accel_range(2)
      allocate(eedf(2, n_bins))
      do ix = 1, n_bins
         eedf(1, ix) = energy_range(1) + (ix-0.5_dp) * (energy_range(2) - energy_range(1)) / n_bins
      end do

      ! Convert to range of velocity**2
      eedf(1,:) = eedf(1,:) * UC_elec_volt / (0.5_dp * PC_get_part_mass())
      accel_range = field_range * UC_elem_charge / UC_elec_mass

      call PC_get_histogram(get_vel_squared, select_part_by_accel, accel_range, eedf(1,:), eedf(2,:))

      ! Convert to energy range
      eedf(1,:) = eedf(1,:) * 0.5_dp * PC_get_part_mass() / UC_elec_volt
   end subroutine PM_get_eedf

   real(dp) function get_vel_squared(my_part)
      type(PC_part_t), intent(in) :: my_part
      get_vel_squared = sum(my_part%v**2)
   end function get_vel_squared

   logical function select_part_by_accel(my_part, accel_range)
      type(PC_part_t), intent(in) :: my_part
      real(dp), intent(in)        :: accel_range(:)
      real(dp)                    :: accel
      accel                = norm2(my_part%a)
      select_part_by_accel = accel >= accel_range(1) .and. accel <= accel_range(2)
   end function select_part_by_accel

end module m_particle_1d
