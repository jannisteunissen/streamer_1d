module m_init_cond_1d

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)
   integer, parameter :: INIT_gaussian_type = 1, &
        INIT_block_type = 2
   integer            :: INIT_type          = -1
   real(dp)           :: INIT_location
   real(dp)           :: INIT_width
   real(dp)           :: INIT_dens
   real(dp)           :: INIT_energy
   real(dp)           :: INIT_small_dens
   logical :: INIT_use_nion

   public :: INIT_init
   public :: INIT_get_elec_dens
   public :: INIT_get_ion_dens
   public :: INIT_get_en_dens
   public :: INIT_get_nion_dens

contains

   subroutine INIT_init()
      use m_config
      character(len=20) :: init_name
      real(dp)          :: domain_length

      ! Get the type of initial condition
      call CFG_get("init_cond_name", init_name)

      select case (init_name)
      case ("gaussian")
         INIT_type = INIT_gaussian_type
      case ("block")
         INIT_type = INIT_block_type
      case default
         print *, "INIT_init_with_config: unknown initial condition specified"
         stop
      end select

      domain_length = CFG_get_real("grid_delta_x") * (CFG_get_int("grid_num_points") - 1)
      INIT_width = CFG_get_real("init_width")
      INIT_location = CFG_get_real("init_rel_pos") * domain_length
      INIT_dens = CFG_get_real("init_dens")
      INIT_energy = CFG_get_real("init_elec_energy")
      INIT_use_nion = CFG_get_logic("init_use_neg_ion")
      INIT_small_dens = 0.0_dp !CFG_get_real("init_cutoff_density")
   end subroutine INIT_init

   real(dp) function get_dens(xx)
      real(dp), intent(in) :: xx

      select case (INIT_type)
      case (INIT_gaussian_type)
         get_dens = INIT_dens * exp(-(xx - INIT_location)**2 / (2 * INIT_width**2))
         if (get_dens < INIT_small_dens) get_dens = 0.0_dp
      case (INIT_block_type)
         if (abs(xx - INIT_location) < INIT_width) then
            get_dens = INIT_dens
         else
            get_dens = 0.0_dp
         end if
      case default
         get_dens = 0.0_dp
         print *, "m_initial_cond is being used without proper initialization"
         stop
      end select
   end function get_dens

   subroutine INIT_get_elec_dens(xx, elec_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: elec_dens
      if (.not. INIT_use_nion) then
         elec_dens = get_dens(xx)
      else
         elec_dens = 0
      end if
   end subroutine INIT_get_elec_dens

   subroutine INIT_get_nion_dens(xx, nion_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: nion_dens
      if (INIT_use_nion) then
         nion_dens = get_dens(xx)
      else
         nion_dens = 0
      end if
   end subroutine INIT_get_nion_dens

   subroutine INIT_get_ion_dens(xx, ion_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: ion_dens
      ion_dens = get_dens(xx)
   end subroutine INIT_get_ion_dens

   subroutine INIT_get_en_dens(xx, en_dens)
      real(dp), intent(in) :: xx
      real(dp), intent(out) :: en_dens
      real(dp) :: elec_dens
      call INIT_get_elec_dens(xx, elec_dens)
      en_dens = INIT_energy * elec_dens
   end subroutine INIT_get_en_dens

end module m_init_cond_1d
