!> Module that provides routines for operations on photons in 1D
module m_photoi_1d
  use m_lookup_table
  use m_phys_domain
  use m_units_constants
  use m_transport_data

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type phmc_tbl_t
     type(LT_table_t) :: tbl           !< The lookup table
     real(dp)         :: frac_in_tbl   !< Fraction photons in table
  end type phmc_tbl_t


  
  type(LT_table_t)       :: FL_lkp
  type(LT_loc_t)         :: loc
  character(len=100)     :: input_file


  ! Public methods
  public :: phmc_initialize
  public :: phmc_set_src_1d
  

contains

  !> Initialize photoionization parameters
  subroutine phmc_initialize(cfg, gas_pressure)
    use m_config
    type(CFG_t), intent(inout) :: cfg
    real(dp), intent(in)       :: gas_pressure
    real(dp), allocatable  :: x_data(:), y_data(:) 

    
    FL_lkp = LT_create(0.0_dp, 0.007_dp, 7000, 0)
    input_file = "input/air_transport_data_siglo.txt"

    call TD_get_td_from_file(input_file, "AIR", &
         trim("abs_func"), x_data, y_data)
    call LT_add_col(FL_lkp, x_data, y_data)
    
    
  end subroutine phmc_initialize



  !> Set the source term due to photoionization for 2D models. At most
  !> phmc_num_photons discrete photons are produced.
  subroutine phmc_set_src_1d(rng, src, photo, dt)
    use m_random
    use m_init_cond_1d
    use omp_lib

    type(RNG_t), intent(inout)  :: rng    !< Random number generator
    real(dp), intent(in)        :: src(:)
    real(dp), intent(inout)     :: photo(:)
    real(dp), intent(in)        :: dt
    integer                     :: i, n, n_create, n_used, i_ph, i_abs, iz_d, j, k
    real(dp)                    :: fac, dist, src_p(PD_grid_size)
    real(dp)                    :: sum_production
    real(dp), allocatable       :: xyz_src(:, :)
    real(dp), allocatable       :: xyz_abs(:, :)
    real(dp), allocatable       :: tmp(:) 
    real(dp), parameter         :: pi = acos(-1.0_dp)
    real(dp), parameter         :: dz = 1e-6_dp


    iz_d = int(INIT_DI/PD_dx+1) 
    photo = 0.0_dp
                    
    
    do i = 1, PD_grid_size
      do j = 1, iz_d
        loc = LT_get_loc(FL_lkp, abs(i-j)*PD_dx)
        photo(i) = photo(i) + src(j) * LT_get_col_at_loc(FL_lkp, 1, loc) * PD_dx
      end do
    end do
    


  end subroutine phmc_set_src_1d
    
end module m_photoi_1d
