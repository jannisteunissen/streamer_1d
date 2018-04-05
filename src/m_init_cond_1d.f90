module m_init_cond_1d
  use m_lookup_table
  use m_transport_data

  implicit none
  private

  integer, parameter  :: dp                 = kind(0.0d0)
  integer, parameter  :: INIT_gaussian_type = 1, &
       INIT_block_type                     = 2
  integer             :: INIT_type          = -1
  real(dp)            :: INIT_location
  real(dp), protected :: eps_DI
  real(dp), protected :: INIT_DI
  real(dp)            :: INIT_width
  real(dp)            :: INIT_dens
  real(dp)            :: INIT_energy
  real(dp)            :: INIT_small_dens
  integer             :: INIT_seed_sign, FL_if_elec, FL_if_ion 
  logical             :: INIT_use_nion
  real(dp)            :: INIT_back
  !character(len=100)  :: input_file
  !type(LT_table_t)    :: FL_lkp
  
  


  public :: INIT_init
  public :: INIT_get_elec_dens
  public :: INIT_get_ion_dens
  public :: INIT_get_en_dens
  public :: INIT_get_nion_dens
  public :: eps_DI
  public :: INIT_DI
  public :: INIT_get_back
  public :: FL_if_elec, FL_if_ion

contains

  subroutine INIT_init(cfg)
    use m_config
    use m_phys_domain

    type(CFG_t), intent(in) :: cfg
    character(len=20)       :: init_name
    real(dp), allocatable   :: x_data(:), y_data(:)                       
    
    !input_file = "input/air_transport_data_siglo.txt"
    !FL_lkp = LT_create(0.0_dp, 0.005_dp, 1000, 0)
    

    ! Get the type of initial condition
    call CFG_get(cfg, "init_cond_name", init_name)

    select case (init_name)
    case ("gaussian")
       INIT_type = INIT_gaussian_type
    case ("block")
       INIT_type = INIT_block_type
    case default
       stop "INIT_init_with_config: unknown initial condition specified"
    end select

    call CFG_get(cfg, "eps_DI", eps_DI)
    call CFG_get(cfg, "init_width", INIT_width)
    call CFG_get(cfg, "init_rel_pos", INIT_location)
    INIT_location = INIT_location * PD_length    
    call CFG_get(cfg, "init_DI_pos", INIT_DI)
    INIT_DI = INIT_DI * PD_length
    call CFG_get(cfg, "init_dens", INIT_dens)
    call CFG_get(cfg, "seed_sign", INIT_seed_sign)
    call CFG_get(cfg, "init_elec_energy", INIT_energy)
    call CFG_get(cfg, "init_use_neg_ion", INIT_use_nion)
    call CFG_get(cfg, "fluid_small_density", INIT_small_dens)
    call CFG_get(cfg, "init_background_density", INIT_back)
    
    
    !call TD_get_td_from_file(input_file, "AIR", &
    !     trim("elec"), x_data, y_data)
    !call LT_add_col(FL_lkp, x_data, y_data)
    !FL_if_elec = LT_get_num_cols(FL_lkp)

    !call TD_get_td_from_file(input_file, 'AIR', &
    !     trim('ion'), x_data, y_data)
    !call LT_add_col(FL_lkp, x_data, y_data)
    !FL_if_ion = LT_get_num_cols(FL_lkp)
    
  end subroutine INIT_init

  subroutine get_dens(xx, dens)
    real(dp), intent(in) :: xx
    real(dp), intent(inout) :: dens
     
     
    if (xx < INIT_DI) then
      select case (INIT_type)
      case (INIT_gaussian_type)
        dens = dens + INIT_dens * exp(-(xx - INIT_location)**2 / (2 * INIT_width**2))
      case (INIT_block_type)
        if ((xx - INIT_location) < INIT_width) dens = dens + INIT_dens
      case default
        stop "m_initial_cond is being used without proper initialization"
      end select
    end if

  end subroutine get_dens
  
  subroutine INIT_get_back(xx, elec_dens, ion_dens)
    real(dp), intent(in) :: xx
    real(dp), intent(out) :: elec_dens, ion_dens
    
    if (xx > INIT_DI) then
      elec_dens = 0.0_dp
      ion_dens = 0.0_dp
    else
      elec_dens = INIT_back
      ion_dens  = INIT_back
    end if
    
  end subroutine INIT_get_back

  subroutine INIT_get_elec_dens(xx, elec_dens)
    real(dp), intent(in)       :: xx
    real(dp), intent(inout)    :: elec_dens
    type(LT_loc_t)             :: loc
    
    if((INIT_DI-xx) > 0) then
      if (.not. INIT_use_nion .and. (INIT_seed_sign == -1 .or. INIT_seed_sign == 0)) call get_dens(xx, elec_dens)
    end if
  end subroutine INIT_get_elec_dens

  subroutine INIT_get_nion_dens(xx, nion_dens)
    real(dp), intent(in) :: xx
    real(dp), intent(inout) :: nion_dens
    
    !if (INIT_use_nion) call get_dens(xx, nion_dens)

  end subroutine INIT_get_nion_dens

  subroutine INIT_get_ion_dens(xx, ion_dens)
    real(dp), intent(in) :: xx
    real(dp), intent(out) :: ion_dens
    type(LT_loc_t)             :: loc
    
    if((INIT_DI-xx) > 0) then 
      if (INIT_seed_sign == 1 .or. INIT_seed_sign == 0) call get_dens(xx, ion_dens)
    end if
  end subroutine INIT_get_ion_dens


  subroutine INIT_get_en_dens(xx, en_dens)
    real(dp), intent(in) :: xx
    real(dp), intent(out) :: en_dens
    real(dp) :: elec_dens
    call INIT_get_elec_dens(xx, elec_dens)
    en_dens = INIT_energy * elec_dens
  end subroutine INIT_get_en_dens

end module m_init_cond_1d
