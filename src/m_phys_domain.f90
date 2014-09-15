module m_phys_domain

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  real(dp), protected         :: PD_length
  real(dp), protected         :: PD_dx, PD_inv_dx
  integer                     :: PD_grid_size

contains

  subroutine PD_set(dx, grid_size)
    real(dp), intent(in) :: dx
    integer, intent(in) :: grid_size

    PD_dx = dx
    PD_inv_dx = 1/dx
    PD_grid_size = grid_size
    PD_length = dx * (grid_size-1)
  end subroutine PD_set

end module m_phys_domain