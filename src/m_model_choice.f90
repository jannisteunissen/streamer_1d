module m_model_choice

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  ! Enumerate the different simulation types
  integer, parameter :: MODEL_part      = 1
  integer, parameter :: MODEL_fluid_lfa = 2
  integer, parameter :: MODEL_fluid_ee  = 3

  integer, protected :: MODEL_type = -1

  ! Public types
  public :: MODEL_type
  public :: MODEL_part, MODEL_fluid_lfa, MODEL_fluid_ee

  ! Public routines
  public :: MODEL_initialize

contains

  subroutine MODEL_initialize(MODEL_type_name)
    character(len=*), intent(in) :: MODEL_type_name

    select case (MODEL_type_name)
    case ("part")
       MODEL_type = MODEL_part
    case ("fluid_lfa")
       MODEL_type = MODEL_fluid_lfa
    case ("fluid_ee")
       MODEL_type = MODEL_fluid_ee
    case default
       print *, "MODEL_initialize: Invalid simulation type given: " // MODEL_type_name
       stop
    end select

  end subroutine MODEL_initialize

end module m_model_choice
