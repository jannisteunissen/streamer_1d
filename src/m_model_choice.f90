module m_model_choice

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  ! Enumerate the different simulation types
  integer, parameter :: MODEL_part      = 1
  integer, parameter :: MODEL_fluid_lfa = 2
  integer, parameter :: MODEL_fluid_ee  = 3

  integer :: MODEL_type = -1

  ! Public types
  public :: MODEL_part, MODEL_fluid_lfa, MODEL_fluid_ee

  ! Public routines
  public :: MODEL_initialize
  public :: MODEL_get_type

contains

  subroutine MODEL_initialize(sim_type_name)
    character(len=*), intent(in) :: sim_type_name

    select case (sim_type_name)
    case ("part")
       MODEL_type = MODEL_part
    case ("fluid_lfa")
       MODEL_type = MODEL_fluid_lfa
    case ("fluid_ee")
       MODEL_type = MODEL_fluid_ee
    case default
       print *, "MODEL_initialize: Invalid simulation type given: " // sim_type_name
       stop
    end select

  end subroutine MODEL_initialize

  integer function MODEL_get_type()
    MODEL_get_type = MODEL_type
    if (MODEL_get_type == -1) then
       print *, "MODEL_get_type: call MODEL_initialize first"
       stop
    end if
  end function MODEL_get_type

end module m_model_choice
