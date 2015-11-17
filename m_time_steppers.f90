!> Explicit time steppers. Variants in 1d and 2d are listed, which are identical, except for the
!! type of argument they expect (1d array of values vs. 2d array of values)
module m_time_steppers

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)

   interface
      subroutine STEP_derivs_1d_type(vars, time, derivs)
         import
         real(dp), intent(in) :: vars(:), time
         real(dp), intent(out) :: derivs(:)
      end subroutine STEP_derivs_1d_type

      subroutine STEP_derivs_2d_type(vars, time, derivs)
         import
         real(dp), intent(in) :: vars(:,:), time
         real(dp), intent(out) :: derivs(:,:)
      end subroutine STEP_derivs_2d_type

      subroutine STEP_derivs_3d_type(vars, time, derivs)
         import
         real(dp), intent(in) :: vars(:,:,:), time
         real(dp), intent(out) :: derivs(:,:,:)
      end subroutine STEP_derivs_3d_type

      subroutine STEP_1d_type(vars, time, dt, pptr_get_derivs)
         import
         real(dp), intent(inout)        :: vars(:), time
         real(dp), intent(in)           :: dt
         procedure(STEP_derivs_1d_type) :: pptr_get_derivs
      end subroutine STEP_1d_type

      subroutine STEP_2d_type(vars, time, dt, pptr_get_derivs)
         import
         real(dp), intent(inout)        :: vars(:,:), time
         real(dp), intent(in)           :: dt
         procedure(STEP_derivs_2d_type) :: pptr_get_derivs
      end subroutine STEP_2d_type

      subroutine STEP_3d_type(vars, time, dt, pptr_get_derivs)
         import
         real(dp), intent(inout)        :: vars(:,:,:), time
         real(dp), intent(in)           :: dt
         procedure(STEP_derivs_3d_type) :: pptr_get_derivs
      end subroutine STEP_3d_type

   end interface

   ! Procedure types
   public :: STEP_1d_type
   public :: STEP_2d_type

   ! Procedures
   public :: STEP_expl_trap_2d
   public :: STEP_forward_euler_1d
   public :: STEP_forward_euler_2d
   public :: STEP_forward_euler_3d
   public :: STEP_rk2_1d
   public :: STEP_rk2_2d
   public :: STEP_rk2_3d
   public :: STEP_rk4_1d
   public :: STEP_rk4_2d
   public :: STEP_rk4_3d
   public :: STEP_rk4a_1d
   public :: STEP_rk4a_2d
   public :: STEP_rk4a_3d

contains

   subroutine STEP_forward_euler_1d(vars, time, dt, pptr_get_derivs)
      real(dp), intent(inout)        :: vars(:), time
      real(dp), intent(in)           :: dt
      procedure(STEP_derivs_1d_type) :: pptr_get_derivs
      real(dp)                       :: derivs(size(vars))

      call pptr_get_derivs(vars, time, derivs)
      vars = vars + dt * derivs
      time = time + dt
   end subroutine STEP_forward_euler_1d

   subroutine STEP_forward_euler_2d(vars, time, dt, pptr_get_derivs)
      real(dp), intent(inout)        :: vars(:,:), time
      real(dp), intent(in)           :: dt
      procedure(STEP_derivs_2d_type) :: pptr_get_derivs
      real(dp)                       :: derivs(size(vars,1), size(vars,2))

      call pptr_get_derivs(vars, time, derivs)
      vars = vars + dt * derivs
      time = time + dt
   end subroutine STEP_forward_euler_2d

   subroutine STEP_forward_euler_3d(vars, time, dt, pptr_get_derivs)
      real(dp), intent(inout)        :: vars(:,:,:), time
      real(dp), intent(in)           :: dt
      procedure(STEP_derivs_3d_type) :: pptr_get_derivs
      real(dp)                       :: derivs(size(vars,1), size(vars,2), size(vars,3))

      call pptr_get_derivs(vars, time, derivs)
      vars = vars + dt * derivs
      time = time + dt
   end subroutine STEP_forward_euler_3d

   subroutine STEP_rk2_1d(vars, time, dt, pptr_get_derivs)
      real(dp), intent(inout)        :: vars(:), time
      real(dp), intent(in)           :: dt
      procedure(STEP_derivs_1d_type) :: pptr_get_derivs
      real(dp)                       :: derivs(size(vars))

      ! Step 1 (at initial time)
      call pptr_get_derivs(vars, time, derivs)

      ! Step 2 (at initial time + dt/2)
      call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)

      ! Result (at initial time + dt)
      vars = vars + dt * derivs
      time = time + dt
   end subroutine STEP_rk2_1d

   subroutine STEP_rk2_2d(vars, time, dt, pptr_get_derivs)
      real(dp), intent(inout)        :: vars(:,:), time
      real(dp), intent(in)           :: dt
      procedure(STEP_derivs_2d_type) :: pptr_get_derivs
      real(dp)                       :: derivs(size(vars,1), size(vars,2))

      ! Step 1 (at initial time)
      call pptr_get_derivs(vars, time, derivs)

      ! Step 2 (at initial time + dt/2)
      call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)

      ! Result (at initial time + dt)
      vars = vars + dt * derivs
      time = time + dt
   end subroutine STEP_rk2_2d

   subroutine STEP_expl_trap_2d(vars, time, dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:,:), time
     real(dp), intent(in)           :: dt
     procedure(STEP_derivs_2d_type) :: pptr_get_derivs
     real(dp)                       :: vars_t0(size(vars,1), size(vars,2))

     vars_t0 = vars

     ! Step 1 (at initial time)
     call STEP_forward_euler_2d(vars, time, dt, pptr_get_derivs)
     call STEP_forward_euler_2d(vars, time, dt, pptr_get_derivs)

     ! Result (at initial time + dt)
     vars = 0.5_dp * (vars + vars_t0)
     time = time - dt
   end subroutine STEP_expl_trap_2d

   subroutine STEP_rk2_3d(vars, time, dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:,:,:), time
     real(dp), intent(in)           :: dt
     procedure(STEP_derivs_3d_type) :: pptr_get_derivs
     real(dp)                       :: derivs(size(vars,1), size(vars,2), size(vars,3))

     ! Step 1 (at initial time)
     call pptr_get_derivs(vars, time, derivs)

     ! Step 2 (at initial time + dt/2)
     call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)

     ! Result (at initial time + dt)
     vars = vars + dt * derivs
     time = time + dt
   end subroutine STEP_rk2_3d

   subroutine STEP_rk4_1d(vars, time, dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:), time
     real(dp), intent(in)           :: dt
     procedure(STEP_derivs_1d_type) :: pptr_get_derivs
     real(dp)                       :: derivs(size(vars))
     real(dp)                       :: sum_derivs(size(vars))
     real(dp), parameter            :: one_sixth = 1 / 6.0_dp

     ! Step 1 (at initial time)
     call pptr_get_derivs(vars, time, derivs)
     sum_derivs = derivs

     ! Step 2 (at initial time + dt/2)
     call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)
     sum_derivs = sum_derivs + 2 * derivs

     ! Step 3 (at initial time + dt/2)
     call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)
     sum_derivs = sum_derivs + 2 * derivs

     ! Step 4 (at initial time + dt)
     call pptr_get_derivs(vars + dt * derivs, time + dt, derivs)
     sum_derivs = sum_derivs + derivs

     ! Combine time derivatives at steps
     vars = vars + one_sixth * dt * sum_derivs
     time = time + dt
   end subroutine STEP_rk4_1d

   subroutine STEP_rk4_2d(vars, time, dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:,:), time
     real(dp), intent(in)           :: dt
     procedure(STEP_derivs_2d_type) :: pptr_get_derivs
     real(dp)                       :: derivs(size(vars,1), size(vars,2))
     real(dp)                       :: sum_derivs(size(vars,1), size(vars,2))
     real(dp), parameter            :: one_sixth = 1 / 6.0_dp

     ! Step 1 (at initial time)
     call pptr_get_derivs(vars, time, derivs)
     sum_derivs = derivs

     ! Step 2 (at initial time + dt/2)
     call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)
     sum_derivs = sum_derivs + 2 * derivs

     ! Step 3 (at initial time + dt/2)
     call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)
     sum_derivs = sum_derivs + 2 * derivs

     ! Step 4 (at initial time + dt)
     call pptr_get_derivs(vars + dt * derivs, time + dt, derivs)
     sum_derivs = sum_derivs + derivs

     ! Combine time derivatives at steps
     vars = vars + one_sixth * dt * sum_derivs
     time = time + dt
   end subroutine STEP_rk4_2d

   subroutine STEP_rk4_3d(vars, time, dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:,:,:), time
     real(dp), intent(in)           :: dt
     procedure(STEP_derivs_3d_type) :: pptr_get_derivs
     real(dp)                       :: derivs(size(vars,1), size(vars,2), size(vars,3))
     real(dp)                       :: sum_derivs(size(vars,1), size(vars,2), size(vars,3))
     real(dp), parameter            :: one_sixth = 1 / 6.0_dp

     ! Step 1 (at initial time)
     call pptr_get_derivs(vars, time, derivs)
     sum_derivs = derivs

     ! Step 2 (at initial time + dt/2)
     call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)
     sum_derivs = sum_derivs + 2 * derivs

     ! Step 3 (at initial time + dt/2)
     call pptr_get_derivs(vars + 0.5_dp * dt * derivs, time + 0.5_dp * dt, derivs)
     sum_derivs = sum_derivs + 2 * derivs

     ! Step 4 (at initial time + dt)
     call pptr_get_derivs(vars + dt * derivs, time + dt, derivs)
     sum_derivs = sum_derivs + derivs

     ! Combine time derivatives at steps
     vars = vars + one_sixth * dt * sum_derivs
     time = time + dt
   end subroutine STEP_rk4_3d

   subroutine STEP_rk4a_1d(vars, max_errs, time, dt, max_dt, new_dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:), time
     real(dp), intent(in)           :: max_errs(:), dt, max_dt
     real(dp), intent(out)          :: new_dt
     procedure(STEP_derivs_1d_type) :: pptr_get_derivs
     real(dp)                       :: new_vars(size(vars))
     real(dp)                       :: cpy_vars(size(vars))
     real(dp)                       :: new_time, cpy_time, max_err_ratio
     real(dp), parameter            :: one_fifteenth = 1 / 15.0_dp
     real(dp), parameter            :: safety_fac = 0.9_dp

     new_dt = min(max_dt, dt)

     do
        new_vars = vars
        cpy_vars = vars
        new_time = time
        cpy_time = time

        ! Two halve steps
        call STEP_rk4_1d(new_vars, new_time, 0.5_dp * new_dt, pptr_get_derivs)
        call STEP_rk4_1d(new_vars, new_time, 0.5_dp * new_dt, pptr_get_derivs)

        ! One big step
        call STEP_rk4_1d(cpy_vars, cpy_time, new_dt, pptr_get_derivs)

        ! cpy_vars now holds the estimated truncation error
        cpy_vars = new_vars - cpy_vars

        ! Check maximum absolute error
        max_err_ratio = maxval(abs(cpy_vars/max_errs))

        ! Guess a new dt
        if (max_err_ratio > 1.0_dp) then
           ! Try again with smaller timestep
           new_dt = safety_fac * new_dt / max_err_ratio**(1/3.0_dp)
           cycle
        else
           ! Routine is done, set results. At most double dt.
           new_dt = safety_fac * new_dt / max(0.5_dp, max_err_ratio**0.25_dp)
           new_dt = min(max_dt, new_dt)
           vars = new_vars ! You can add local extrapolation with + cpy_vars * one_fifteenth
           time = new_time
           exit
        end if
     end do

   end subroutine STEP_rk4a_1d

   subroutine STEP_rk4a_2d(vars, max_errs, time, dt, max_dt, new_dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:,:), time
     real(dp), intent(in)           :: max_errs(:,:), dt, max_dt
     real(dp), intent(out)          :: new_dt
     procedure(STEP_derivs_2d_type) :: pptr_get_derivs
     real(dp)                       :: new_vars(size(vars,1), size(vars,2))
     real(dp)                       :: cpy_vars(size(vars,1), size(vars,2))
     real(dp)                       :: new_time, cpy_time, max_err_ratio
     real(dp), parameter            :: one_fifteenth = 1 / 15.0_dp
     real(dp), parameter            :: safety_fac = 0.9_dp

     new_dt = min(max_dt, dt)

     do
        new_vars = vars
        cpy_vars = vars
        new_time = time
        cpy_time = time

        ! Two halve steps
        call STEP_rk4_2d(new_vars, new_time, 0.5_dp * new_dt, pptr_get_derivs)
        call STEP_rk4_2d(new_vars, new_time, 0.5_dp * new_dt, pptr_get_derivs)

        ! One big step
        call STEP_rk4_2d(cpy_vars, cpy_time, new_dt, pptr_get_derivs)

        ! cpy_vars now holds the estimated truncation error
        cpy_vars = new_vars - cpy_vars

        ! Check maximum absolute error
        max_err_ratio = maxval(abs(cpy_vars/max_errs))

        ! Guess a new dt
        if (max_err_ratio > 1.0_dp) then
           ! Try again with smaller timestep
           new_dt = safety_fac * new_dt / max_err_ratio**(1/3.0_dp)
           cycle
        else
           ! Routine is done, set results. At most double dt.
           new_dt = safety_fac * new_dt / max(0.5_dp, max_err_ratio**0.25_dp)
           new_dt = min(max_dt, new_dt)
           vars = new_vars ! You can add local extrapolation with + cpy_vars * one_fifteenth
           time = new_time
           exit
        end if
     end do

   end subroutine STEP_rk4a_2d

   subroutine STEP_rk4a_3d(vars, max_errs, time, dt, max_dt, new_dt, pptr_get_derivs)
     real(dp), intent(inout)        :: vars(:,:,:), time
     real(dp), intent(in)           :: max_errs(:,:,:), dt, max_dt
     real(dp), intent(out)          :: new_dt
     procedure(STEP_derivs_3d_type) :: pptr_get_derivs
     real(dp)                       :: new_vars(size(vars,1), size(vars,2), size(vars,3))
     real(dp)                       :: cpy_vars(size(vars,1), size(vars,2), size(vars,3))
     real(dp)                       :: new_time, cpy_time, max_err_ratio
     real(dp), parameter            :: one_fifteenth = 1 / 15.0_dp
     real(dp), parameter            :: safety_fac = 0.9_dp

     new_dt = min(max_dt, dt)

     do
        new_vars = vars
        cpy_vars = vars
        new_time = time
        cpy_time = time

        ! Two halve steps
        call STEP_rk4_3d(new_vars, new_time, 0.5_dp * new_dt, pptr_get_derivs)
        call STEP_rk4_3d(new_vars, new_time, 0.5_dp * new_dt, pptr_get_derivs)

        ! One big step
        call STEP_rk4_3d(cpy_vars, cpy_time, new_dt, pptr_get_derivs)

        ! cpy_vars now holds the estimated truncation error
        cpy_vars = new_vars - cpy_vars

        ! Check maximum absolute error
        max_err_ratio = maxval(abs(cpy_vars/max_errs))

        ! Guess a new dt
        if (max_err_ratio > 1.0_dp) then
           ! Try again with smaller timestep
           new_dt = safety_fac * new_dt / max_err_ratio**(1/3.0_dp)
           cycle
        else
           ! Routine is done, set results. At most double dt.
           new_dt = safety_fac * new_dt / max(0.5_dp, max_err_ratio**0.25_dp)
           new_dt = min(max_dt, new_dt)
           vars = new_vars ! You can add local extrapolation with + cpy_vars * one_fifteenth
           time = new_time
           exit
        end if
     end do

   end subroutine STEP_rk4a_3d

 end module m_time_steppers
