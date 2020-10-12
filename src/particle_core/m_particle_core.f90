! Authors: Jannis Teunissen, concepts based on work of Chao Li, Margreet Nool
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Core particle module. The user has to supply his own routines for customization.
module m_particle_core
  use m_lookup_table
  use m_linked_list
  use m_random
  use m_cross_sec

  implicit none
  private

  integer, parameter  :: dp               = kind(0.0d0)

  !> Special weight value indicating a particle has been removed
  real(dp), parameter :: PC_dead_weight   = -1e100_dp

  !> The maximum number of collisions
  !> \todo Consider making this a variable again (but check OpenMP performance)
  integer, parameter, public :: PC_max_num_coll = 100

  !> Maximum number of particles created by one collision (e.g., ionization)
  integer, parameter :: PC_coll_max_part_out = 2

  !> Buffer size for advancing particle in parallel (not important for users)
  integer, parameter :: PC_advance_buf_size = 1000

  !> Initial size of an event list (will be automatically increased when required)
  integer, parameter :: PC_event_list_init_size = 10*1000

  !> Integer indicating the 'collision type' of events in which particles went
  !> outside the domain
  integer, parameter, public :: PC_particle_went_out = -1

  !> The particle type
  type, public :: PC_part_t
     integer  :: ptype  = 0 !< Type of particle (not used yet)
     integer  :: id     = 0 !< Can be used to e.g. optimize code
     real(dp) :: x(3)   = 0 !< Position
     real(dp) :: v(3)   = 0 !< Velocity
     real(dp) :: a(3)   = 0 !< Acceleration
     real(dp) :: w      = 0 !< Weight
     real(dp) :: t_left = 0 !< Propagation time left
  end type PC_part_t

  !> An event (particle collision)
  type, public :: PC_event_t
     type(PC_part_t) :: part  !< Particle that had collision
     integer         :: cix   !< Collision index
     integer         :: ctype !< Collision type
  end type PC_event_t

  !> Particle buffer for parallel simulations
  type, public        :: PC_buf_t
     !> Index for newly created particles
     integer          :: i_part  = 0
     !> Buffer for newly created particles
     type(PC_part_t)  :: part(PC_advance_buf_size)
     !> Index for removed particles
     integer          :: i_rm    = 0
     !> Buffer for removed particles
     integer          :: rm(PC_advance_buf_size)
     !> Index for events
     integer          :: i_event = 0
     !> Buffer for events
     type(PC_event_t) :: event(PC_advance_buf_size)
  end type PC_buf_t

  !> Type for directly specifying collision rates
  type, public :: rate_func_t
     type(CS_coll_t)                       :: coll
     procedure(rate_func), pointer, nopass :: ptr => null()
  end type rate_func_t

  !> Particle core type, storing the particles and the collisions
  type, public                    :: PC_t
     !> Array storing the particles
     type(PC_part_t), allocatable :: particles(:)

     !> Number of particles
     integer                      :: n_part

     !> List of collisions
     type(CS_coll_t), allocatable :: colls(:)

     !> Number of collisions
     integer                      :: n_colls

     !> Whether a collision should be stored as an event
     logical, allocatable         :: coll_is_event(:)

     !> Lookup table with collision rates
     type(LT_t)                   :: rate_lt
     real(dp)                     :: max_rate, inv_max_rate

     !> List of particles to be removed
     type(LL_int_head_t)          :: clean_list

     !> Fixed mass for the particles
     real(dp)                     :: mass

     !> State of random number generator
     type(RNG_t)                  :: rng

     !> Maximum time step for particle mover
     real(dp)                     :: dt_max = huge(1.0_dp)

     !> Magnetic field (for Boris mover)
     real(dp)                     :: B_vec(3) = [0.0_dp, 0.0_dp, 0.0_dp]

     !> Number of stored events
     integer                       :: n_events = 0

     !> List of events
     type(PC_event_t), allocatable :: event_list(:)

     !> If assigned call this method after moving particles to check whether
     !> they are outside the computational domain (then it should return a
     !> positive value)
     procedure(p_to_int_f), pointer, nopass :: outside_check  => null()

     !> The particle mover
     procedure(subr_mover), pointer :: particle_mover => null()

     !> Called after the particle mover
     procedure(subr_after), pointer :: after_mover => null()

     !> The method to get particle accelerations
     procedure(p_to_r3_f), pointer, nopass :: accel_function => null()

     !> To set the velocity of tracer particles
     procedure(p_to_r3_f), pointer, nopass :: tracer_velocity => null()

   contains

     ! A list of methods
     procedure, non_overridable :: initialize
     procedure, non_overridable :: resize_part_list
     procedure, non_overridable :: remove_particles
     procedure, non_overridable :: advance
     procedure, non_overridable :: advance_openmp
     procedure, non_overridable :: clean_up
     procedure, non_overridable :: create_part
     procedure, non_overridable :: add_part
     procedure, non_overridable :: remove_part
     procedure, non_overridable :: periodify
     procedure, non_overridable :: translate
     procedure, non_overridable :: get_num_sim_part
     procedure, non_overridable :: get_num_real_part
     procedure, non_overridable :: set_accel
     procedure, non_overridable :: get_max_coll_rate
     procedure, non_overridable :: get_colls_of_type
     procedure, non_overridable :: loop_iopart
     procedure, non_overridable :: loop_ipart
     procedure, non_overridable :: compute_scalar_sum
     procedure, non_overridable :: compute_vector_sum
     procedure, non_overridable :: move_and_collide
     procedure, non_overridable :: use_cross_secs
     procedure, non_overridable :: use_rate_funcs
     procedure, non_overridable :: get_mean_energy
     procedure, non_overridable :: get_coll_rates
     procedure, non_overridable :: check_space

     procedure, non_overridable :: sort
     procedure, non_overridable :: merge_and_split
     procedure, non_overridable :: merge_and_split_range
     procedure, non_overridable :: histogram

     procedure, non_overridable :: get_num_colls
     procedure, non_overridable :: get_colls
     procedure, non_overridable :: get_coeffs

     procedure, non_overridable :: init_from_file
     procedure, non_overridable :: to_file
  end type PC_t

  type, abstract, public :: PC_bin_t
     integer :: n_bins
   contains
     procedure(bin_f), deferred :: bin_func
  end type PC_bin_t

  abstract interface

     subroutine p_to_r3_p(my_part, my_vec)
       import
       type(PC_part_t), intent(in) :: my_part
       real(dp), intent(out)       :: my_vec(3)
     end subroutine p_to_r3_p

     function p_to_r3_f(my_part) result(my_vec)
       import
       type(PC_part_t), intent(inout) :: my_part
       real(dp)                    :: my_vec(3)
     end function p_to_r3_f

     real(dp) function p_to_r_f(my_part)
       import
       type(PC_part_t), intent(in) :: my_part
     end function p_to_r_f

     logical function p_to_logic_f(my_part)
       import
       type(PC_part_t), intent(in) :: my_part
     end function p_to_logic_f

     integer function p_to_int_f(my_part)
       import
       type(PC_part_t), intent(inout) :: my_part
     end function p_to_int_f

     integer function bin_f(binner, my_part)
       import
       class(PC_bin_t), intent(in) :: binner
       type(PC_part_t), intent(in) :: my_part
     end function bin_f

     subroutine subr_mover(self, part, dt)
       import
       class(PC_t), intent(in)        :: self
       type(PC_part_t), intent(inout) :: part
       real(dp), intent(in)           :: dt
     end subroutine subr_mover

     subroutine subr_after(self, dt)
       import
       class(PC_t), intent(inout)     :: self
       real(dp), intent(in)           :: dt
     end subroutine subr_after

     subroutine sub_merge(part_a, part_b, part_out, rng)
       import
       type(PC_part_t), intent(in)  :: part_a, part_b
       type(PC_part_t), intent(out) :: part_out
       type(RNG_t), intent(inout)   :: rng
     end subroutine sub_merge

     subroutine sub_split(part_a, w_ratio, part_out, n_part_out, rng)
       import
       type(PC_part_t), intent(in)    :: part_a
       real(dp), intent(in)           :: w_ratio
       type(PC_part_t), intent(inout) :: part_out(:)
       integer, intent(inout)         :: n_part_out
       type(RNG_t), intent(inout)     :: rng
     end subroutine sub_split

     function rate_func(v) result(rate)
       import
       real(dp), intent(in) :: v
       real(dp)             :: rate
     end function rate_func
  end interface

  ! Public types
  public :: PC_particle_went_out

  ! Public procedures
  public :: PC_merge_part_rxv
  public :: PC_split_part
  public :: PC_v_to_en
  public :: PC_speed_to_en

  public :: PC_verlet_advance
  public :: PC_verlet_correct_accel
  public :: PC_boris_advance
  public :: PC_tracer_advance_midpoint
  public :: PC_after_dummy

contains

  !> Initialization routine for the particle module
  subroutine initialize(self, mass, n_part_max, rng_seed)
    use m_units_constants
    class(PC_t), intent(inout)      :: self
    real(dp), intent(in)            :: mass
    integer, intent(in)             :: n_part_max
    integer, intent(in), optional   :: rng_seed(4)
    integer, parameter              :: i8 = selected_int_kind(18)
    integer(i8)                     :: rng_seed_8byte(2)

    allocate(self%particles(n_part_max))
    self%mass   = mass
    self%n_part = 0

    if (present(rng_seed)) then
       rng_seed_8byte = transfer(rng_seed, rng_seed_8byte)
       call self%rng%set_seed(rng_seed_8byte)
    else
       call self%rng%set_random_seed()
    end if

    ! Set default particle mover
    if (.not. associated(self%particle_mover)) &
         self%particle_mover => PC_verlet_advance

    if (.not. associated(self%after_mover)) &
         self%after_mover => PC_after_dummy

  end subroutine initialize

  subroutine check_methods(self)
    class(PC_t), intent(inout) :: self

    if (associated(self%particle_mover, PC_tracer_advance_midpoint)) then
       if (.not. associated(self%tracer_velocity)) &
            error stop "Set the tracer_velocity method for tracers"
    end if

  end subroutine check_methods

  !> Initialization routine for the particle module
  subroutine init_from_file(self, param_file, lt_file, rng_seed)
    use m_cross_sec
    class(PC_t), intent(inout)    :: self
    character(len=*), intent(in)  :: param_file, lt_file
    integer, intent(in), optional :: rng_seed(4)

    integer, parameter            :: i8 = selected_int_kind(18)
    integer(i8)                   :: rng_seed_8byte(2)
    integer                       :: my_unit
    integer                       :: n_part_max

    open(newunit=my_unit, file=trim(param_file), form='UNFORMATTED', &
         access='STREAM', status='OLD')

    read(my_unit) n_part_max
    read(my_unit) self%n_colls
    allocate(self%coll_is_event(self%n_colls))
    read(my_unit) self%coll_is_event
    read(my_unit) self%mass
    read(my_unit) self%max_rate
    allocate(self%colls(self%n_colls))
    read(my_unit) self%colls

    close(my_unit)

    self%inv_max_rate = 1 / self%max_rate
    self%n_part       = 0
    allocate(self%particles(n_part_max))

    call LT_from_file(self%rate_lt, lt_file)

    if (present(rng_seed)) then
       rng_seed_8byte = transfer(rng_seed, rng_seed_8byte)
       call self%rng%set_seed(rng_seed_8byte)
    else
       call self%rng%set_random_seed()
    end if
  end subroutine init_from_file

  subroutine to_file(self, param_file, lt_file)
    use m_cross_sec
    class(PC_t), intent(inout)   :: self
    character(len=*), intent(in) :: param_file, lt_file
    integer                      :: my_unit

    open(newunit=my_unit, file=trim(param_file), form='UNFORMATTED', &
         access='STREAM', status='REPLACE')
    write(my_unit) size(self%particles)
    write(my_unit) self%n_colls
    write(my_unit) self%coll_is_event
    write(my_unit) self%mass
    write(my_unit) self%max_rate
    write(my_unit) self%colls
    close(my_unit)

    call LT_to_file(self%rate_lt, lt_file)
  end subroutine to_file

  subroutine get_colls_of_type(pc, ctype, ixs)
    class(PC_t), intent(in) :: pc
    integer, intent(in) :: ctype
    integer, intent(inout), allocatable :: ixs(:)
    integer :: nn, i

    nn = count(pc%colls(:)%type == ctype)
    allocate(ixs(nn))

    nn = 0
    do i = 1, pc%n_colls
       if (pc%colls(i)%type == ctype) then
          nn = nn + 1
          ixs(nn) = i
       end if
    end do
  end subroutine get_colls_of_type

  subroutine remove_particles(self)
    class(PC_t), intent(inout) :: self
    self%n_part = 0
    call LL_clear(self%clean_list)
  end subroutine remove_particles

  subroutine resize_part_list(self, new_size)
    class(PC_t), intent(inout) :: self
    integer, intent(in)          :: new_size
    type(PC_part_t), allocatable :: parts_copy(:)

    allocate(parts_copy(self%n_part))
    deallocate(self%particles)
    allocate(self%particles(new_size))
    self%particles(1:self%n_part) = parts_copy
  end subroutine resize_part_list

  subroutine check_events_allocated(self)
    class(PC_t), intent(inout) :: self

    if (.not. allocated(self%event_list)) then
       allocate(self%event_list(PC_event_list_init_size))
    end if
  end subroutine check_events_allocated

  subroutine ensure_events_storage(self, n_new)
    class(PC_t), intent(inout) :: self
    integer, intent(in)        :: n_new
    integer                    :: n_stored
    type(PC_event_t), allocatable :: cpy(:)

    n_stored = self%n_events
    if (size(self%event_list) < n_stored + n_new) then
       allocate(cpy(self%n_events))
       cpy = self%event_list(1:self%n_events)
       deallocate(self%event_list)
       allocate(self%event_list(2 * (n_stored + n_new)))
       self%event_list(1:n_stored) = cpy(1:n_stored)
    end if
  end subroutine ensure_events_storage

  !> Limit advance time so that not too many collision happen per particle per
  !> time step (otherwise the buffer for new particles could overflow)
  subroutine limit_advance_dt(self, dt, n_steps, dt_step)
    class(PC_t), intent(in) :: self
    real(dp), intent(in)    :: dt
    integer, intent(out)    :: n_steps
    real(dp), intent(out)   :: dt_step
    real(dp)                :: max_dt

    if (dt < 0) error stop "dt < 0"
    max_dt = 0.25_dp * self%inv_max_rate * PC_advance_buf_size
    n_steps = max(1, ceiling(dt / max_dt))
    dt_step = dt / n_steps
  end subroutine limit_advance_dt

  !> Advance particles over dt, using one or more step
  subroutine advance(self, dt)
    class(PC_t), intent(inout)       :: self
    real(dp), intent(in)             :: dt
    integer                          :: n, n_steps
    real(dp)                         :: dt_step

    call check_methods(self)
    call limit_advance_dt(self, dt, n_steps, dt_step)
    call check_events_allocated(self)

    ! Advance the particles in one or more steps, which makes sure the buffers
    ! for newly created particles and events are large enough
    do n = 1, n_steps
       call advance_step(self, dt_step)
    end do

    call self%clean_up()
  end subroutine advance

  !> Advance particles over dt in a single step
  subroutine advance_step(self, dt)
    class(PC_t), intent(inout)       :: self
    real(dp), intent(in)             :: dt
    integer                          :: n
    type(PC_buf_t)                   :: buffer

    self%particles(1:self%n_part)%t_left = dt
    n = 1

    do while (n <= self%n_part)
       call self%move_and_collide(n, self%rng, buffer)
       call handle_buffer(self, buffer, 0)
       n = n + 1
    end do
  end subroutine advance_step

  subroutine init_buffer(buffer)
    type(PC_buf_t), intent(inout) :: buffer
    buffer%i_part  = 0
    buffer%i_rm    = 0
    buffer%i_event = 0
  end subroutine init_buffer

  !> If the buffers for a thread are getting too full, empty them
  subroutine handle_buffer(self, buffer, max_size)
    class(PC_t), intent(inout)       :: self
    type(PC_buf_t), intent(inout)    :: buffer
    !> Keep at most this many buffered items
    integer, intent(in)              :: max_size
    integer                          :: i

    ! The buffer for new particles
    if (buffer%i_part > max_size) then
       !$omp critical
       i = self%n_part
       self%n_part = self%n_part + buffer%i_part
       !$omp end critical
       call self%check_space(i + buffer%i_part)
       self%particles(i+1:i+buffer%i_part) = buffer%part(1:buffer%i_part)
       buffer%i_part = 0
    end if

    ! The buffer for removed particles
    if (buffer%i_rm > max_size) then
       !$omp critical
       do i = 1, buffer%i_rm
          call self%remove_part(buffer%rm(i))
       end do
       !$omp end critical
       buffer%i_rm = 0
    end if

    ! The buffer for events
    if (buffer%i_event > max_size) then
       !$omp critical
       i = self%n_events
       call ensure_events_storage(self, buffer%i_event)
       self%n_events = self%n_events + buffer%i_event
       !$omp end critical
       self%event_list(i+1:i+buffer%i_event) = buffer%event(1:buffer%i_event)
       buffer%i_event = 0
    end if
  end subroutine handle_buffer

  !> Advance the particles over dt in parallel, using one or more steps
  subroutine advance_openmp(self, dt)
    use omp_lib
    class(PC_t), intent(inout)       :: self
    real(dp), intent(in)             :: dt
    integer                          :: n, n_steps
    real(dp)                         :: dt_step
    type(prng_t)                     :: prng

    call check_methods(self)
    call check_events_allocated(self)
    call limit_advance_dt(self, dt, n_steps, dt_step)
    call prng%init_parallel(omp_get_max_threads(), self%rng)

    ! Advance the particles in one or more steps, which makes sure the buffers
    ! for newly created particles and events are large enough
    do n = 1, n_steps
       call advance_openmp_step(self, dt_step, prng)
       call self%clean_up()
    end do

    ! Update self%rng (otherwise it would always stay the same)
    call prng%update_seed(self%rng)

  end subroutine advance_openmp

  !> Advance the particles over dt in parallel, in a single step
  subroutine advance_openmp_step(self, dt, prng)
    use omp_lib
    class(PC_t), intent(inout)       :: self
    real(dp), intent(in)             :: dt
    type(prng_t), intent(inout)      :: prng
    type(PC_buf_t)                   :: buffer
    integer                          :: n, tid, n_lo, n_hi

    self%particles(1:self%n_part)%t_left = dt

    !$omp parallel private(n, tid, n_lo, n_hi, buffer)
    call init_buffer(buffer)    ! Make sure private copies are initialized
    tid     = omp_get_thread_num() + 1
    n_lo    = 1

    do
       n_hi = self%n_part
       !$omp barrier
       !$omp do
       do n = n_lo, n_hi
          call self%move_and_collide(n, prng%rngs(tid), buffer)
          ! Make sure buffers are not getting too full
          call handle_buffer(self, buffer, PC_advance_buf_size/2)
       end do
       !$omp end do

       ! Ensure buffers are empty at the end of the loop
       call handle_buffer(self, buffer, 0)

       ! Ensure all particles have been added before the test below
       !$omp barrier

       ! Check if more loops are required
       if (self%n_part > n_hi) then
          n_lo = n_hi + 1
       else
          exit
       end if
    end do
    !$omp end parallel

  end subroutine advance_openmp_step

  !> Advance particles and collide them with neutrals.
  subroutine move_and_collide(self, ix, rng, buffer)
    use m_cross_sec
    class(PC_t), intent(inout)    :: self
    integer, intent(in)           :: ix !< Index of particle
    type(RNG_t), intent(inout)    :: rng
    type(PC_buf_t), intent(inout) :: buffer

    integer            :: i, cIx, cType, n_coll_out
    real(dp)           :: coll_time, new_vel
    type(PC_part_t)    :: coll_out(PC_coll_max_part_out)

    associate(part => self%particles(ix))
      do
         ! If the particles are treated as a tracer, advance without collision
         if (associated(self%particle_mover, PC_tracer_advance_midpoint)) exit

         ! Get the next collision time
         coll_time = sample_coll_time(rng%unif_01(), self%inv_max_rate)

         ! If larger than t_left, advance the particle without a collision
         if (coll_time > part%t_left) exit

         ! Ensure we don't move the particle over more than dt_max
         do while (coll_time > self%dt_max)
            call self%particle_mover(part, self%dt_max)
            coll_time = coll_time - self%dt_max
         end do

         ! Move particle to collision time
         call self%particle_mover(part, coll_time)
         call handle_particles_outside(self, part, ix, buffer)
         if (part%w <= PC_dead_weight) return

         new_vel = norm2(part%v)
         cIx     = get_coll_index(self%rate_lt, self%n_colls, self%max_rate, &
              new_vel, rng%unif_01())

         if (cIx > 0) then
            ! Perform the corresponding collision
            cType    = self%colls(cIx)%type

            if (self%coll_is_event(cIx)) then
               call add_event(buffer, part, cIx, cType)
            end if

            select case (cType)
            case (CS_attach_t)
               n_coll_out = 0
            case (CS_elastic_t)
               call elastic_collision(part, coll_out, &
                    n_coll_out, self%colls(cIx), rng)
            case (CS_excite_t)
               call excite_collision(part, coll_out, &
                    n_coll_out, self%colls(cIx), rng)
            case (CS_ionize_t)
               call ionization_collision(part, coll_out, &
                    n_coll_out, self%colls(cIx), rng)
            case default
               error stop "Wrong collision type"
            end select

            ! Store the particles returned
            if (n_coll_out == 0) then
               part%w = PC_dead_weight
               buffer%i_rm = buffer%i_rm + 1
               buffer%rm(buffer%i_rm) = ix
               return
            else if (n_coll_out == 1) then
               part = coll_out(1)
            else
               part = coll_out(1)
               i = buffer%i_part
               buffer%part(i+1:i+n_coll_out-1) = &
                    coll_out(2:n_coll_out)
               buffer%i_part = i+n_coll_out-1
            end if
         end if
      end do

      ! Move particle to end of the time step
      call self%particle_mover(part, part%t_left)
      call handle_particles_outside(self, part, ix, buffer)
    end associate

  end subroutine move_and_collide

  !> Check if particles went out of the domain, and remove them if so
  subroutine handle_particles_outside(self, part, ix, buffer)
    class(PC_t), intent(inout)     :: self
    type(PC_part_t), intent(inout) :: part
    integer, intent(in)            :: ix !< Particle index
    type(PC_buf_t), intent(inout)  :: buffer
    integer                        :: cIx

    if (associated(self%outside_check)) then
       cIx = self%outside_check(part)
       if (cIx > 0) then
          call add_event(buffer, part, cIx, PC_particle_went_out)
          buffer%i_rm = buffer%i_rm + 1
          buffer%rm(buffer%i_rm) = ix
          part%w = PC_dead_weight
       end if
    end if
  end subroutine handle_particles_outside

  !> Create event for a collision
  subroutine add_event(buffer, part, cIx, cType)
    type(PC_buf_t), intent(inout) :: buffer !< Buffer for events
    type(PC_part_t), intent(in)   :: part   !< Particle creating the event
    integer, intent(in)           :: cIx    !< Index of collision
    integer, intent(in)           :: cType  !< Collision type

    buffer%i_event = buffer%i_event + 1
    buffer%event(buffer%i_event)%part  = part
    buffer%event(buffer%i_event)%cix   = cIx
    buffer%event(buffer%i_event)%ctype = cType
  end subroutine add_event

  !> Returns a sample from the exponential distribution of the collision times
  ! RNG_uniform() is uniform on [0,1), but log(0) = nan, so we take 1 - RNG_uniform()
  real(dp) function sample_coll_time(unif_01, inv_max_rate)
    real(dp), intent(in) :: unif_01, inv_max_rate
    sample_coll_time = -log(1 - unif_01) * inv_max_rate
  end function sample_coll_time

  integer function get_coll_index(rate_lt, n_colls, max_rate, &
       velocity, rand_unif)
    type(LT_t), intent(in) :: rate_lt
    real(dp), intent(IN)   :: velocity, rand_unif, max_rate
    integer, intent(in)    :: n_colls
    real(dp)               :: rand_rate
    real(dp)               :: buffer(PC_max_num_coll)

    ! Fill an array with interpolated rates
    buffer(1:n_colls) = LT_get_mcol(rate_lt, velocity)
    rand_rate         = rand_unif * max_rate
    get_coll_index    = find_index_adaptive(buffer(1:n_colls), rand_rate)

    ! If there was no collision, the index exceeds the list and is set to 0
    if (get_coll_index == n_colls+1) get_coll_index = 0
  end function get_coll_index

  real(dp) function get_max_coll_rate(self)
    class(PC_t), intent(in) :: self
    get_max_coll_rate = self%max_rate
  end function get_max_coll_rate

  !> Perform an elastic collision for particle 'll'
  subroutine elastic_collision(part_in, part_out, n_part_out, coll, rng)
    type(PC_part_t), intent(in)    :: part_in
    type(PC_part_t), intent(inout) :: part_out(:)
    integer, intent(out)           :: n_part_out
    type(CS_coll_t), intent(in)    :: coll
    type(RNG_t), intent(inout)     :: rng
    real(dp)                       :: bg_vel(3), com_vel(3)

    ! TODO: implement random bg velocity
    bg_vel      = 0.0_dp
    n_part_out  = 1
    part_out(1) = part_in

    ! Compute center of mass velocity
    com_vel = (coll%rel_mass * part_out(1)%v + bg_vel) / (1 + coll%rel_mass)

    ! Scatter in center of mass coordinates
    part_out(1)%v = part_out(1)%v - com_vel
    call scatter_isotropic(part_out(1), norm2(part_out(1)%v), rng)
    part_out(1)%v = part_out(1)%v + com_vel
  end subroutine elastic_collision

  !> Perform an excitation-collision for particle 'll'
  subroutine excite_collision(part_in, part_out, n_part_out, coll, rng)
    use m_units_constants
    type(PC_part_t), intent(in)    :: part_in
    type(PC_part_t), intent(inout) :: part_out(:)
    integer, intent(out)           :: n_part_out
    type(CS_coll_t), intent(in)    :: coll
    type(RNG_t), intent(inout)     :: rng

    real(dp)             :: energy, old_en, new_vel

    old_en  = PC_v_to_en(part_in%v, coll%part_mass)
    energy  = max(0.0_dp, old_en - coll%en_loss)
    new_vel = PC_en_to_vel(energy, coll%part_mass)

    n_part_out  = 1
    part_out(1) = part_in
    call scatter_isotropic(part_out(1), new_vel, rng)
  end subroutine excite_collision

  !> Perform an ionizing collision for particle 'll'
  subroutine ionization_collision(part_in, part_out, n_part_out, coll, rng)
    use m_units_constants
    type(PC_part_t), intent(in)    :: part_in
    type(PC_part_t), intent(inout) :: part_out(:)
    integer, intent(out)           :: n_part_out
    type(CS_coll_t), intent(in)    :: coll
    type(RNG_t), intent(inout)     :: rng
    real(dp)                       :: energy, old_en, velocity

    old_en      = PC_v_to_en(part_in%v, coll%part_mass)
    energy      = max(0.0_dp, old_en - coll%en_loss)
    velocity    = PC_en_to_vel(0.5_dp * energy, coll%part_mass)

    n_part_out  = 2
    part_out(1) = part_in
    part_out(2) = part_in
    call scatter_isotropic(part_out(1), velocity, rng)
    call scatter_isotropic(part_out(2), velocity, rng)
  end subroutine ionization_collision

  subroutine scatter_isotropic(part, vel_norm, rng)
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: vel_norm
    type(RNG_t), intent(inout)     :: rng
    real(dp)                       :: sum_sq, tmp_sqrt, rands(2)

    ! Marsaglia method for uniform sampling on sphere
    do
       rands(1) = 2 * rng%unif_01() - 1
       rands(2) = 2 * rng%unif_01() - 1
       sum_sq   = rands(1)**2 + rands(2)**2
       if (sum_sq <= 1) exit
    end do

    tmp_sqrt    = sqrt(1 - sum_sq)
    part%v(1:2) = 2 * rands(1:2) * tmp_sqrt
    part%v(3)   = 1 - 2 * sum_sq
    part%v      = part%v * vel_norm ! Normalization
  end subroutine scatter_isotropic

  !> Use a Verlet scheme to advance the particle position and velocity over time
  !> dt, and update t_left.
  subroutine PC_verlet_advance(self, part, dt)
    class(PC_t), intent(in)         :: self
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: dt

    part%x      = part%x + part%v * dt + &
         0.5_dp * part%a * dt**2
    part%v      = part%v + part%a * dt
    part%t_left = part%t_left - dt
  end subroutine PC_verlet_advance

  !> Advance the particle position and velocity over time dt taking into account
  !> a constant magnetic field using Boris method.
  subroutine PC_boris_advance(self, part, dt)
    use m_units_constants
    class(PC_t), intent(in)        :: self
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: dt
    real(dp)                       :: t_vec(3), tmp(3)

    ! Push the particle over dt/2
    part%x = part%x + 0.5_dp * dt * part%v ! Use the previous velocity
    part%v = part%v + 0.5_dp * dt * part%a

    ! Rotate the velocity
    tmp    = 0.5_dp * dt * self%B_vec * UC_elec_q_over_m
    t_vec  = 2 * tmp / (1.d0 + sum(tmp**2))
    tmp    = part%v + cross_product(part%v, tmp)
    tmp    = cross_product(tmp, t_vec)
    part%v = part%v + tmp

    ! Push the particle over dt/2
    part%v = part%v + 0.5_dp * dt * part%a
    part%x = part%x + 0.5_dp * dt * part%v ! Use the new velocity

    ! Update time left
    part%t_left = part%t_left - dt
  end subroutine PC_boris_advance

  !> Use a midpoint scheme to advance a particle tracing a velocity field
  subroutine PC_tracer_advance_midpoint(self, part, dt)
    class(PC_t), intent(in)        :: self
    type(PC_part_t), intent(inout) :: part
    real(dp), intent(in)           :: dt
    real(dp)                       :: x(3)

    x           = part%x
    part%x      = part%x + 0.5_dp * part%v * dt
    part%v      = self%tracer_velocity(part)
    part%x      = x + part%v * dt
    part%v      = self%tracer_velocity(part)
    part%t_left = part%t_left - dt
  end subroutine PC_tracer_advance_midpoint

  subroutine PC_after_dummy(self, dt)
    class(PC_t), intent(inout)     :: self
    real(dp), intent(in)           :: dt

    error stop "No procedure after_mover specified"
  end subroutine PC_after_dummy

  !> Return the cross product of vectors a and b
  pure function cross_product(a, b) result(vec)
    real(dp), intent(in) :: a(3)
    real(dp), intent(in) :: b(3)
    real(dp)             :: vec(3)

    vec(1) = a(2) * b(3) - a(3) * b(2)
    vec(2) = a(3) * b(1) - a(1) * b(3)
    vec(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  !> Perform the velocity correction of a Verlet scheme
  !!
  !! During the timestep x,v have been advanced to:
  !! x(t+1) = x(t) + v(t)*dt + 0.5*a(t)*dt^2,
  !! v(t+1) = v(t) + a(t)*dt
  !! But the velocity at t+1 should be v(t+1) = v(t) + 0.5*(a(t) + a(t+1))*dt
  !! to have a second order scheme, which is corrected here.
  subroutine PC_verlet_correct_accel(pc, dt)
    use m_units_constants
    class(PC_t), intent(inout) :: pc
    real(dp), intent(in)      :: dt
    integer                   :: ll
    real(dp)                  :: new_accel(3)

    if (.not. associated(pc%accel_function)) &
         stop "particle_core error: accel_func is not set"

    !$omp parallel do private(ll, new_accel)
    do ll = 1, pc%n_part
       new_accel = pc%accel_function(pc%particles(ll))
       pc%particles(ll)%v = pc%particles(ll)%v + &
            0.5_dp * (new_accel - pc%particles(ll)%a) * dt
       pc%particles(ll)%a = new_accel
    end do
    !$omp end parallel do
  end subroutine PC_verlet_correct_accel

  subroutine set_accel(self)
    class(PC_t), intent(inout) :: self
    integer                    :: ll

    !$omp parallel do
    do ll = 1, self%n_part
       self%particles(ll)%a = self%accel_function(self%particles(ll))
    end do
    !$omp end parallel do
  end subroutine set_accel

  !> Remove dead particles from the list
  subroutine clean_up(self)
    class(PC_t), intent(inout) :: self
    integer :: ix_end, ix_clean, n_part
    logical :: success

    do
       ! Get an index that has to be cleaned from the list
       call LL_pop(self%clean_list, ix_clean, success)
       if (.not. success) exit

       n_part      = self%n_part
       ! This is overridden if a replacement is found
       self%n_part = min(self%n_part, ix_clean-1)

       ! Find the last "alive" particle in the list
       do ix_end = n_part, ix_clean+1, -1
          if (self%particles(ix_end)%w /= PC_dead_weight) then
             self%particles(ix_clean) = self%particles(ix_end)
             self%n_part              = ix_end-1
             exit
          end if
       end do
    end do
  end subroutine clean_up

  !> Add a particle
  subroutine add_part(self, part)
    class(PC_t), intent(inout)  :: self
    type(PC_part_t), intent(in) :: part

    call self%check_space(self%n_part + 1)
    self%n_part                 = self%n_part + 1
    self%particles(self%n_part) = part
  end subroutine add_part

  !> Create a particles
  subroutine create_part(self, x, v, a, weight, t_left, id, ptype)
    class(PC_t), intent(inout) :: self
    real(dp), intent(IN)       :: x(3), v(3), a(3), weight, t_left
    type(PC_part_t)            :: my_part
    integer, intent(in), optional :: id, ptype
    my_part%x      = x
    my_part%v      = v
    my_part%a      = a
    my_part%w      = weight
    my_part%t_left = t_left

    if (present(id)) my_part%id = id
    if (present(ptype)) my_part%ptype = ptype

    call self%add_part(my_part)
  end subroutine create_part

  !> Mark particle for removal
  subroutine remove_part(self, ix_to_remove)
    class(PC_t), intent(inout) :: self
    integer, intent(in) :: ix_to_remove

    call LL_add(self%clean_list, ix_to_remove)
    self%particles(ix_to_remove)%w = PC_dead_weight
  end subroutine remove_part

  !> Map particle coordinates to a periodic domain
  subroutine periodify(self, is_periodic, lengths)
    class(PC_t), intent(inout) :: self
    logical, intent(in) :: is_periodic(3)
    real(dp), intent(in) :: lengths(3)
    integer :: i_dim, n_part

    n_part = self%n_part
    do i_dim = 1, 3
       if (is_periodic(i_dim)) then
          self%particles(1:n_part)%x(i_dim) = &
               modulo(self%particles(1:n_part)%x(i_dim), lengths(i_dim))
       end if
    end do
  end subroutine periodify

  !> Translate all particle coordinates
  subroutine translate(self, delta_x)
    class(PC_t), intent(inout) :: self
    real(dp), intent(in) :: delta_x(3)
    integer              :: ll
    do ll = 1, self%n_part
       self%particles(ll)%x = self%particles(ll)%x + delta_x
    end do
  end subroutine translate

  !> Return the number of real particles
  real(dp) function get_num_real_part(self)
    class(PC_t), intent(in) :: self

    get_num_real_part = sum(self%particles(1:self%n_part)%w)
  end function get_num_real_part

  !> Return the number of simulation particles
  integer function get_num_sim_part(self)
    class(PC_t), intent(in) :: self
    get_num_sim_part = self%n_part
  end function get_num_sim_part

  !> Loop over all the particles and call pptr for each of them
  subroutine loop_iopart(self, pptr)
    class(PC_t), intent(inout) :: self
    interface
       subroutine pptr(my_part)
         import
         type(PC_part_t), intent(inout) :: my_part
       end subroutine pptr
    end interface

    integer :: ll

    do ll = 1, self%n_part
       call pptr(self%particles(ll))
    end do
  end subroutine loop_iopart

  !> Loop over all the particles and call pptr for each of them
  subroutine loop_ipart(self, pptr)
    class(PC_t), intent(in) :: self
    interface
       subroutine pptr(my_part)
         import
         type(PC_part_t), intent(in) :: my_part
       end subroutine pptr
    end interface

    integer :: ll

    do ll = 1, self%n_part
       call pptr(self%particles(ll))
    end do
  end subroutine loop_ipart

  subroutine compute_vector_sum(self, pptr, my_sum)
    class(PC_t), intent(inout) :: self
    interface
       subroutine pptr(my_part, my_reals)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_reals(:)
       end subroutine pptr
    end interface
    real(dp), intent(out)      :: my_sum(:)

    integer                    :: ll
    real(dp)                   :: temp(size(my_sum))

    my_sum = 0.0_dp
    do ll = 1, self%n_part
       call pptr(self%particles(ll), temp)
       my_sum = my_sum + temp
    end do
  end subroutine compute_vector_sum

  subroutine compute_scalar_sum(self, pptr, my_sum)
    class(PC_t), intent(in) :: self
    interface
       subroutine pptr(my_part, my_real)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(out)       :: my_real
       end subroutine pptr
    end interface
    real(dp), intent(out)     :: my_sum

    integer                   :: ll
    real(dp)                  :: temp

    my_sum = 0.0_dp
    do ll = 1, self%n_part
       call pptr(self%particles(ll), temp)
       my_sum = my_sum + temp
    end do
  end subroutine compute_scalar_sum

  !> Get energy from velocity vector
  pure real(dp) function PC_v_to_en(v, mass)
    real(dp), intent(in) :: v(3), mass
    PC_v_to_en = 0.5_dp * mass * sum(v**2)
  end function PC_v_to_en

  !> Get energy from velocity (scalar)
  real(dp) elemental function PC_speed_to_en(vel, mass)
    real(dp), intent(in) :: vel, mass
    PC_speed_to_en = 0.5_dp * mass * vel**2
  end function PC_speed_to_en

  !> Get speed (scalar) from energy
  real(dp) elemental function PC_en_to_vel(en, mass)
    real(dp), intent(in) :: en, mass
    PC_en_to_vel = sqrt(2 * en / mass)
  end function PC_en_to_vel

  !> Check whether enough space is available
  subroutine check_space(self, n_req)
    class(PC_t), intent(in) :: self
    integer, intent(in) :: n_req

    if (n_req > size(self%particles)) then
       print *, "Error: particle list too small"
       print *, "Size of list:                 ", size(self%particles)
       print *, "Number of particles required: ", n_req
       stop
    end if
  end subroutine check_space

  !> Create a lookup table with cross collision rates from cross sections
  subroutine use_cross_secs(self, max_ev, table_size, cross_secs)
    use m_units_constants
    use m_cross_sec
    use m_lookup_table
    class(PC_t), intent(inout) :: self
    real(dp), intent(in)       :: max_ev
    integer, intent(in)        :: table_size
    type(CS_t), intent(in)     :: cross_secs(:)

    real(dp)                  :: vel_list(table_size), rate_list(table_size)
    real(dp)                  :: sum_rate_list(table_size)
    integer                   :: ix, i_c, i_row, n_colls
    real(dp)                  :: max_vel, en_eV

    max_vel      = PC_en_to_vel(max_ev * UC_elec_volt, self%mass)
    n_colls      = size(cross_secs)
    self%n_colls = n_colls
    allocate(self%colls(n_colls))
    allocate(self%coll_is_event(n_colls))
    self%coll_is_event(:) = .false.

    ! Set a range of velocities
    do ix = 1, table_size
       vel_list(ix) = (ix-1) * max_vel / (table_size-1)
       sum_rate_list(ix) = 0
    end do

    ! Create collision rate table
    self%rate_lt = LT_create(0.0_dp, max_vel, table_size, n_colls)

    do i_c = 1, n_colls
       self%colls(i_c) = cross_secs(i_c)%coll

       ! Linear interpolate cross sections by energy
       do i_row = 1, table_size
          en_eV = PC_speed_to_en(vel_list(i_row), self%mass) / UC_elec_volt
          call LT_lin_interp_list(cross_secs(i_c)%en_cs(1, :), &
               cross_secs(i_c)%en_cs(2, :), en_eV, rate_list(i_row))
          rate_list(i_row) = rate_list(i_row) * vel_list(i_row)
       end do

       ! We store the sum of the collision rates with index up to i_c, this
       ! makes the null-collision method straightforward.
       sum_rate_list = sum_rate_list + rate_list
       call LT_set_col(self%rate_lt, i_c, vel_list, sum_rate_list)
    end do

    self%max_rate = maxval(sum_rate_list)
    self%inv_max_rate = 1 / self%max_rate
  end subroutine use_cross_secs

  !> Create a lookup table with cross collision rates from a list of functions
  subroutine use_rate_funcs(self, max_ev, table_size, rate_funcs)
    use m_units_constants
    use m_cross_sec
    use m_lookup_table
    class(PC_t), intent(inout)    :: self
    real(dp), intent(in)          :: max_ev
    integer, intent(in)           :: table_size
    type(rate_func_t), intent(in) :: rate_funcs(:)

    real(dp)                  :: vel_list(table_size), rate_list(table_size)
    real(dp)                  :: sum_rate_list(table_size)
    integer                   :: ix, i_c, i_row, n_colls
    real(dp)                  :: max_vel

    max_vel      = PC_en_to_vel(max_ev * UC_elec_volt, self%mass)
    n_colls      = size(rate_funcs)
    self%n_colls = n_colls
    allocate(self%colls(n_colls))
    allocate(self%coll_is_event(n_colls))
    self%coll_is_event(:) = .false.

    ! Set a range of velocities
    do ix = 1, table_size
       vel_list(ix) = (ix-1) * max_vel / (table_size-1)
       sum_rate_list(ix) = 0
    end do

    ! Create collision rate table
    self%rate_lt = LT_create(0.0_dp, max_vel, table_size, n_colls)

    do i_c = 1, n_colls
       self%colls(i_c) = rate_funcs(i_c)%coll

       ! Linear interpolate cross sections by energy
       do i_row = 1, table_size
          rate_list(i_row) = rate_funcs(i_c)%ptr(vel_list(i_row))
       end do

       ! We store the sum of the collision rates with index up to i_c, this
       ! makes the null-collision method straightforward.
       sum_rate_list = sum_rate_list + rate_list
       call LT_set_col(self%rate_lt, i_c, vel_list, sum_rate_list)
    end do

    self%max_rate = maxval(sum_rate_list)
    self%inv_max_rate = 1 / self%max_rate
  end subroutine use_rate_funcs

  !> Sort the particles according to sort_func
  subroutine sort(self, sort_func)
    use m_mrgrnk
    class(PC_t), intent(inout)   :: self
    procedure(p_to_r_f)          :: sort_func

    integer                      :: ix, n_part
    integer, allocatable         :: sorted_ixs(:)
    real(dp), allocatable        :: part_values(:)
    type(PC_part_t), allocatable :: part_copies(:)

    n_part = self%n_part
    allocate(sorted_ixs(n_part))
    allocate(part_values(n_part))
    allocate(part_copies(n_part))

    do ix = 1, n_part
       part_copies(ix) = self%particles(ix)
       part_values(ix) = sort_func(part_copies(ix))
    end do

    call mrgrnk(part_values, sorted_ixs)

    do ix = 1, n_part
       self%particles(ix) = part_copies(sorted_ixs(ix))
    end do
  end subroutine sort

  subroutine histogram(self, hist_func, filter_f, &
       filter_args, x_values, y_values)
    use m_mrgrnk
    class(PC_t), intent(in) :: self
    procedure(p_to_r_f)   :: hist_func
    real(dp), intent(in)        :: x_values(:)
    real(dp), intent(out)       :: y_values(:)
    interface
       logical function filter_f(my_part, real_args)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(in) :: real_args(:)
       end function filter_f
    end interface
    real(dp), intent(in)        :: filter_args(:)

    integer                     :: ix, p_ix, o_ix, n_used, num_bins, n_part
    integer, allocatable        :: sorted_ixs(:)
    real(dp)                    :: boundary_value
    real(dp), allocatable       :: part_values(:)
    logical, allocatable        :: part_mask(:)

    n_part = self%n_part
    allocate(part_mask(n_part))
    do ix = 1, n_part
       part_mask(ix) = filter_f(self%particles(ix), filter_args)
    end do

    n_used = count(part_mask)

    allocate(part_values(n_used))
    allocate(sorted_ixs(n_used))

    p_ix = 0
    do ix = 1, self%n_part
       if (part_mask(ix)) then
          p_ix = p_ix + 1
          part_values(p_ix) = hist_func(self%particles(p_ix))
       end if
    end do

    call mrgrnk(part_values, sorted_ixs)

    num_bins = size(x_values)
    p_ix     = 1
    y_values = 0

    outer: do ix = 1, num_bins-1
       boundary_value = 0.5_dp * (x_values(ix) + x_values(ix+1))
       do
          if (p_ix == n_used + 1) exit outer
          o_ix = sorted_ixs(p_ix) ! Index in the 'old' list
          if (part_values(o_ix) > boundary_value) exit

          y_values(ix) = y_values(ix) + self%particles(o_ix)%w
          p_ix         = p_ix + 1
       end do
    end do outer

    ! Fill last bin
    y_values(num_bins) = sum(self%particles(sorted_ixs(p_ix:n_used))%w)
  end subroutine histogram

  ! Routine to merge and split particles. Input arguments are the coordinate
  ! weights, used to determine the 'distance' between particles. The first three
  ! elements of the array are the weights of the xyz position coordinates, the
  ! next three the weights of the xyz velocity coordinates. Max_distance is the
  ! maxium Euclidean distance between particles to be merged. The weight_func
  ! returns the desired weight for a particle, whereas the pptr_merge and
  ! pptr_split procedures merge and split particles.
  subroutine merge_and_split(self, x_mask, v_fac, use_v_norm, weight_func, &
       max_merge_distance, pptr_merge, pptr_split)
    use m_mrgrnk
    use kdtree2_module
    class(PC_t), intent(inout) :: self
    real(dp), intent(in)       :: v_fac, max_merge_distance
    logical, intent(in)        :: x_mask(3), use_v_norm
    procedure(sub_merge)       :: pptr_merge
    procedure(sub_split)       :: pptr_split

    procedure(p_to_r_f)    :: weight_func

    integer, parameter     :: num_neighbors  = 1
    integer, parameter     :: n_part_out_max = 2
    real(dp), parameter    :: large_ratio    = 1.5_dp
    real(dp), parameter    :: small_ratio    = 1 / large_ratio
    type(kdtree2), pointer :: kd_tree
    type(kdtree2_result)   :: kd_results(num_neighbors)

    integer               :: n_x_coord, n_coords
    integer               :: num_part, num_merge, num_split
    integer               :: p_min, p_max, n_too_far
    integer               :: o_ix, o_nn_ix
    integer               :: i, ix, neighbor_ix
    integer               :: n_part_out
    logical, allocatable  :: already_merged(:)
    integer, allocatable  :: sorted_ixs(:)
    real(dp), allocatable :: coord_data(:, :), weight_ratios(:)
    type(PC_part_t)       :: part_out(n_part_out_max)
    real(dp)              :: dist

    p_min                                    = 1
    p_max                                    = self%n_part
    num_part                                 = p_max - p_min + 1

    allocate(weight_ratios(num_part))
    allocate(sorted_ixs(num_part))

    do ix = 1, num_part
       weight_ratios(ix) = self%particles(p_min+ix-1)%w / &
            weight_func(self%particles(p_min+ix-1))
    end do

    num_merge = count(weight_ratios <= small_ratio)
    n_x_coord = count(x_mask)
    n_coords  = n_x_coord + 3
    if (use_v_norm) n_coords = n_coords - 2

    allocate(coord_data(n_coords, num_merge))
    allocate(already_merged(num_merge))
    already_merged = .false.
    n_too_far      = 0

    ! print *, "num_merge", num_merge
    ! Sort particles by their relative weight
    call mrgrnk(weight_ratios, sorted_ixs)

    ! Only create a k-d tree if there are enough particles to be merged
    if (num_merge > n_coords) then
       do ix = 1, num_merge
          o_ix = sorted_ixs(ix)
          coord_data(1:n_x_coord, ix) = pack(self%particles(o_ix)%x, x_mask)
          if (use_v_norm) then
             coord_data(n_x_coord+1, ix) = v_fac * norm2(self%particles(o_ix)%v)
          else
             coord_data(n_x_coord+1:, ix) = v_fac * self%particles(o_ix)%v
          end if
       end do

       ! Create k-d tree
       kd_tree => kdtree2_create(coord_data)

       ! Merge particles
       do ix = 1, num_merge
          if (already_merged(ix)) cycle

          call kdtree2_n_nearest_around_point(kd_tree, idxin=ix, &
               nn=num_neighbors, correltime=1, results=kd_results)
          neighbor_ix = kd_results(1)%idx

          if (already_merged(neighbor_ix)) cycle

          dist = norm2(coord_data(:, ix)-coord_data(:, neighbor_ix))
          if (dist > max_merge_distance) cycle

          ! Get indices in the original particle list
          o_ix = sorted_ixs(ix)
          o_nn_ix = sorted_ixs(neighbor_ix)

          ! Merge, then remove neighbor
          call pptr_merge(self%particles(o_ix), self%particles(o_nn_ix), &
               part_out(1), self%rng)
          self%particles(o_ix) = part_out(1)
          call self%remove_part(o_nn_ix)
          already_merged((/ix, neighbor_ix/)) = .true.
       end do

       call kdtree2_destroy(kd_tree)
    end if

    ! Split particles.
    num_split = count(weight_ratios >= large_ratio)
    ! print *, "num_split:", num_split

    do ix = num_part - num_split + 1, num_part
       o_ix = sorted_ixs(ix)
       call pptr_split(self%particles(o_ix), weight_ratios(o_ix), part_out, &
            n_part_out, self%rng)
       self%particles(o_ix) = part_out(1)
       do i = 2, n_part_out
          call self%add_part(part_out(i))
       end do
    end do

    ! print *, "clean up"
    call self%clean_up()
  end subroutine merge_and_split

  !> Should be thread safe, so it can be called in parallel as long as the
  !> ranges [i0:i1] don't overlap.
  subroutine merge_and_split_range(self, i0, i1, x_mask, v_fac, use_v_norm, &
       weight_func, max_merge_distance, pptr_merge, pptr_split)
    use m_mrgrnk
    use kdtree2_module
    class(PC_t), intent(inout) :: self
    integer, intent(in)        :: i0, i1
    real(dp), intent(in)       :: v_fac, max_merge_distance
    logical, intent(in)        :: x_mask(3), use_v_norm
    procedure(sub_merge)       :: pptr_merge
    procedure(sub_split)       :: pptr_split
    procedure(p_to_r_f)        :: weight_func

    integer, parameter     :: num_neighbors  = 1
    integer, parameter     :: n_part_out_max = 2
    real(dp), parameter    :: large_ratio    = 1.5_dp
    real(dp), parameter    :: small_ratio    = 1 / large_ratio
    integer, parameter     :: min_merge      = 20
    type(kdtree2), pointer :: kd_tree
    type(kdtree2_result)   :: kd_results(num_neighbors)

    integer                      :: n_x_coord, n_coords
    integer                      :: num_part, num_merge, num_split
    integer                      :: n_too_far, n_free, i_pbuf
    integer                      :: o_ix, o_nn_ix
    integer                      :: i, ix, neighbor_ix
    integer                      :: n_part_out
    logical, allocatable         :: already_merged(:)
    integer, allocatable         :: sorted_ixs(:), free_ixs(:)
    real(dp), allocatable        :: coord_data(:, :), weight_ratios(:)
    type(PC_part_t), allocatable :: p_buf(:)
    type(PC_part_t)              :: part_out(n_part_out_max)
    real(dp)                     :: dist

    num_part = i1 - i0 + 1
    allocate(weight_ratios(i0:i1))

    do ix = i0, i1
       weight_ratios(ix) = self%particles(ix)%w / &
            weight_func(self%particles(ix))
    end do

    num_merge = count(weight_ratios <= small_ratio)
    num_split = count(weight_ratios >= large_ratio)

    ! Exit if there is nothing to do
    if (num_merge < min_merge .and. num_split == 0) return

    ! Sort particles by their relative weight
    allocate(sorted_ixs(num_part))
    call mrgrnk(weight_ratios, sorted_ixs)
    sorted_ixs = sorted_ixs + i0 - 1

    n_free = 0                  ! Will be checked at the end

    ! Create a k-d tree if there are enough particles to be merged
    if (num_merge > min_merge) then
       n_x_coord         = count(x_mask)
       n_coords          = n_x_coord + 3
       if (use_v_norm) n_coords = n_coords - 2

       allocate(free_ixs(num_merge/2))
       allocate(coord_data(n_coords, num_merge))
       allocate(already_merged(num_merge))

       already_merged(:) = .false.
       n_too_far         = 0

       do ix = 1, num_merge
          o_ix = sorted_ixs(ix)
          coord_data(1:n_x_coord, ix) = pack(self%particles(o_ix)%x, x_mask)
          if (use_v_norm) then
             coord_data(n_x_coord+1, ix) = v_fac * norm2(self%particles(o_ix)%v)
          else
             coord_data(n_x_coord+1:, ix) = v_fac * self%particles(o_ix)%v
          end if
       end do

       ! Create k-d tree
       kd_tree => kdtree2_create(coord_data)

       ! Merge particles
       do ix = 1, num_merge
          if (already_merged(ix)) cycle

          call kdtree2_n_nearest_around_point(kd_tree, idxin=ix, &
               nn=num_neighbors, correltime=1, results=kd_results)
          neighbor_ix = kd_results(1)%idx

          if (already_merged(neighbor_ix)) cycle

          dist = norm2(coord_data(:, ix)-coord_data(:, neighbor_ix))
          if (dist > max_merge_distance) cycle

          ! Get indices in the original particle list
          o_ix = sorted_ixs(ix)
          o_nn_ix = sorted_ixs(neighbor_ix)

          ! Merge, then remove neighbor
          call pptr_merge(self%particles(o_ix), self%particles(o_nn_ix), &
               part_out(1), self%rng)
          self%particles(o_ix) = part_out(1)
          n_free               = n_free + 1
          free_ixs(n_free)     = o_nn_ix
          already_merged((/ix, neighbor_ix/)) = .true.
       end do

       call kdtree2_destroy(kd_tree)
    end if

    ! Split particles
    if (num_split > 0) then
       allocate(p_buf(num_split))
       i_pbuf            = 0

       do ix = num_part - num_split + 1, num_part
          o_ix = sorted_ixs(ix)
          call pptr_split(self%particles(o_ix), weight_ratios(o_ix), part_out, &
               n_part_out, self%rng)
          self%particles(o_ix) = part_out(1)

          if (n_free > 0) then
             i = free_ixs(n_free)
             n_free = n_free - 1
             self%particles(i) = part_out(2)
          else
             i_pbuf = i_pbuf + 1
             p_buf(i_pbuf) = part_out(2)
          end if
       end do

       if (i_pbuf > 0) then
          !$omp critical
          i = self%n_part
          self%n_part = self%n_part + i_pbuf
          !$omp end critical
          call self%check_space(i + i_pbuf)
          self%particles(i+1:i+i_pbuf) = p_buf(1:i_pbuf)
       end if
    end if

    if (n_free > 0) then
       !$omp critical
       do i = 1, n_free
          call self%remove_part(free_ixs(i))
       end do
       !$omp end critical
    end if

    ! Users have to call clean_up afterwards!
  end subroutine merge_and_split_range

  ! Merge two particles into part_a, should remove part_b afterwards
  subroutine PC_merge_part_rxv(part_a, part_b, part_out, rng)
    type(PC_part_t), intent(in)    :: part_a, part_b
    type(PC_part_t), intent(out) :: part_out
    type(RNG_t), intent(inout) :: rng

    if (rng%unif_01() > part_a%w / (part_a%w + part_b%w)) then
       part_out%x      = part_b%x
       part_out%v      = part_b%v
       part_out%a      = part_b%a
       part_out%t_left = part_b%t_left
    else
       part_out%x      = part_a%x
       part_out%v      = part_a%v
       part_out%a      = part_a%a
       part_out%t_left = part_a%t_left
    end if
    part_out%w = part_a%w + part_b%w
  end subroutine PC_merge_part_rxv

  subroutine PC_split_part(part_a, w_ratio, part_out, n_part_out, rng)
    type(PC_part_t), intent(in)    :: part_a
    real(dp), intent(in)           :: w_ratio
    type(PC_part_t), intent(inout) :: part_out(:)
    integer, intent(inout)         :: n_part_out
    type(RNG_t), intent(inout)     :: rng

    n_part_out    = 2
    part_out(1)   = part_a
    part_out(1)%w = 0.5_dp * part_out(1)%w
    part_out(2)   = part_out(1)
  end subroutine PC_split_part

  integer function get_num_colls(self)
    class(PC_t), intent(in) :: self
    get_num_colls = self%n_colls
  end function get_num_colls

  function get_mean_energy(self) result(mean_en)
    class(PC_t), intent(in) :: self
    real(dp)                :: mean_en, weight
    integer                 :: ll
    mean_en = 0
    weight  = 0
    do ll = 1, self%n_part
       weight  = weight + self%particles(ll)%w
       mean_en = mean_en + self%particles(ll)%w * &
            PC_v_to_en(self%particles(ll)%v, self%mass)
    end do
    mean_en = mean_en / weight
  end function get_mean_energy

  subroutine get_coll_rates(self, velocity, coll_rates)
    class(PC_t), intent(in) :: self
    real(dp), intent(in) :: velocity
    real(dp), intent(inout) :: coll_rates(:)
    integer :: i
    coll_rates = LT_get_mcol(self%rate_lt, velocity)

    do i = size(coll_rates), 2, -1
       coll_rates(i) = coll_rates(i) - coll_rates(i-1)
    end do
  end subroutine get_coll_rates

  subroutine get_colls(self, out_colls)
    class(PC_t), intent(in) :: self
    type(CS_coll_t), intent(inout), allocatable :: out_colls(:)
    allocate(out_colls(self%n_colls))
    out_colls = self%colls
  end subroutine get_colls

  subroutine get_coeffs(self, coeff_data, coeff_names, n_coeffs)
    use m_cross_sec
    use m_lookup_table
    class(PC_t), intent(in)                    :: self
    real(dp), intent(out), allocatable         :: coeff_data(:,:)
    character(len=*), intent(out), allocatable :: coeff_names(:)
    integer, intent(out)                       :: n_coeffs
    type(CS_coll_t), allocatable               :: coll_data(:)
    integer                                    :: nn, n_rows

    call self%get_colls(coll_data)
    n_coeffs = self%n_colls + 1
    n_rows = self%rate_lt%n_points
    allocate(coeff_data(n_coeffs, n_rows))
    allocate(coeff_names(n_coeffs))

    call LT_get_data(self%rate_lt, coeff_data(1, :), coeff_data(2:,:))
    coeff_names(1) = "velocity (m/s)"
    do nn = 1, self%n_colls
       select case (self%colls(nn)%type)
       case (CS_ionize_t)
          coeff_names(1+nn) = "ionization"
       case (CS_attach_t)
          coeff_names(1+nn) = "attachment"
       case (CS_elastic_t)
          coeff_names(1+nn) = "elastic"
       case (CS_excite_t)
          coeff_names(1+nn) = "excitation"
       case default
          coeff_names(1+nn) = "unknown"
       end select
    end do
  end subroutine get_coeffs

end module m_particle_core
