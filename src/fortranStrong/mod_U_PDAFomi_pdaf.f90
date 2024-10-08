! this file is part of pypdaf

! copyright (c) 2022 university of reading and
! national centre for earth observation

! this program is free software: you can redistribute it and/or modify
! it under the terms of the gnu general public license as published by
! the free software foundation, either version 3 of the license, or
! (at your option) any later version.

! this program is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the
! gnu general public license for more details.

! you should have received a copy of the gnu general public license
! along with this program.  if not, see <http://www.gnu.org/licenses/>.

module mod_u_pdafomi_pdaf
use mod_kind_pdaf, only: wp
use mod_observations_pdaf, only: n_obs, obs
implicit none

integer :: timer_getobs_start, timer_getobs_end, t_rate
integer :: timer_dimomi_start, timer_dimomi_end
integer :: timer_op_start, timer_op_end
real(wp) :: getobs_dur, dimomi_dur, op_dur

contains
   ! write the synthetic observations into netcdf file
   subroutine get_obs_f(step, dim_obs_f, observation_f)
      use mod_obswriter_pdaf, only: writeobs
      ! arguments
      integer,  intent(in)  :: step                     ! current time step
      integer,  intent(in)  :: dim_obs_f                ! dimension of obs. vector
      real(wp), intent(in)  :: observation_f(dim_obs_f) ! observation vector
      ! local variables
      integer :: i_obs, istart, iend

      call system_clock(timer_getobs_start)

      istart = 1
      do i_obs = 1, n_obs
         iend = istart + obs(i_obs)%dim_obs - 1
         call writeobs(i_obs, step, observation_f(istart:iend))
         istart = iend + 1
      end do
      call system_clock(timer_getobs_end, t_rate)
      getobs_dur = getobs_dur + &
          (real(timer_getobs_end, wp) - real(timer_getobs_start, wp))/real(t_rate, wp)
   end subroutine get_obs_f

   !! set up the observation information and obtain the
   !! size of the observation vector for synthetic obs.
   subroutine init_dim_obs_gen_pdafomi(step, dim_obs)
      use mod_observations_pdaf, only: init_dim_obs_gen
      ! arguments
      integer, intent(in)  :: step     !< current time step
      integer, intent(out) :: dim_obs  !< dimension of full observation vector
      ! local variables
      integer :: this_dim_obs
      integer :: i_obs

      call system_clock(timer_dimomi_start)

      dim_obs = 0
      do i_obs = 1, n_obs
         call init_dim_obs_gen(i_obs, this_dim_obs)
         dim_obs = dim_obs + this_dim_obs
      end do

      call system_clock(timer_dimomi_end, t_rate)
      dimomi_dur = dimomi_dur + &
          (real(timer_dimomi_end, wp) - real(timer_dimomi_start, wp))/real(t_rate, wp)
   end subroutine init_dim_obs_gen_pdafomi

   !! set up the observation information and obtain the
   !! size of the observation vector for assimilation
   subroutine init_dim_obs_pdafomi(step, dim_obs)
      use mod_observations_pdaf, only: init_dim_obs

      integer, intent(in)  :: step     !< current time step
      integer, intent(out) :: dim_obs  !< dimension of full observation vector

      integer :: this_dim_obs
      integer :: i_obs

      call system_clock(timer_dimomi_start)

      dim_obs = 0
      do i_obs = 1, n_obs
         if (obs(i_obs)%doassim == 0) cycle
         call init_dim_obs(i_obs, this_dim_obs)
         dim_obs = dim_obs + this_dim_obs
      end do

      call system_clock(timer_dimomi_end, t_rate)
      dimomi_dur = dimomi_dur + &
          (real(timer_dimomi_end, wp) - real(timer_dimomi_start, wp))/real(t_rate, wp)
   end subroutine init_dim_obs_pdafomi

   ! observation operator for obs. on gridpoints
   subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)
      ! include functions for different observations
      use mod_observations_pdaf, only: obs_op_gridpoint
      implicit none
      ! *** arguments ***
      integer,  intent(in)    :: step                 !< current time step
      integer,  intent(in)    :: dim_p                !< pe-local state dimension
      integer,  intent(in)    :: dim_obs              !< dimension of full observed state
      real(wp), intent(in)    :: state_p(dim_p)       !< pe-local model state
      real(wp), intent(inout) :: ostate(dim_obs)      !< pe-local full observed state
      ! local variables
      integer :: i_obs

      ! ******************************************************
      ! *** apply observation operator h on a state vector ***
      ! ******************************************************
      call system_clock(timer_op_start)

      do i_obs = 1, n_obs
         if (obs(i_obs)%doassim == 0) cycle
         call obs_op_gridpoint(i_obs, dim_p, dim_obs, state_p, ostate)
      end do

      call system_clock(timer_op_end, t_rate)
      op_dur = op_dur + &
          (real(timer_op_end, wp) - real(timer_op_start, wp))/real(t_rate, wp)
   end subroutine obs_op_pdafomi

   SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)
      use mod_observations_pdaf, only: init_dim_obs_l
      use mod_model_pdaf, only: nx, ny, pi, maooam_model

       ! *** Arguments ***
      integer, intent(in)  :: domain_p     !< Index of current local analysis domain
      integer, intent(in)  :: step         !< Current time step
      integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
      integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector

      real(wp) :: coords_l(2)
      real(wp) :: dx, dy, n
      integer :: state_index
      integer :: i_obs
      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************
      state_index = domain_p - ((domain_p - 1)/nx/ny)*nx*ny
      n = maooam_model%model_configuration%physics%n
      dx = 2*pi/n/(nx - 1)
      dy = pi/(ny - 1)
      coords_l(1) = ceiling(real(state_index)/real(ny))
      coords_l(2) = state_index - (coords_l(1) - 1)*ny
      coords_l(1) = (coords_l(1) - 1)*dx
      coords_l(2) = (coords_l(2) - 1)*dy

      do i_obs = 1, n_obs
         if (obs(i_obs)%doassim == 0) cycle
         call init_dim_obs_l(i_obs, coords_l, dim_obs_l)
      end do
      
   END SUBROUTINE init_dim_obs_l_pdafomi
end module mod_u_pdafomi_pdaf

