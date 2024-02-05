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
         call obs_op_gridpoint(i_obs, dim_p, dim_obs, state_p, ostate)
      end do

      call system_clock(timer_op_end, t_rate)
      op_dur = op_dur + &
          (real(timer_op_end, wp) - real(timer_op_start, wp))/real(t_rate, wp)
   end subroutine obs_op_pdafomi
end module mod_u_pdafomi_pdaf

